# ADR 0005: hybrid vertical coordinate and column state

**Status:** accepted for Phase 4

## Context

Phase 4 adds a vertical one-dimensional column before the cubed-sphere horizontal
dynamics and the vertical operators are combined in Phase 5. The column must define
one pressure, staggering, mass, thermodynamic, and flux convention that can later be
used independently at every horizontal cell.

The vertical coordinate moves when surface pressure changes. Treating that movement
as an occasional interpolation would make constant-state preservation and the mass
budget depend on a separate remapper. It would also hide the relation between
horizontal mass convergence, surface-pressure tendency, and vertical mass flux that
the three-dimensional model needs.

## Decision

Use a fixed-in-time hybrid sigma-pressure coordinate with downward-increasing index.
For `nz` layers, interface pressure is

```text
p_half[k] = A_pa[k] + B[k] * ps,       k = 0 .. nz
```

Production atmospheric columns require

```text
A_pa[0] = p_top,  B[0] = 0,
A_pa[nz] = 0,     B[nz] = 1,
p_top > 0.
```

The configured `A/B` arrays are immutable during a run. Their lengths, finiteness,
endpoints, and pressure monotonicity are validated before a state is created. Because
each interface pressure is affine in `ps`, monotonicity and the minimum pressure
thickness over a configured surface-pressure interval are checked at both interval
endpoints.

Layer pressure thickness and air mass per unit horizontal area are diagnostic:

```text
delta_p[k] = p_half[k+1] - p_half[k]
mass[k] = delta_p[k] / g
```

The prognostic column state stores `ps`, mass-weighted potential temperature, and one
mass-weighted passive tracer. It does not independently store `delta_p` or layer air
mass, so coordinate and state mass cannot silently diverge. The checkpoint layout is
versioned as `vertical_column_hybrid_v1`.

Potential temperature uses the dimensionless Exner function

```text
Pi(p) = (p / p0)^kappa,       kappa = Rd / cp
T = theta * Pi.
```

Full-level Exner pressure is the log-pressure mean

```text
Pi_full[k] = (Pi_lower - Pi_upper)
             / (kappa * log(p_lower / p_upper)),
p_full[k] = p0 * Pi_full[k]^(1/kappa).
```

A numerically stable limiting branch is used for nearly equal interface pressures.
This choice places the full level strictly inside its layer and makes the hydrostatic
layer integral

```text
Phi_upper = Phi_lower + cp * theta[k] * (Pi_lower - Pi_upper)
```

exact for layer-constant potential temperature. It is also exact for an isothermal
column when the full-level potential temperature is sampled using `Pi_full`.
Interface geopotential is integrated upward from the configured surface geopotential.

For the future three-dimensional continuity equation, let `H[k]` be the horizontal
mass tendency of a layer in `kg m^-2 s^-1`, positive for net horizontal convergence.
With impermeable top and surface boundaries,

```text
ps_dot = g * sum(H[k])
mass_dot[k] = (B[k+1] - B[k]) * ps_dot / g
F[0] = 0
F[k+1] = F[k] + H[k] - mass_dot[k]
```

where `F` is coordinate-relative vertical mass flux, positive downward. The lower
boundary residual `F[nz]` must be zero to rounding tolerance; it is reported and never
silently corrected. Mass-weighted potential temperature and tracers use the same
interface mass flux. A constant scalar is preserved when its supplied horizontal
mass-weighted tendency is the scalar value times `H`.

Phase 4 keeps `A/B` fixed and therefore does not use conservative remapping. If a later
phase introduces time-varying coefficients, vertical grid regeneration, or model-level
adaptation, that feature requires a separate ADR and a conservative remap with its own
mass, boundedness, and constant-preservation gates.

Phase 4 also does not add a column-only Web protocol or visualizer. The user-facing
visualizer is updated in Phase 5 after the column is coupled to the cubed-sphere and the
three-dimensional state, staggering, and frame-size requirements are known.

## Limits fixed for Phase 4

- The column has no horizontal geometry and no horizontal momentum.
- Top and surface mass fluxes are zero in production column runs.
- Dry potential temperature is advected adiabatically; pressure-work and total-energy
  closure belong to the coupled primitive equations in Phase 5.
- One passive tracer is sufficient to validate consistency and boundedness. A
  multi-tracer container may be added later without changing the coordinate contract.
- `p_top = 0` is rejected because Exner/log-pressure formulas would be singular. A
  sigma-like limit uses a positive model-top pressure.

## Consequences

The column mass is determined by one scalar `ps`, and the layer sum is structurally
`(ps - p_top) / g`. Surface-pressure changes, vertical mass flux, and constant-state
preservation can be tested without an interpolation path. Phase 5 must supply the
horizontal mass and scalar tendencies consumed by the flux diagnosis.

The pressure-coordinate limit with `B=0` remains an operator fixture rather than a
production atmospheric column because its bottom pressure is independent of `ps`.
The production sigma-like fixture uses
`A[k] = p_top * (1 - B[k])`, so both required endpoints remain valid.
