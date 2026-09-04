# ADR 0006: dry hydrostatic state and horizontal-vertical coupling

**Status:** accepted for Phase 5

## Context

Phase 5 combines the Phase 2 cubed-sphere dynamics with the Phase 4 hybrid-pressure
column. The coupling must preserve the Phase 4 relation between horizontal layer-mass
convergence, surface-pressure tendency, and coordinate-relative vertical mass flux. It
must also keep Cartesian tangent vectors at cubed-sphere seams and provide one state from
which checkpoints, diagnostics, and later visualization are derived.

The first three-dimensional core is a reference implementation. Its purpose is to expose
coupling and budget errors before adding terrain, physics, implicit solvers, multiple
tracers, or parallel storage. Phase 2 selected the collocated Rusanov scheme over the
compatible candidate on the current non-orthogonal grid. Phase 5 does not reopen that
decision without one of the re-evaluation conditions in ADR 0003.

## Decision

### Prognostic state

For `C = 6*N*N` horizontal cells and `K` layers, store

```text
DryHydrostaticStateV1
  time_s
  step
  surface_pressure_pa[C]
  horizontal_momentum_mass[C*K]       # M*u, Cartesian tangent vector
  potential_temperature_mass[C*K]     # M*theta
  tracer_mass[C*K]                    # M*q
```

where

```text
p_half[c,k] = A[k] + B[k] * ps[c]
M[c,k]      = (p_half[c,k+1] - p_half[c,k]) / g
offset(c,k) = c*K + k
```

`M`, pressure, Exner function, temperature, density, geopotential, velocity, potential
temperature, and tracer mixing ratio are diagnostic. Layer air mass is not stored a
second time. The state owns one passive tracer because that is sufficient to test
mass consistency; a general tracer registry is not part of Phase 5.

Momentum uses the planet-fixed Cartesian representation from ADR 0002. At every stage
boundary it is projected onto the tangent plane at the horizontal cell centre. A state is
invalid if any value is non-finite, pressure is outside the configured valid range, a
layer is not positive, temperature is below its floor, a nonnegative tracer is negative,
or radial momentum exceeds the stated roundoff tolerance. Invalid states stop the run;
they are not repaired by clipping.

The checkpoint layout is `dry_hydrostatic_cell_column_v1`. Its field order and
`cell-major, then downward-increasing k` indexing are fixed independently of any browser
frame layout.

### Layer equations and flux signs

Let `H[c,k]` be horizontal air-mass tendency per unit horizontal area, positive for net
convergence, as in ADR 0005. At each Runge--Kutta stage and horizontal cell, diagnose

```text
ps_dot = g * sum_k(H[k])
M_dot_target[k] = (B[k+1] - B[k]) * ps_dot / g
F[0] = 0
F[k+1] = F[k] + H[k] - M_dot_target[k]
```

`F` is coordinate-relative vertical mass flux in `kg m^-2 s^-1`, positive downward.
`F[K]` must close to zero without correction. For any mass-specific quantity `chi`,

```text
d(M*chi)[k]/dt = H_chi[k] + F_chi[k] - F_chi[k+1]
```

and the same diagnosed `F` and vertical reconstruction policy are used for momentum,
potential temperature, and tracer. The Phase 4 donor-cell and bounded linear operators
are reused component-wise; no independent three-dimensional remapper is introduced.
The corresponding coordinate-relative pressure velocity is `omega_relative = g*F` in
`Pa s^-1`. It is not labelled as geometric vertical velocity `w`.

### Horizontal reference discretization

Each model level uses the current unique-edge cubed-sphere traversal, linear
reconstruction, and Barth--Jespersen limiter. Face primitives are `M`, tangent `u`,
`theta`, and `q`. For an edge-normal velocity `u_n`, the advective physical fluxes are

```text
f_M     = M*u_n
f_P     = M*u_n*u
f_theta = M*u_n*theta
f_q     = M*u_n*q
```

The Rusanov flux uses one edge wave-speed estimate

```text
a = max(|u_n,L| + sqrt(gamma*Rd*T_L),
        |u_n,R| + sqrt(gamma*Rd*T_R))
gamma = cp/(cp-Rd)
```

and is evaluated once per unique edge and layer. This is the explicit reference estimate
for the hydrostatic external mode, not a claim that acoustic modes are being integrated.
The scalar flux is constructed from the same left/right mass states and the same `a`, so
constant `theta` or `q` gives exactly that constant times the air-mass flux. Vector face
values are represented in the shared edge tangent plane before the left/right jump is
formed.

The horizontal pressure-gradient acceleration at full level is defined by the generalized
pressure-coordinate identity

```text
G[k] = grad_eta(Phi_full[k]) + specific_volume[k] * grad_eta(p_full[k])
specific_volume = Rd*T_full/p_full
```

using the existing least-squares spherical gradient for both terms and the Phase 4
log-pressure full level and hydrostatic integral. Momentum receives `-M*G`. Using the same
gradient operator for `Phi` and `p` makes the cancellation testable; Phase 5 does not claim
exact discrete energy conservation from this reference pressure-gradient form.

Coriolis is

```text
-2 * M * P_tangent(Omega cross u)
```

and is retained separately from advection, pressure gradient, vertical transport, and
optional diffusion in diagnostics. Existing horizontal Laplacian/biharmonic machinery
may be applied level by level by a benchmark preset, but Phase 5 adds no vertical mixing
or new filter family.

### Time integration

Use the existing SSP-RK3 method of lines. Every stage performs this complete sequence:

```text
state -> hybrid/thermodynamic/hydrostatic diagnostics
      -> unique-edge horizontal tendencies
      -> per-column ps_dot and vertical mass flux
      -> pressure-gradient/Coriolis/diffusion sources
      -> shared vertical transport
      -> stage update, tangent projection, invariant validation
```

One step size satisfies the horizontal Rusanov CFL, Phase 4 vertical transport CFL,
surface-pressure bounds, positivity conditions, and the requested maximum. The final
step lands exactly on the requested end time. Split-explicit, semi-implicit, IMEX, and
subcycled barotropic solvers require a later ADR if measurements show that the explicit
reference method cannot meet the Phase 5 validation gate.

### Lower and upper boundaries

Phase 5 has a spherical, flat, impermeable lower boundary with constant
`surface_geopotential = 0`. The model top is the fixed positive pressure from ADR 0005.
Top and bottom coordinate-relative mass fluxes are zero. Horizontally varying surface
geopotential, terrain-following pressure-gradient corrections, mountain torque, and
surface exchange belong to Phase 6.

### Diagnostics and claims

Dry mass, potential-temperature mass, and tracer mass must close to rounding in unforced,
impermeable runs. Total dry energy is diagnosed as

```text
sum_cell area * sum_k M * (0.5*|u|^2 + cv*T + Phi_full)
```

and absolute axial angular momentum includes both relative and planetary components.
Their component tendencies and end-to-end residuals are reported, but exact conservation
is not claimed for the Rusanov/least-squares reference scheme. Acceptance is based on
bounded drift, convergence under refinement, and attribution of dissipation rather than
on hiding it in a residual.

## Visualizer boundary

The C++ state, checkpoint, benchmark, and numerical gates are completed before the Web
protocol changes. The browser never reconstructs hybrid pressure, integrates hydrostatic
geopotential, or advances a column. A later `FrameV2` is a derived, read-only view of a
validated C++ state and does not define the solver layout.

The first Phase 5 UI supports fixed three-dimensional presets with no new initial-condition
edit type. Temperature painting, tracer painting, wind impulses, arbitrary vertical
subsetting, remote execution, and a general analysis dashboard are separate features.

## Consequences

- Phase 4's continuity recurrence is the only horizontal/vertical mass-coupling path.
- Horizontal and vertical scalar transport share air-mass fluxes, making constant-field
  and global-budget tests structural.
- The collocated reference remains diffusive and the explicit external-mode CFL can be
  restrictive, but both are measurable without adding a second solver in the same phase.
- A no-topography baroclinic benchmark must use a genuinely flat-surface formulation.
  Phase 5 uses the Ullrich--Melvin--Jablonowski--Staniforth 2014 case; the smooth-topography
  Jablonowski--Williamson 2006 variant remains in Phase 6.

## Implementation audit

The Phase 7 entry gate re-read this decision against `DryHydrostaticDriver`. The first
implementation took one forward Euler step per interval, limited it by the horizontal
Rusanov CFL only, and used piecewise-constant face values. All three now follow the text
above: the complete stage sequence is re-evaluated at each of the three SSP-RK3 stages, the
step size additionally honours the Phase 4 vertical transport CFL and the surface-pressure
bounds, and each level is reconstructed linearly with the Barth--Jespersen limiter. The
tests and the registered unforced flat baseline are recorded in
[Phase 5 validation](../validation/phase-5.md). The decision itself is unchanged.
