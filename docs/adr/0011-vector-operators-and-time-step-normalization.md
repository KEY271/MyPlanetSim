# ADR 0011: Spherical vector operators and per-cell time-step normalization

**Status:** accepted for bounded Phase 9 implementation

## Context

A full-code review after Phase 8 (`c3246a9`) found two defects that share one property:
both are silent. Neither raises, neither fails an existing test, and both change results
of any run that exercises them.

The finite-volume vector Laplacian on the sphere is assembled from a divergence gradient
and a vorticity gradient. Its rotational term carries the wrong sign, and the curvature
term of the spherical vector Laplacian is absent. Measured on a unit sphere with a rigid
rotation `u = z_hat cross r`, the implementation returns `<lap(u).u>/<u.u> = +1.996`
where the exact value is `-1`. The defect is hidden because the only registered consumer
is the shallow-water biharmonic path, where the sign error is applied twice and cancels,
and because the sole operator test uses a zero field.

The dry hydrostatic core limits its horizontal time step per edge,

```text
dt = min over faces of ( cfl * min(A_left, A_right) / (L_f * lambda_f) )
```

while `spherical_transport.cpp` and `shallow_water_driver.cpp` limit it per cell,

```text
dt * sum over faces of ( lambda_f * L_f ) / A_cell <= cfl
```

On a quadrilateral cell the edge form is about four times weaker. The registered Phase 7
production candidate reports CFL 0.45 while running at an effective cell CFL near 1.02.
The same configuration key `cfl` therefore means two different things depending on which
solver reads it.

The Phase 8 surface reservoir integrates `C dT_s/dt = ...` explicitly with no stability
constraint at all. Only positivity is checked, so a small `C_surface` may oscillate and
still pass every gate.

## Decision

### Spherical vector Laplacian

For a tangent field `v` on a sphere of radius `a`, with `D = div(v)`,
`zeta = k_hat . curl(v)`, and `k_hat` the outward unit radial direction:

```text
lap(v) = grad(D) + k_hat cross grad(zeta) + v / a^2
```

The second term takes the opposite sign from the current implementation. The third term
is the curvature correction and is currently absent. It is `+v/a^2` for both the
rotational and the divergent part, so for a field built from a degree-`n` spherical
harmonic,

```text
v = k_hat cross grad(psi_n)  ->  lap(v) = (1 - n(n+1)) / a^2 * v
v = grad(chi_n)              ->  lap(v) = (1 - n(n+1)) / a^2 * v
```

and rigid rotation (`n = 1`) gives `lap(v) = -v/a^2`.

`kBiharmonic` currently produces correct decay because the sign error cancels; it must
keep producing decay after the fix. `kLaplacian` momentum diffusion must reduce the
kinetic energy of a smooth rotational field, which it presently does not.

### One CFL definition for three solvers

`transport.cfl`, `shallow_water.cfl`, and `dry_hydrostatic.cfl` denote the same quantity:
the per-cell sum of face wave speeds times face lengths, divided by cell area, times the
time step. The dry core is brought to that form. With
`lambda = |u_n| + sqrt(gamma * Rd * T)` and `gamma = cp / (cp - Rd)`,

```text
dt <= cfl * A_cell / sum over faces of ( lambda_f * L_f )
```

evaluated on both cells of each face. A configured `cfl` value is portable between the
three solvers from this point on.

### Surface reservoir time step

The explicit surface integration is constrained before the fact rather than retried after
a positivity violation:

```text
tau_surface(c) = C_surface[c] / (4 * epsilon * sigma_SB * T_s[c]^3)
dt <= surface.cfl * min over cells of tau_surface(c)
```

`surface.cfl` defaults to 0.5. The existing positivity retry stays as a backstop but is
no longer the only guard.

### Preset consequences

The dry-core change tightens the stable step for presets whose configured
`run.time_step_s` exceeds the per-cell limit. Measured effect:

| preset | `N`/`K` | configured step | per-cell limit | affected |
|---|---|---:|---:|---|
| `phase5_*` (8 presets) | 4 / 8 | 10 s | ~780 s | no |
| `phase6_dcmip_2_0_0` | 8 / 15 | 60 s | re-evaluated after ADR 0013 | — |
| `phase6_jw06_steady`, `phase6_jw06_baroclinic` | 16 / 8 | 300 s | ~188 s | **yes** |
| `phase6_linear_mountain_wave` | 16 / 8 | 30 s | ~204 s | no |
| `phase7_held_suarez_short` | 2 / 2 | 10 s | ~1700 s | no |
| `phase7_held_suarez`, `phase8_earth_like`, `phase8_rapid_rotator`, `phase8_slow_rotator`, `phase8_tidally_locked` | 12 / 20 | 600 s | ~264 s | **yes** |
| `phase8_earth_geography` | 4 / 4 | 60 s | ~840 s | no |
| all shallow-water and transport presets | — | — | — | no (already per cell) |

None of the seven affected presets has a registered production run, so no existing
baseline is invalidated and every CI gate stays byte-exact. Doing this after the Phase 7
production matrix has run would instead cost six re-runs.

## Consequences

- Explicit horizontal diffusion becomes usable for the first time (ADR 0006, wired in
  Phase 9), because it would otherwise be wired to a wrong operator.
- The Phase 7 production candidate's step drops from 600 s to roughly 264 s, which is the
  cost driver recorded in the ADR 0009 amendment.
- `cfl` is comparable across solvers and across presets, so a stability question has one
  answer instead of three.
- The surface reservoir cannot silently oscillate under a small heat capacity.

## Rejected alternatives

- **Keep the edge-wise dry-core limit and rename the key.** Two names for the same
  physical condition preserves the trap rather than removing it, and leaves the three
  solvers incomparable.
- **Fix only the Laplacian sign and skip the curvature term.** The operator would still
  disagree with the analytic spherical-harmonic eigenvalues, so no exact gate could be
  registered.
- **Change the shallow-water and transport solvers to the edge-wise form instead.** The
  edge form is not the finite-volume Courant condition; it would loosen two solvers that
  are currently correct.
- **Clip or retry the surface temperature.** Retrying after a violation hides oscillation
  that stays inside the positive range.

## Validation obligations

- Rigid rotation `u = Omega cross r` reproduces `-u/a^2` within a registered tolerance.
- Rotational and divergent analytic fields at `n = 1, 2` have observed order at least 1.
- A constant vector field leaves only the curvature term.
- `kLaplacian` momentum diffusion decreases kinetic energy; `kBiharmonic` still decays.
- The dry-core stable step matches `cfl * A / sum(lambda*L)` on an analytic cell, and the
  three solvers agree on one CFL definition.
- The surface step honours `surface.cfl * min tau_surface`.

Commit order and preset re-registration are recorded in
[the Phase 9 implementation plan](../phase-9-plan.md).
