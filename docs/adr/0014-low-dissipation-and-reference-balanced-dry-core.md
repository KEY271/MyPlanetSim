# ADR 0014: low-dissipation transport and reference-balanced terrain pressure force

**Status:** accepted

## Context

The dry hydrostatic flux used the horizontal Lamb speed
`sqrt(gamma*Rd*T)` as a scalar Rusanov coefficient on mass, both momentum
components, potential temperature, and tracer. At the Phase 7 `N=12` resolution the
piecewise-constant equivalent coefficient was `1.259e8 m2/s`, while the configured
explicit Laplacian coefficient was `1e5 m2/s`.

The DCMIP 2-0-0 rest state also evaluated `grad(Phi)` and
`alpha*grad(p)` separately. Its registered initial acceleration did not decrease over
`N=8/16/32`, although extended refinement showed that this was a pre-asymptotic terrain
resolution effect rather than a non-convergent gradient: the L2 residual decreased from
`1.221e-4` at `N=48` to `3.728e-5 m/s2` at `N=96`.

## Decision

The dry horizontal physical flux is an advective block because pressure is applied by a
separate source. Its jump dissipation therefore uses

```text
a_adv = max(abs(u_n,L), abs(u_n,R))
```

The fast estimate

```text
a_fast = max(abs(u_n,L) + sqrt(gamma*Rd*T_L),
             abs(u_n,R) + sqrt(gamma*Rd*T_R))
```

is retained unchanged for the per-cell SSPRK3 CFL. This is a split by equation block,
not a constant reduction of the wave speed.

The named DCMIP hydrostatic rest benchmark constructs one immutable reference field at
driver construction. Its pressure force is evaluated as

```text
G = D(Phi - Phi_ref)
  + alpha * D(p - p_ref)
  + (alpha - alpha_ref) * D(p_ref)
```

so the reference state has an exactly zero semi-discrete pressure force. Reference-state
balancing is not silently applied to jets or forced cases. Restart reconstructs the same
reference from configuration; it is not prognostic or serialized.

## Validation

For a smooth `40*cos(latitude)^3` jet with linear Barth--Jespersen reconstruction, the
Laplacian-equivalent kinetic-energy dissipation changed as follows:

| N | before (m2/s) | after (m2/s) |
|---:|---:|---:|
| 8 | 1.124e7 | 1.727e5 |
| 12 | 2.935e6 | 4.166e4 |
| 16 | 1.158e6 | 1.592e4 |
| 24 | 3.205e5 | 4.292e3 |
| 32 | 1.306e5 | 1.731e3 |
| 48 | 3.733e4 | 4.913e2 |
| 64 | 1.546e4 | 2.030e2 |

The corresponding tracer values are of the same order. The `N=12` jet reduction is a
factor of 70.5, making the `1e5 m2/s` explicit Laplacian larger than the diagnosed smooth
numerical viscosity. Mean fast speed remains `338.9 m/s` in this test.

DCMIP initial acceleration is exactly zero at `N=8/16/32`. A six-day `N=8`, `K=15`
integration completed 8,640 steps with maximum wind `7.734e-11 m/s`; the maximum raw
prognostic-state change was `3.866e-8` in its stored units. The flat linear-wave preset
completed 600 s at `N=4/8/16/32`; its surface-pressure mode converged to `0.999259 Pa`
while retaining the acoustic CFL.

The release suite passes all 75 tests, including RHS zero-allocation and checkpoint
restart gates.

## Consequences

- Explicit diffusion now controls smooth resolved jet and tracer damping at the coarse
  production resolution used by the Phase 7 candidate.
- A stationary terrain-following reference has neither Rusanov mass/thermal diffusion
  nor pressure-gradient acceleration.
- The Lamb mode still constrains the explicit time step; this ADR changes dissipation,
  not the physical fast-wave estimate.
- General local hydrostatic reconstruction and a vertically coupled characteristic
  external-mode flux remain possible later extensions.
