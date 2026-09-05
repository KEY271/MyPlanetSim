# ADR 0013: Balanced benchmark initial states

**Status:** accepted for bounded Phase 9 implementation

## Context

Three presets claim a property they do not have.

`umjs14_steady` initializes a constant isothermal temperature field, a solid-body wind
`u = 35 cos(phi)`, a uniform surface pressure, and flat terrain. That is neither the
Ullrich--Melvin--Jablonowski--Staniforth 2014 state nor any steady state. Its initial
momentum residual does not converge:

```text
Linf|du/dt| [m/s2]     N=8       N=16      N=32
isothermal rest        0.0e+00   0.0e+00   0.0e+00
JW06 steady            7.19e-04  4.14e-04  3.03e-04
UMJS14 steady          2.66e-03  2.66e-03  2.65e-03
```

`2.66e-3 m/s2` is 230 m/s per day. `umjs14_baroclinic` is bit-identical to
`umjs14_steady`: `max|delta M| = max|delta(theta M)| = max|delta ps| = 0`, so the
perturbation that defines the baroclinic case does not exist.

`surface_orography.cpp` hardcodes `R = 6.371229e6`, `Omega = 7.29212e-5`, `u0 = 35`, and
`eta0 = 0.252` for the JW06 terrain, duplicating constants that the benchmark module also
defines and ignoring `config.planet`. Selecting `orography.kind=jw06` on a non-Earth
planet silently builds terrain that is not in balance with the wind.

`configs/phase6_dcmip_2_0_0.cfg` uses `a_half = [20544.8, 0, ..., 0]` with a uniform
`b_half`, so `p_half[1] = ps/15 = 6667 Pa` is below `p_half[0] = 20544.8 Pa` and the
preset fails validation at startup. The Phase 6 pressure-gradient-error measurement was
therefore taken on a substitute coordinate, where it does not converge:

```text
 N     Linf|du/dt| [m/s2]   L2 [m/s2]
 8       7.39e-4            1.12e-4    (terrain unresolved: 1 cell per half wavelength)
16       1.67e-3            1.89e-4
32       1.64e-3            3.18e-4
```

## Decision

### What "steady" means

A preset may be named `steady` only if its **initial residual decreases with
resolution**. The gate is convergence, not smallness: an absolute threshold alone accepts
a state that is uniformly wrong at every resolution, which is exactly the present UMJS14
failure. Absolute thresholds are registered after the Phase 9 pilot measurement and are
not relaxed afterwards.

### UMJS14

Implement the Ullrich--Melvin--Jablonowski--Staniforth 2014 analytic balanced temperature
and surface-pressure fields together with the prescribed baroclinic perturbation, so that
`umjs14_steady` is steady and `umjs14_baroclinic` differs from it. This closes the Phase 5
production gate.

If that implementation is not completed, the fallback is to stop claiming the name: rename
the case to `solid_body_unbalanced` and demote it to a non-balanced smoke test, and at the
same time withdraw the corresponding gate from the Phase 5 validation report and from
[docs/plan.md](../plan.md). Silently keeping the name is not an option.

### JW06 terrain and the configured planet

The JW06 surface geopotential keeps its analytic form, but its constants live in one place
and `R` and `Omega` come from `config.planet`. JW06 is defined against Earth constants, so
a config whose planet constants differ from the registered JW06 values and which selects
`orography.kind=jw06` is **rejected explicitly** rather than given unbalanced terrain.

### DCMIP 2-0-0 coordinate

The hybrid coefficients are replaced by the same decreasing ramp that
`uniform_sigma_coefficients(20544.8, 15)` produces, so interfaces are monotone. That every
shipped preset actually starts is added as a smoke gate, since a preset that cannot start
cannot be measured. The pressure-gradient-error convergence is then measured on the
repaired coordinate.

Changing the vertical discretization itself — coordinate reconstruction, a
Simmons--Burridge pressure-gradient force, or a reference-state subtraction — is out of
scope. Phase 9 produces the measurement; a later ADR may act on it.

The Phase 9 pilot registered the following absolute `Linf` ceilings for the two steady
states at `N=8/16/32`: UMJS14 `6.0e-4/3.5e-4/3.0e-4 m/s2`, and JW06
`8.0e-4/4.5e-4/3.3e-4 m/s2`. The observed DCMIP values
`7.335e-4/1.653e-3/1.625e-3 m/s2` do not converge; that signature is retained as a
regression measurement, not accepted as a balanced-state threshold.

## Consequences

- A `steady` preset that stops being steady fails a test, at every resolution triple.
- The Phase 5 UMJS14 gate and the Phase 6 DCMIP 2-0-0 gate become measurable instead of
  blocked.
- Non-Earth planets lose the JW06 terrain option, which they never legitimately had.
- The terrain-following pressure-gradient error is quantified on the intended coordinate,
  giving a later discretization ADR real evidence.

## Rejected alternatives

- **Keep UMJS14 as it is and register its residual as a tolerance.** That converts a bug
  into a contract and makes the Phase 5 convergence gate unfalsifiable.
- **Register an absolute residual threshold instead of a convergence rate.** A uniformly
  unbalanced state passes any threshold loose enough to admit the coarsest resolution.
- **Rescale the JW06 constants to an arbitrary planet.** JW06's balance is derived for its
  own constants; rescaling produces a state that is neither JW06 nor balanced.
- **Adjust the DCMIP `a_half[0]` value alone.** The failure is the non-monotone ramp, not
  the top pressure.
- **Fix the pressure-gradient discretization in Phase 9.** It would land in the same phase
  as the measurement that is supposed to justify it.

## Validation obligations

- Initial residual of every `steady` preset decreases across `N = 8, 16, 32`.
- `umjs14_baroclinic` and `umjs14_steady` have different initial states.
- `orography.kind=jw06` with non-registered planet constants is rejected.
- Every shipped preset starts, and `phase6_dcmip_2_0_0` has monotone interfaces.
- Resting-atmosphere spurious wind `Linf` and `L2` for DCMIP 2-0-0 are recorded against
  resolution on the repaired coordinate.

## References

- P. A. Ullrich, T. Melvin, C. Jablonowski, A. Staniforth (2014), [A proposed baroclinic
  wave test case for deep- and shallow-atmosphere dynamical cores](https://doi.org/10.1002/qj.2241).
- C. Jablonowski and D. L. Williamson (2006), [A baroclinic instability test case for
  atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12).
- P. A. Ullrich et al. (2012), [DCMIP 2012 Test Case Document v1.7](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf).

Commit order is recorded in [the Phase 9 implementation plan](../phase-9-plan.md).
