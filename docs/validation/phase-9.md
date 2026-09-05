# Phase 9 validation

## Scope and result

Phase 9 is complete at the bounded scope registered in the implementation plan. The
incorrect spherical vector Laplacian, per-edge dry CFL, ignored dry diffusion setting,
invalid DCMIP coordinate, unbalanced UMJS14 initialization, tautological residuals, hot
RHS allocation/recomputation, and unbounded gateway run retention are corrected. No new
physics, parallel runtime, Riemann solver, or pressure-gradient discretization was added.

The Phase 5--7 production integrations remain deliberately unrun. Phase 9 removes their
code and cost blockers and measures a 300-step pilot; the 1200-day cost below is an
extrapolation, not a climate-validation result.

## Reproduction commands

```sh
cmake --preset dev
cmake --build --preset dev -j8
ctest --preset dev --output-on-failure
cmake --build build/dev --target format-check

cmake --preset release
cmake --build --preset release -j8 --target benchmark_dry_core
./build/release/tools/benchmark_dry_core --steps 300

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan -j8
ctest --preset asan-ubsan --output-on-failure \
  -R 'phase9.balanced_benchmarks|unit.vector_operators|unit.time_step_normalization|unit.dry_hydrostatic_diffusion|unit.rhs_allocation|unit.surface_energy_balance|phase7.forced_dry_core|phase8.planetary_core'

cmake -S . -B build/gcc-phase9 -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 -DMPS_WARNINGS_AS_ERRORS=OFF
cmake --build build/gcc-phase9 -j8
ctest --test-dir build/gcc-phase9 --output-on-failure \
  -R 'phase9.balanced_benchmarks|unit.vector_operators|unit.time_step_normalization|unit.dry_hydrostatic_diffusion|unit.rhs_allocation|unit.surface_energy_balance|phase7.forced_dry_core|phase8.planetary_core|smoke.presets'

npm test --prefix web
npm run typecheck --prefix web
```

## Gate matrix

| Gate | Result | Observed cost |
|---|---:|---:|
| Clang development suite | 75/75 passed | 55.60 s |
| GCC 15.2 focused Phase 9/regression set | 9/9 passed | 8.70 s |
| ASan/UBSan focused Phase 9/regression set | 8/8 passed | 13.66 s |
| Gateway | 16/16 passed | 0.22 s |
| UI | 28/28 passed | 0.99 s |
| Protocol | 23/23 passed | 0.44 s |
| TypeScript typecheck | passed | — |
| `format-check` | passed | <1 s |

The development suite increased from 74 tests / 55.48 s before P9.12 to 75 tests /
55.60 s, an observed increase of 0.12 s. GCC retains the pre-existing aggregate
initializer warnings, so warnings-as-errors is disabled on that secondary compiler path.
The first Web attempt failed only because the sandbox denied loopback `listen`; the same
suite passed outside that restriction.

## Correctness evidence

The vector-operator gate checks degree 1 and 2 rotational and divergent harmonic fields at
`N=8/16/32`. Both refinement pairs must have observed order at least 1.0, and the finest
Rayleigh-quotient error must be below 1% of the analytic eigenvalue. Laplacian and
biharmonic diffusion both have negative kinetic-energy tendency. The dry, shallow-water,
and transport CFL tests use the same per-cell sum over faces.

Initial momentum residuals, in `m/s2`, are:

| case | norm | N=8 | N=16 | N=32 |
|---|---|---:|---:|---:|
| UMJS14 steady | L2 | 2.135e-4 | 8.539e-5 | 6.243e-5 |
|  | Linf | 5.524e-4 | 3.180e-4 | 2.712e-4 |
| JW06 steady | L2 | 3.177e-4 | 1.243e-4 | 8.600e-5 |
|  | Linf | 7.186e-4 | 4.117e-4 | 3.016e-4 |
| DCMIP 2-0-0 | L2 | 8.564e-5 | 1.416e-4 | 2.203e-4 |
|  | Linf | 7.335e-4 | 1.653e-3 | 1.625e-3 |

UMJS14 and JW06 decrease at both refinement steps and are protected by absolute ceilings.
The repaired DCMIP preset now starts, but its terrain-following pressure-gradient error
does not converge. The gate records that failure signature instead of weakening the
steady-state criterion. A later ADR must address the pressure-gradient discretization.

Dry Laplacian diffusion at `1e5 m2/s` reduces vorticity and kinetic energy without changing
mass; `diffusion_kind=none` produces exactly zero diffusion tendency and preserves the
registered baselines. Hydrostatic and surface-energy residual tests perturb their stored
state so the diagnostics demonstrably become nonzero. Global dry mass, theta mass, tracer
mass, energy, and angular momentum use compensated accumulation.

## Performance and pilot cost

Apple Silicon, one thread, Release with interprocedural optimization, fixed
`phase7_held_suarez` (`N=12`, `K=20`), 300 steps:

| metric | observed | registered gate |
|---|---:|---:|
| cell-level updates/s | 833,382 | >=800,000 |
| seconds/step | 0.020735 | <=0.022 |
| allocations/RHS after warm-up | 0 | 0 |
| peak RSS | 21,676,032 bytes | reported |
| peak RSS/cell-level | 1,254.4 bytes | reported |

At the amended 200 s production step, 1200 days is 518,400 steps. Linear extrapolation
from the pilot is 10,749 s, or about 2.99 hours per run. Six serial runs are about 17.9
core-hours if all six retain the baseline size. Accounting for the preregistered half-step
member and the `N=16`, `K=30` member gives about 25.9 core-hours under linear scaling;
because the runs are independent, process-level parallel execution is sufficient.
Checkpoint and output costs are not included in these compute-only extrapolations. Scaling
the observed RSS/cell-level to `N=48`, `K=30` gives about 520 MB before fixed overhead,
which is practical on the measured host.

For the uniform 264 K initial state, the measured Rusanov wave speed is 325.724 m/s.
Using each cached edge-center distance, its Laplacian-equivalent coefficient
`nu_num = 0.5 lambda Delta-x` is `1.022e8` to `1.358e8 m2/s`, mean `1.259e8 m2/s`.
That mean is 1,259 times the explicit `1e5 m2/s` Laplacian setting used by the diffusion
gate. This is a scale comparison, not an assertion that the two operators have identical
spectra; it supports a later Riemann-solver decision while Phase 9 keeps Rusanov unchanged.

## Unblocked and deferred work

- Every shipped preset parses, validates, and constructs its initial driver state.
- The dry diffusion member is no longer a duplicate of `none`; the Phase 7 six-run matrix
  has no implementation blocker.
- UMJS14/JW06 balanced-state gates and the intended DCMIP coordinate are now measurable.
- The Phase 5 quantitative long integrations, Phase 6 six-day comparisons, and Phase 7
  1200-day climate matrix remain their respective phases' production-validation work.
- DCMIP pressure-gradient convergence remains a known numerical defect. New physics,
  HLL/HLLC, PGF redesign, and OpenMP/MPI remain outside Phase 9.
