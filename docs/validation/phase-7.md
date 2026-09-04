# Phase 7 validation

## Scope and result

The bounded Phase 7 implementation is complete at CI scale. It adds exactly one optional
dry-physics selection, the standard Held--Suarez temperature relaxation and lower-level
drag kernel, unsplit source coupling at every SSP-RK3 stage, source-budget attribution,
deterministic seeded initialization, and bounded online climate statistics. The atmospheric
state and checkpoint layouts, FrameV2, control protocol, and `web/` are unchanged.

This is not a completed production climate-validation claim. On 2026-09-05 the user
explicitly deferred the Phase 5/6 production entry gate because of its cost. ADR 0008
therefore remains proposed. The 1200-day, six-member Phase 7 climate matrix is also not run:
the production candidate has not passed its pilot cost/activity gate, and the dry core's
configured horizontal-diffusion selection is parsed but not applied. Reporting six runs in
that state would either be prohibitively expensive or make the diffusion sensitivity a
duplicate, so both are left as explicit known gaps.

## Implemented commit sequence

| Commit | Responsibility |
|---|---|
| `4e1aabe` | ADR 0009 equations, constants, coupling, statistics, seed, and preregistered envelopes |
| `07f0b2c` | additive `none|held_suarez` selection, strict pairing, preset, and deterministic initializer |
| `50a978a` | pure Held--Suarez tendency kernel and formula/property tests |
| `bba4645` | forcing in every complete SSP-RK3 RHS, stage-weighted energy contributions, and budget CSV |
| `b8889e9` | equal-area latitude bins, online zonal/eddy moments, streamfunction, and climate CSV |
| `972032a` | coupled short-run, output, invariant, restart, compiler, and sanitizer gates |
| P7.07 (this commit) | corrected discrete energy attribution and recorded validation evidence plus the production-matrix deferral |

Validation of the budget output exposed that temperature forcing changes both internal
energy and the diagnosed hydrostatic geopotential above each heated layer. The final
implementation therefore evaluates the exact directional derivative of both terms. A
central-difference property test fixes that discrete attribution.

## Reproduction commands

The final validation used:

```sh
cmake --preset dev
cmake --build --preset dev -j4
ctest --preset dev --output-on-failure
cmake --build build/dev --target format-check

cmake -S . -B build/gcc-phase7 -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DMPS_WARNINGS_AS_ERRORS=OFF
cmake --build build/gcc-phase7 -j4
ctest --test-dir build/gcc-phase7 -L phase7_gate --output-on-failure

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan -j4
ctest --preset asan-ubsan -L phase7_gate --output-on-failure
```

| Gate | Result | Observed cost |
|---|---:|---:|
| Clang development suite | 62/62 passed | 67.13 s |
| Clang `phase7_gate` | 4/4 passed | 1.57 s |
| GCC 15.2 `phase7_gate` | 4/4 passed | 1.66 s |
| ASan/UBSan `phase7_gate` | 4/4 passed | 1.76 s |
| `format-check` | passed | less than 1 s |

GCC uses `MPS_WARNINGS_AS_ERRORS=OFF` because the pre-existing aggregate-initializer warning
gap documented by Phase 6 remains open. No new Phase 7 runtime or sanitizer failure was
observed.

## Evidence by contract

| Contract | CI-scale evidence |
|---|---|
| Formula and constants | equator, pole, sigma 0.7, upper-atmosphere drag, and 200 K floor values are checked directly |
| Conserved-variable conversion | source arrays have `C*K` shape; temperature uses full-level Exner; momentum remains tangent |
| Source signs | equilibrium rest is exact zero and Rayleigh work is non-positive |
| Discrete energy attribution | thermal power matches a central directional derivative of total internal plus hydrostatic potential energy |
| Stage coupling | observer receives a nonzero stage-weighted thermal contribution; source-only SSP-RK3 error ratios exceed 7 under step halving |
| `none` compatibility | source path is bypassed and repeated unforced RHS arrays are exactly equal; the Phase 5 baseline remains green |
| Initialization | same seed is bitwise repeatable, another seed differs, zero wind is exact, and each level's area-weighted temperature perturbation mean is below `1e-12 K` |
| Online moments | constant analytic fields recover exact means and zero eddy moments; chunked accumulation agrees; spin-up clipping and CSV order are fixed |
| Coupled health | the 30 s forced run moves temperature toward equilibrium while preserving mass, tracer, pressure bounds, temperature floor, and finiteness |
| Restart | a restart at 10 s, before a shortened 20 s statistics boundary, reproduces the uninterrupted 30 s state, post-boundary budget, and statistics exactly |
| Output | the native short run writes finite, deterministically ordered budget and climate CSV files with one header plus `2*N*K = 8` climate rows |

For the native 20 s output fixture, dry mass changed by about `6e-16` relatively, tracer
mass remained exactly zero, Rayleigh work was non-positive, and all reported fields were
finite. After attributing hydrostatic potential energy, the 20 s model-energy residual was
`7.15e12 J`, about `6.1e-6` of the magnitude of the accumulated physics contribution. This
residual includes the existing dynamics and time-discretization contribution and is not
treated as a conservation error of the externally forced system.

The CI fixture starts at day 200 solely to exercise CSV accumulation at the fixed boundary.
It is initialized there and has not undergone a 200-day spin-up, so its climate values are
format and reduction evidence only.

## Output contract

`physics_diagnostics.csv` separates instantaneous rates from accumulated energy. Its
production accumulation starts at day 200 and includes applied potential-temperature mass,
eastward/northward momentum, thermal energy, Rayleigh work, their total, measured model
energy change, residual, and health fields.

`climate_statistics.csv` contains one south-to-north, top-to-bottom row per native
`(latitude_bin, level)`. The `2*N` latitude bins are uniform in `sin(latitude)`. Reported
fields are time/zonal means of pressure, temperature, eastward and northward wind;
instantaneous-zonal eddy momentum and heat flux; eddy kinetic energy; temperature variance;
top-down meridional mass streamfunction; and accumulated time. No full state history is
retained.

## Known gaps

- Phase 5 quantitative 3D convergence and UMJS14 production envelopes remain open.
- Phase 6 DCMIP 2-0-0, Williamson 5, linear mountain-wave, and JW06 production gates remain
  open, and ADR 0008 remains proposed.
- The `N=12`, `K=20`, 600 s production candidate has not completed a cost/stability/eddy
  pilot and is not frozen.
- The 3-seed baseline, half-step, higher-resolution, and horizontal-diffusion six-run matrix
  has not been run. The last member additionally requires the existing dry-diffusion config
  gap to be resolved outside the deliberately narrow Phase 7 implementation.
- Consequently, published temperature/wind envelopes, Hadley cells, four eddy metrics,
  hemispheric symmetry, stationarity, and ensemble spread have not been accepted.
- Restart equality is gated with a shortened pre-window CI case, not from a production
  atmospheric checkpoint through day 1200.

## References

The forcing and comparison contract is preregistered in
[ADR 0009](../adr/0009-held-suarez-coupling-and-statistics.md), with references to Held and
Suarez (1994), Wan et al. (2008), Heng et al. (2011), Mayne et al. (2014), and Sergeev et
al. (2023).
