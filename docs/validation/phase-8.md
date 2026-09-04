# Phase 8 validation

## Scope and result

The bounded Phase 8 implementation is complete at CI scale. It adds signed-rotation
orbital geometry, safe planetary-scale calculations, immutable fractional land--ocean
boundaries, two deliberately separate forcing paths, an optional one-temperature surface
reservoir, and surface-aware checkpoint/CSV output. Existing configurations keep the old
dry state and checkpoint layout when no surface model is selected.

This is not a production planetary-climate validation claim. The checked Earth data is an
explicit analytic CI fixture, not a derivative of NOAA ETOPO 2022 and GSHHG 2.3.7. The
long slow-, rapid-, and synchronous-rotator comparisons were not run, and the inherited
Phase 5--7 production gaps remain open. The implementation therefore establishes bounded
contracts and regression gates without assigning physical credibility to the fixture runs.

## Implemented commit sequence

| Commit | Responsibility |
|---|---|
| `f37f7e0` | P8.01 forcing, fractional-surface, pairing, and provenance contracts |
| `c438053` | P8.02 Kepler/orbit geometry and safe planetary-scale calculations |
| `067bae5` | P8.03 immutable `SurfaceBoundary` and exact fractional heat capacity |
| `1738aa5` | P8.04 deterministic Earth-format CI fixture, generator, manifest, and remap |
| `c35d003` | P8.05 axisymmetric/substellar analytic planetary forcing and presets |
| `9a2165e` | P8.06 zero-depth surface energy balance and SSP-RK3 coupling |
| `9462510` | P8.07 optional surface checkpoint layout and standalone CSV output |
| `d1ac44a` | P8.08 coupled Earth-fixture, restart, similarity, and toolchain gates |
| P8.09 (this commit) | bounded validation evidence and explicit production deferrals |

## Reproduction commands

```sh
cmake --preset dev
cmake --build --preset dev -j4
ctest --preset dev --output-on-failure
ctest --preset dev -L phase8_gate --output-on-failure
cmake --build build/dev --target format-check

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan -j4
ctest --preset asan-ubsan -L phase8_gate --output-on-failure

cmake -S . -B build/gcc-phase8 -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DMPS_WARNINGS_AS_ERRORS=OFF
cmake --build build/gcc-phase8 -j4
ctest --test-dir build/gcc-phase8 -L phase8_gate --output-on-failure
```

| Gate | Result | Observed cost |
|---|---:|---:|
| Clang development suite | 69/69 passed | 77.48 s |
| Clang `phase8_gate` | 7/7 passed | 14.14 s |
| GCC 15.2 `phase8_gate` | 7/7 passed | 15.69 s |
| ASan/UBSan `phase8_gate` | 7/7 passed | 36.42 s |
| `format-check` | passed | less than 1 s |

GCC uses `MPS_WARNINGS_AS_ERRORS=OFF` because the pre-existing aggregate-initializer
warning gap remains visible in old code and tests. No runtime or sanitizer failure was
observed in the Phase 8 gate.

## Evidence by contract

| Contract | CI-scale evidence |
|---|---|
| Orbit geometry | initial longitude, `n-Omega` motion, synchronous fixed direction, eccentric flux ratio, night-side clipping, and global incoming power are checked |
| Scale singularities | Rossby, Burger, deformation-radius, thermal-Rossby, and time-scale quantities carry an explicit `defined` state instead of numerical NaN/inf |
| Old-path compatibility | axisymmetric planetary forcing produces the exact Held--Suarez rates and driver RHS; surface fields are absent unless selected |
| Fractional boundary | pure ocean, pure land, and `0.25/0.5/0.75` mixed capacities retain exact linear interpolation without thresholding |
| Surface physics | flux signs, exact evaluated-flux budget identity, radiative equilibrium, the 20:1 land/ocean response ratio, and coupled atmosphere/surface advancement are checked |
| Stage coupling | the surface and lowest atmospheric layer both change under the complete SSP-RK3 source path while pressure and other state components retain their contracts |
| Earth-format fixture | all target fractions are bounded; the cubed-sphere remap contains pure ocean, pure land, and mixed coast cells; water-only height is exactly zero |
| Restart | flattening and restoring at 60 s reproduces the uninterrupted 120 s coupled state exactly |
| Similarity | the registered length/time transformation preserves four dimensionless scale values to `1e-14` |
| Output | the native run emits finite, ordered `surface_state.csv` and `surface_diagnostics.csv`; checkpoint continuation was also compared byte-for-byte with an uninterrupted run |

The checked 1-degree fixture has 64,800 source cells, 1,078 fractional cells, a mean raw
land fraction of `0.40938657407407408`, and heights from exact zero to
`2738.67953633 m`. Its SHA-256 is
`82e706ae2d4dfd4610d0c85d8790913d6ca312bdaa44f1098248dff5d7b0a14f`.
These numbers test file integrity and coupling behavior only; they are not Earth-data
validation statistics.

## Output contract

Surface-enabled runs select checkpoint layout `dry_hydrostatic_surface_v1`, appending one
surface temperature per horizontal cell after the unchanged atmospheric state. Runs with
no surface reservoir retain the old layout. `surface_state.csv` reports geometry, height,
land fraction, mixed heat capacity, and final surface temperature in panel-major order.
`surface_diagnostics.csv` reports instantaneous absorbed stellar, internal, outgoing
longwave, sensible, storage, and residual powers at diagnostic samples.

The scale API distinguishes undefined quantities, including zero-rotation Rossby/Burger
values and a synchronous stellar day. Phase 8 does not yet serialize configured or
diagnosed characteristic scales into run metadata; callers can compute them through the
typed API without manufacturing infinities.

## Known gaps

- Phase 5 quantitative 3D convergence and UMJS14 gates, Phase 6 production mountain/JW06
  gates, and the Phase 7 1200-day climate matrix remain open. ADR 0008 is therefore still
  proposed rather than accepted.
- The repository does not contain the NOAA ETOPO 2022 and GSHHG 2.3.7 raw inputs. The
  checked product is deliberately marked `ci_fixture_not_production_noaa_derivative`; raw
  hashes, polygon-area intersection, geographic land-area acceptance, rotation histogram,
  and production Earth remap remain unverified.
- Slow-, rapid-, and tidally locked preset files parse and their kernels are gated, but no
  long integrations, jet statistics, day/night pattern, ensemble spread, or published
  benchmark-envelope comparisons were performed.
- The coupled gate is 120 s at small resolution. It checks finiteness, positivity, restart,
  and data shape, not diurnal amplitude/phase or climate stationarity.
- Surface diagnostics are instantaneous sampled rates. Stage-weighted interval
  accumulation, the atmospheric discrete thermal-energy contribution, and a combined
  interval energy residual are not yet emitted.
- A dedicated source-only SSP-RK3 third-order convergence test and the registered
  `Omega -> 0` coupled integration remain open; orbital zero/synchronous singularities are
  covered at unit level.
- Characteristic scale inputs and diagnosed `U_rms` Rossby number are not yet part of
  canonical run metadata/output.
- Phase 8 adds no ocean depth, ocean dynamics, moisture, clouds, realistic radiation,
  FrameV3, or Web UI surface fields, as required by the scope boundary.

## Data provenance boundary

The production source choices remain NOAA NCEI ETOPO 2022 60 arc-second Ice Surface and
GSHHG 2.3.7, as preregistered by ADR 0010. The exact status, generator command, product
fingerprints, source URLs, and restriction are recorded in
[`data/earth/manifest.txt`](../../data/earth/manifest.txt). Replacing the analytic fixture
with a production derivative requires obtaining and hashing both raw datasets and passing
the deferred geographic gates; it must not silently reuse the fixture product identifier.
