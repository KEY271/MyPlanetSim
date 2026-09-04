# Phase 6 validation

## Scope and result

The fixed-terrain C++ core is implemented and passes the CI-scale regression gate. The
validated scope is deliberately narrower than full Phase 6 conformance: the production
six-day and published-reference comparisons listed under **Known gaps** have not been run,
and the Phase 5 entry gate remains open. This report therefore does not mark Phase 6 as
complete.

The implementation adds one immutable, driver-owned cell-centred surface-geopotential
field; cellwise hydrostatic lower boundaries; the existing generalized pressure-gradient
evaluation over sloping hybrid surfaces; one Rusanov shallow-water terrain source; four
named analytic profiles; one strict lat-lon CSV path; one bounded smoothing family; and
terrain budget diagnostics. Terrain remains outside prognostic state and checkpoint
payloads. `web/`, FrameV2, the gateway, and the control protocol are unchanged.

## Commands and CI-scale gates

Validated on 2026-09-04 with:

```sh
cmake --preset dev
cmake --build --preset dev -j4
ctest --preset dev --output-on-failure
cmake --build build/dev --target format-check

cmake -S . -B build/gcc-phase6 -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 \
  -DMPS_WARNINGS_AS_ERRORS=OFF
cmake --build build/gcc-phase6 -j4
ctest --test-dir build/gcc-phase6 -L phase6_gate --output-on-failure

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan -j4
ctest --preset asan-ubsan -L phase6_gate --output-on-failure
```

| Gate | Result | Cost |
|---|---:|---:|
| Clang development suite | 56/56 passed | 45.87 s |
| Clang `phase6_gate` | 3/3 passed | 1.06 s |
| GCC 15.2 `phase6_gate` | 3/3 passed | 1.65 s |
| ASan/UBSan `phase6_gate` | 3/3 passed | 2.10 s |
| `format-check` | passed | less than 1 s |

The three `phase6_gate` executables cover bounded orographic responses, JW06 analytic
state/profile regression, and rigidly rotated DCMIP terrain norms. The full development
suite additionally covers the exact flat path, a uniform geopotential offset, sloping
manufactured pressure-gradient refinement, DCMIP 2-0-0 initialization, Williamson test 5
source sign and tangency, terrain torque/work and absolute pressure velocity, checkpoint
fingerprint checks, strict CSV interpolation, and smoothing invariants.

## Evidence by responsibility

| Responsibility | CI-scale evidence |
|---|---|
| Fixed field and lower boundary | shape/finiteness checks; flat exact identity; uniform-offset hydrostatic and energy properties |
| Sloping pressure surfaces | separate geopotential and pressure components; decreasing manufactured horizontal error |
| Analytic profiles | DCMIP 2-0-0, Williamson 5, linear bell, and JW06 values and bounded state initialization |
| Shallow-water coupling | exact-zero flat source; non-flat sign/unit/tangent checks; compatible-path rejection |
| Dry response and budgets | bounded short-response amplitude test; pressure-gradient, wind, torque, work, and `omega_abs` diagnostics |
| Imported terrain | strict tensor-product CSV, byte fingerprint, periodic bilinear interpolation, config-relative path |
| Smoothing | constant, area-mean, extrema, and cubed-sphere unique-edge properties |
| Rotation | DCMIP profile mean, L2 norm, and maximum agree within registered N=24 tolerances after three rigid rotations |
| Restart and compatibility | unchanged checkpoint layout; terrain regeneration/fingerprint checks; all Phase 0--5 C++ tests pass |

The DCMIP and JW06 analytic constants and equations follow the DCMIP test-case document and
the Jablonowski--Williamson paper. Williamson test 5 follows the standard shallow-water test
set. The tests are regressions of this implementation at CI scale; they are not substitutes
for the published long-duration reference comparisons.

## Toolchain note

The ordinary GCC 15.2 build and Phase 6 tests pass. Enabling `MPS_WARNINGS_AS_ERRORS=ON`
with GCC currently stops on `-Wmissing-field-initializers`, first in Phase 5
`dry_hydrostatic_state.cpp` and then in existing aggregate fixtures. Those warnings are
not a terrain runtime failure, but the strict GCC warning gate is not claimed as passing.
They are left untouched here to avoid expanding Phase 6 into a broad aggregate-initializer
cleanup.

## Known gaps

- The Phase 5 quantitative 3D convergence and UMJS14 published reference envelope required
  by the Phase 6 plan are still open.
- DCMIP 2-0-0 has not been integrated for six days at three horizontal resolutions and
  both K=15 and K=30. No production observed-order or absolute-error ceiling is claimed.
- Williamson test 5 has not been compared at days 10 and 15 with a bundled, citable
  spectral reference dataset. The current gate covers initialization, source behavior,
  bounded short evolution, and rotation sensitivity only.
- The linear-bell gate checks bounded short response and approximate amplitude linearity;
  it does not yet validate published wavelength, phase, or pressure-coordinate momentum
  flux over a reflection-free production run.
- JW06 steady and perturbed production runs have not been compared with published-day
  pressure, low-level temperature, or vorticity envelopes.
- The strict GCC warnings-as-errors gap described above remains open. The ordinary GCC
  build, Clang suite, and sanitizer gate are recorded separately rather than treating that
  warning-policy failure as a numerical result.

## References

- P. A. Ullrich et al., [DCMIP 2012 Test Case Document v1.7](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf), test 2-0-0.
- D. L. Williamson et al., [A standard test set for numerical approximations to the shallow water equations in spherical geometry](https://doi.org/10.1016/S0021-9991(05)80016-6), test 5.
- C. Jablonowski and D. L. Williamson, [A baroclinic instability test case for atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12).
