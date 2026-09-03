# Phase 0 validation report

## Result

Phase 0 is complete for local development. The scientific/software foundation builds
without warnings, passes all 16 CTest entries, passes ASan/UBSan, and demonstrates the
designed temporal orders for Forward Euler and SSP-RK3. Checkpoint serialization and a
restarted ODE integration reproduce the uninterrupted final state bit for bit.

Validation was run on commit `724dd31` on 2026-09-04 in the Asia/Tokyo timezone. The
documentation-only commit containing this report does not alter the validated binaries.

The repository has no configured remote, so the GitHub Actions workflow has not yet run
on a hosted runner. Equivalent GCC, Clang, Release, and sanitizer jobs were executed
locally; the hosted status must be confirmed on the first push.

## Toolchains and test matrix

| Toolchain | Configuration | Result |
|---|---|---:|
| GNU C++ 15.2.0 | Debug, warnings as errors | 16/16 passed |
| GNU C++ 15.2.0 | Release, warnings as errors | 16/16 passed |
| Homebrew Clang 22.1.4 | Debug, warnings as errors | 16/16 passed |
| Homebrew Clang 22.1.4 | Release, warnings as errors | 16/16 passed |
| Apple Clang 21.0.0 | `dev` preset | 16/16 passed |
| Apple Clang 21.0.0 | `release` preset | 16/16 passed |
| Apple Clang 21.0.0 | `asan-ubsan` preset | 16/16 passed |
| clang-format 22.1.4 | `format-check` | passed |

No AddressSanitizer or UndefinedBehaviorSanitizer diagnostics were emitted.

## Temporal convergence

The manufactured equation is

```text
dy/dt = -y,  y(0) = 1,  0 <= t <= 1
```

and the exact final value is `exp(-1)`. The automated gate uses the following results.

| Integrator | `dt` | Absolute error | Observed order |
|---|---:|---:|---:|
| Forward Euler | 0.2 | 4.01994e-2 | — |
| Forward Euler | 0.1 | 1.92010e-2 | 1.06599 |
| Forward Euler | 0.05 | 9.39352e-3 | 1.03144 |
| Forward Euler | 0.025 | 4.64700e-3 | 1.01537 |
| SSP-RK3 | 0.2 | 1.43957e-4 | — |
| SSP-RK3 | 0.1 | 1.66068e-5 | 3.11579 |
| SSP-RK3 | 0.05 | 1.99429e-6 | 3.05783 |
| SSP-RK3 | 0.025 | 2.44345e-7 | 3.02889 |

The gates require monotonically decreasing error and observed order `>= 0.7` for
Forward Euler and `>= 2.7` for SSP-RK3. Both pass at every refinement.

## Restart and reproducibility

- A ten-step SSP-RK3 run with `dt=0.1` was compared with a run stopped after step 5,
  serialized to the version-1 checkpoint format, deserialized, and continued.
- Final time, step, and every bit of the binary64 state match the uninterrupted run.
- `max_digits10` serialization preserves the tested finite binary64 values, including
  signed zero and the smallest/largest finite magnitudes.
- Canonical configuration text and its FNV-1a fingerprint are stable for identical
  inputs and change when a configuration value changes.
- Equal explicit `std::mt19937_64` seeds produce identical sampled sequences.
- Build metadata records model version, Git revision/dirty state, compiler, build type,
  seed, configuration fingerprint, and canonical configuration.

## Commands used

Primary preset validation:

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev

cmake --preset release
cmake --build --preset release
ctest --preset release

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan

cmake --build build/dev --target format-check
```

The GCC and Homebrew Clang checks used separate `build/gcc-*` and `build/clang-*`
directories with `MPS_WARNINGS_AS_ERRORS=ON` and both Debug and Release build types.

The convergence table can be regenerated with:

```sh
ctest --test-dir build/dev -V -L phase0_gate
```

## Delivered foundation

- CMake library/application/test target separation and reproducible presets
- dependency-free unit test support and GCC/Clang CI matrix
- scientific units, coordinates, signs, indices, and error norm conventions
- validated `Real`, validation helpers, and `PlanetParameters`
- strict deterministic experiment configuration
- checked contiguous `Field2D` and periodic halo filling
- allocation-stable Forward Euler and SSP-RK3 steppers
- compensated serial reductions and weighted error norms
- canonical run metadata, explicit random seed, and versioned checkpoints
- runnable Phase 0 ODE application with checkpoint/restart

## Known limitations carried into Phase 1

- Hosted CI remains to be observed on the first push; no Git remote exists locally.
- Windows/MSVC configuration is declared but was not tested on this macOS host.
- Reproducible reductions are currently serial; thread/MPI decomposition invariance is
  intentionally deferred.
- The flat configuration format and text checkpoint are Phase 0 formats, not the final
  large-field I/O design.
- `Field2D::operator()` uses Debug assertions while `at()` is always checked. Numerical
  kernels must use tested bounds before relying on unchecked Release access.
