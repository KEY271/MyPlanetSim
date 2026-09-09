# Phase 15 validation report

Status: complete (2026-09-10)

Phase 15 replaces the interactive, instantaneous-frame visualizer with a read-only
terrain and period-mean dataset workflow. The implementation sequence is P15.01--P15.08,
based on `742f11f`. Historical Phase 2/3 plans and reports remain in `docs/`; their
shallow-water and control-path executables are no longer present.

## Implemented contract

- `PeriodMeanAccumulator` integrates accepted-step end states with overlap seconds over
  half-open fixed-length periods. Variable steps, spin-up, multi-period crossings,
  duplicates, invalid order, and non-finite values have unit coverage.
- A normal dry-hydrostatic CLI run always writes `viewer/manifest.json` and
  `viewer/terrain.bin`. When `statistics.enabled = true`, it also writes one
  `viewer/means/period_NNNNNN.bin` per non-empty period.
- `MPSTERR1` and `MPSMEAN1` files have a 24-byte little-endian header: 8-byte magic,
  uint32 schema version, uint32 field count, and uint64 total float64 value count.
  Atmospheric payloads are field-major, then cell-major, then downward-increasing model
  level. The manifest carries all offsets, shapes, units, period bounds, coverage,
  fingerprints, and declared byte lengths.
- Checkpoints with active statistics write a hash-matched `.statistics` sidecar. The
  browserless `phase15.visual_dataset_cli` gate compares a continuous run with a stop at
  step 6 and restart: both final checkpoint and both period binaries are SHA-256
  identical.
- The TypeScript reader validates safe relative paths, schema, cubed-sphere shape,
  offsets, bounds, byte budgets, binary headers, and finite values. It reads one selected
  period at a time with a bounded cache.
- The UI opens a local dataset folder and shares field, period, level, and cell selection
  across the 2D map, 3D globe, inspector, and vertical profile. It has no solver process,
  gateway, token, initial-condition editor, run buttons, or instantaneous timeline.

The final C++-to-TypeScript integration gate exposed one contract defect that isolated
fixtures had missed: the terrain writer stored cell count in the header where the reader
correctly expected total value count. P15.08 changes the writer to store
`field_count * cell_count` and adds a C++ header assertion plus a real CLI decoder test.

## Removed paths

- The shallow-water model, benchmarks, diagnostics, configuration enums/keys, dedicated
  diffusion/reconstruction, compatible operators, dual topology, Williamson/Galewsky
  presets, and their tests.
- Native control requests, initial-condition edits, FrameV1/FrameV2 I/O, NDJSON events,
  gateway process management, browser client/mock state machines, and live tests.
- Playwright/Chromium dependencies and CI steps, browser test artifacts, interactive
  preset lists, `just gateway`, and the former combined gateway/UI launcher.

ADR 0002 and 0003 are marked retired. ADR 0004 and 0007 are superseded by ADR 0028.
Their detailed text remains as historical design evidence.

## Reproduction

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev --output-on-failure
cmake --build build/dev --target format-check

cd web
npm run lint
npm run typecheck
npm run test
npm run build
cd ..

./build/dev/my_planet_sim \
  --config configs/phase15_viewer_rest_n4.cfg \
  --progress-interval-s 0
MPS_VISUAL_DATASET="$PWD/output/phase15-viewer-rest-n4/viewer" \
  npm run test:cli-dataset --prefix web --workspace @myplanetsim/protocol
```

The checked fixture has `N=4`, `K=8`, twelve accepted 10-second steps, and two complete
60-second periods. The CTest gate produces the same fixture under `build/dev/output/` and
also exercises the step-6 restart.

## Gate results

The 2026-09-10 local Debug build used GCC 15.2.0. The complete CTest run passed all 94
tests, including the new `phase15_gate`. C++ format-check passed. Web lint and typecheck
passed; Vitest passed 16 UI tests and 2 protocol unit tests, and the separate real CLI
dataset decoder passed. The Vite production build completed without errors. The remaining
bundle-size message is a warning, not a failed budget.

A repository search found no shallow-water, control request, event stream, old frame,
gateway, Playwright, or interactive-preset references in live C++, configuration, Web,
Justfile, or CI paths. References remain only where historical documents explain their
former behavior or retirement.

## Accumulation cost sample

This bounded cost sample used the same Debug binary and an enlarged `N=24`, `K=8`
isothermal-rest fixture for 120 seconds at requested `dt=10 s`. Both variants completed
12 accepted steps with the same reported dry mass (`5.14981e18 kg`) and total energy
(`1.50911e24 J`). The enabled case selected seven fields and wrote two 60-second periods.
`/usr/bin/time -l` supplied wall time and maximum RSS.

| measurement | statistics off | statistics on | difference |
|---|---:|---:|---:|
| wall time | 7.05 s | 7.27 s | +3.1% |
| maximum RSS | 40,878,080 B | 46,628,864 B | +5,750,784 B (+14.1%) |
| viewer output | 139,719 B | 2,849,634 B | +2,709,915 B |
| accepted steps | 12 | 12 | 0 |

The off case still writes the static terrain and manifest. The on-case increment is two
1,354,776-byte mean binaries plus period metadata. This is a single local sizing sample,
not a portable performance threshold. The accumulator retains only the current window;
period count increases output files, not accumulator storage. The accepted-step hook calls
diagnosis for requested visual fields but does not evaluate a full RHS or physics process.

## Manual display and remaining limits

The production bundle and all data/state logic were validated without Chromium, as
intended. Canvas/WebGL pixels were not asserted by automation; opening the generated
`viewer` folder in `just dev` remains the manual rendering smoke test.

The first schema uses fixed SI-second periods and model levels. It does not provide
calendar climatologies, arbitrary UI-side re-averaging, isobaric interpolation, variance,
extrema, anomalies, instantaneous animation, old-frame conversion, or remote dataset
hosting. Phase 14 long-duration statistical and performance work remains separate and is
not claimed complete by this cleanup.
