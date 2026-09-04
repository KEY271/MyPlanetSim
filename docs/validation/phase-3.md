# Phase 3 validation report

## Scope

Phase 3 adds an independent TypeScript/React visualizer, a loopback-only Node gateway,
and a versioned control/event/frame boundary for the existing native C++ simulator.
The browser preview remains usable without C++; live execution is owned by the gateway
and the native process remains authoritative for frame data.

## Reproducibility

From a clean checkout with the Phase 3 dependencies installed:

```sh
cmake --build --preset dev
ctest --preset dev
cd web
npm ci
npm run lint
npm run test
npm run build
npm run test:live --workspace @myplanetsim/gateway
```

The C++ suite completed 37/37 tests. The Web workspace completed TypeScript lint,
gateway tests (3/3), UI tests (6/6), protocol tests (5/5), and the Vite production build.
The optional native live test completed 1/1 test in about 0.25 s.

## Scientific and protocol gates

- The native live test uses `configs/phase3_rest_n4.cfg`, submits a Gaussian depth edit,
  waits for `run.accepted`, `run.started`, and `run.completed`, and fetches the published
  binary frame only after `frame.ready`.
- Frame 0 has the `MPSFRAM1` binary header. A later frame is fetched and asserted to differ
  from frame 0, so the test exercises C++ time evolution rather than a synthetic animation.
- A second native run is cancelled through the gateway and must reach `cancelled`.
- C++ unit tests cover strict control parsing, validated initial edits, frame writing, and
  observer scheduling. Existing Phase 0–2 tests remain in the 37-test regression suite.
- The request, translated control text, event log, diagnostics, and frame references are
  returned by the bundle endpoint.

## Security and resource boundaries

The gateway binds to `127.0.0.1` or `::1`, checks loopback host/origin and an in-memory
Bearer token, uses a configured preset allowlist, rejects oversized/invalid requests, and
spawns the configured executable with an argument array and `shell: false`. Frame paths are
resolved below the isolated run directory, and unpublished or symlinked files are rejected.
The gateway regression tests cover unauthorized access, path-free control translation,
shell-free argv, reconnect sequencing, and frame isolation.

## Geometry, UI, and performance

Protocol/UI tests cover frame decoding, CSV adapters, edit validation/undo, the shared
cubed-sphere geometry, 3D picking, 2D projection round trips, grid edge modes, lifecycle,
diagnostics, and bundle handling. The structural performance gate constructs the expected
`6*N*N` cells and `12*N*N` unique edges for `N=48` and `N=96`; it passed in about 0.4 s.
The UI has responsive layout and keyboard-accessible native controls. Explicit WebGL failure
messaging and a browser runner are not committed in this phase.

## Known gaps

The current browser entrypoint deliberately uses the deterministic mock client by default.
The gateway API and `HttpSimulationClient` contract are present, but the built UI does not yet
serve from the gateway or automatically receive its session token. Vite proxy wiring and a
Chromium/Firefox/WebKit offline/live matrix therefore remain follow-up work. The native
gateway live gate is still independent of that UI wiring and verifies the scientific control
path end to end.
