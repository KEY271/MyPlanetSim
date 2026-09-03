# ADR 0004: interactive simulator control boundary

**Status:** accepted for Phase 3

## Context

Phase 3 adds a local browser visualizer for the existing native shallow-water
simulator. The browser needs to edit a small set of initial-condition and run
parameters, observe progress, and retrieve binary fields without becoming a
second implementation of the numerical model.

The simulator is a native C++ application with an existing human-oriented CLI.
The browser must not select arbitrary commands, binaries, configuration paths,
output paths, or environment variables. A local control gateway is therefore
needed to validate browser requests, own run directories, and manage the native
child process.

## Decision

Use three independent components with versioned file and stream contracts:

```text
browser UI -- RunRequestV1 / EventV1 / FrameV1 -- localhost gateway
                                                     |
                                  ControlRequestV1 + child process
                                                     |
                                              native C++ simulator
```

- The browser depends only on the Web protocol package and `VisualDatasetV1`.
- The gateway binds to loopback, owns the configured binary/preset/run-root
  paths, and starts the binary with an argument array and `shell: false`.
- The C++ application parses and validates `ControlRequestV1` independently;
  gateway validation is not a substitute for the native validation.
- `FrameV1` frame zero, written by C++ after initial-condition edits and before
  integration, is the authoritative initial field. Any browser-side preview is
  explicitly non-scientific and is discarded when frame zero arrives.
- The simulator emits lifecycle and frame metadata as versioned NDJSON events.
  Binary field payloads are little-endian Float64 data in files published by
  atomic rename; event metadata identifies the published relative path.
- The gateway exposes only opaque run IDs and server-resolved frame paths. It
  never forwards arbitrary shell text, filesystem paths, or environment data.
- Phase 3 is single-user and local-only. Authentication is a memory-only
  session token plus Origin/Host checks; remote deployment, multi-user queues,
  TLS/OAuth, and resource isolation are outside this ADR.

## Protocol ownership

The protocol version is `1`. The shared Web package owns browser-facing schema
types and validation. The C++ application owns the native control parser and
frame encoding. The gateway translates between the two and rejects startup if
the simulator's reported capabilities are incompatible with its allowlist.

The first production edit is the tagged `gaussian_depth` edit. It stores a
unit-vector origin, amplitude, angular width, and an explicit mass policy. The
edit is applied in C++ to the benchmark initial state; the gateway and browser
do not compute the authoritative field.

## Failure and reproducibility

Each run has an isolated server-owned directory containing the original request,
resolved control request, ordered events, stderr, frame manifest, binary
frames, build metadata, and terminal status. Cancellation is requested at a
solver step boundary and is represented by a terminal event. A gateway timeout
may force termination only after the graceful interval expires.

Unknown protocol versions, malformed or non-finite values, unsupported presets,
path escape, published-file violations, and child-process failures are explicit
errors. Human CLI output and existing CSV/checkpoint behavior remain unchanged
when machine mode is not requested.

## Consequences

The UI can run offline against deterministic fixtures and the gateway can be
tested with a fake simulator. The native solver remains the sole source of
scientific state, while the cost is a small duplicated protocol model and a
gateway process boundary that must be tested for lifecycle and local security.
