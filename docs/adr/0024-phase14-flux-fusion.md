# ADR 0024: Phase 14 reconstruction and flux fusion

## Status

Accepted and integrated.

## Context

The production RHS materialized left and right primitive states for every edge and
level, then scanned the arrays again to compute flux. Multi-tracer flux was evaluated
twice, once while scattering to each adjacent cell.

## Decision

Gradient and limiter preparation is stored per cell and level. The edge-owned kernel
reconstructs its two faces and immediately computes the dynamics and all tracer fluxes.
Only the edge flux buffer crosses into the deterministic cell gather. The retained
`reconstruct_dry_hydrostatic_face_states` API remains a fixture/debug path and is not
part of production RHS storage.

The prepared and edge scratch byte counts, plus the eliminated face-buffer size, are
reported by profiled benchmarks. Runtime adoption does not depend on a target speedup,
and no limiter or numerical tolerance is changed.

## Consequences

Every tracer flux is evaluated once per edge and level. The production workspace no
longer owns global face states. Prepared and retained reconstruction are required to
match exactly in the unit fixture for piecewise-constant and limited-linear paths.
