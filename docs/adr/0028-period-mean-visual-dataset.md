# ADR 0028: Period-mean visual dataset

Status: accepted (Phase 15)

## Context

The interactive gateway and instantaneous FrameV1/FrameV2 stream tied visualization to
a running simulator and retained one full grid for every displayed time. The supported
atmospheric workflow instead needs reproducible inspection of terrain and precomputed
period means without a simulator process or browser automation.

## Decision

- The model accumulates accepted-step end states with a time-weighted rectangle rule.
  A step contributes the duration of its overlap with each half-open period
  `[start_time_s + j*period_s, start_time_s + (j+1)*period_s)`.
- Period boundaries never alter the dynamics time step. Rejected steps, stages, and
  zero-duration notifications do not contribute.
- Each dataset contains `manifest.json`, one `terrain.bin`, and one binary file per
  non-empty period. Numeric arrays are little-endian float64 and use the model's
  cell-major, level-minor ordering. Every binary has a magic, schema version, shape,
  and byte count.
- The manifest records the grid and vertical coordinate, field identifiers and units,
  scheduled and actual period bounds, accumulated seconds, coverage, completeness,
  fingerprints, and the averaging rule. Terrain-only datasets have no periods.
- Accumulator state is checkpointed in a sidecar and matched to the model checkpoint.
  A legacy checkpoint without a sidecar may restart, but its first period is partial.
- The viewer reads only this dataset contract. It does not start or control simulations,
  edit initial conditions, retain instantaneous histories, or reinterpret a model-level
  mean as an isobaric mean.

## Consequences

Memory use is independent of the number of model steps and periods. Output grows with
the number of fields, cells, levels, and completed periods. Existing diagnostics remain
the numerical-validation interface; the visual dataset is an additional read-only
product. FrameV1/FrameV2 and the control protocol are intentionally not compatible with
this schema.
