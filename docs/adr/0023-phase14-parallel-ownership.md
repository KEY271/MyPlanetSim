# ADR 0023: Phase 14 parallel ownership

## Status

Accepted and integrated.

## Context

The dry RHS used edge scatter and shared column scratch. Adding worksharing directly
would race, while unordered floating-point reductions would make thread-count effects
hard to distinguish from numerical changes.

## Decision

- Every shared edge is evaluated once into an edge-owned slot. Cells gather their four
  incident slots in increasing edge-ID order, preserving the former serial addition
  order. Tracer flux is also evaluated once per edge, tracer and level.
- State diagnosis, vertical coupling and moist physics own disjoint cell output.
  Pressure/source and scalar-to-momentum fast-operator work own disjoint level output.
- Scratch is allocated per OpenMP worker. Exceptions are captured in indexed slots and
  rethrown in cell, edge or level order outside the parallel region.
- Maximum and minimum diagnostics are reduced after worksharing in deterministic index
  order. Dependent GMRES iterations and their norm/dot reduction order remain serial.
- `schedule(static)` is the default. The runtime thread count is not capped in code;
  deployment controls it with `OMP_NUM_THREADS` and affinity variables.

## Consequences

OpenMP OFF keeps the same code ownership and gather order. OpenMP ON can parallelize
the RHS without atomics in hot loops. Full face reconstruction arrays remain until the
P14.08 fusion decision. Scaling measurements are an explicit foreground operation via
`tools/compare_phase14_openmp.py`; they are not part of configure, build or tests.
