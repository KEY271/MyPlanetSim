# Phase 5 validation

## C++ core gate

Validated with:

```sh
cmake --preset dev
cmake --build --preset dev -j4
ctest --preset dev --output-on-failure
build/dev/my_planet_sim --config configs/phase5_isothermal_rest.cfg
```

The unit gate covers cell-major flattening, strict configuration, derived hybrid-pressure
columns, constant-scalar Rusanov consistency, uniform pressure-gradient cancellation,
the Phase 4 continuity recurrence, tangent projection, checkpoint layout, and deterministic
budget reductions. The coupled rest preset completes 60 steps to 600 s with unchanged
state and finite dry-mass and energy diagnostics.

The DCMIP and UMJS presets in this phase are bounded initialization/regression cases. They
do not claim full published-test conformance or a reference envelope; quantitative long-run
convergence at production resolution remains a known validation gap. The reference core is
explicit and intentionally has no implicit solver, terrain, physics, multi-tracer registry,
or vertical mixing.

## Visualizer gate

The browser-facing gate is recorded separately below when FrameV2, gateway, level slicing,
and selected-column tests are run. Passing this section alone is not a Phase 5 numerical
validation claim.
