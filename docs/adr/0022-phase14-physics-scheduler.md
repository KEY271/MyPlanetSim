# ADR 0022: Phase 14 physics event scheduler

## Status

Accepted and integrated.

## Context

Phase 13 applies boundary-layer, dry-convective, SBM and saturation processes on one
common substep. Phase 14 needs independent maximum intervals without scheduling by
step number or leaving accepted increments pending across dynamics-step boundaries.

## Decision

`physics.schedule = legacy` remains the implicit default and preserves existing config
fingerprints. `process_intervals` enables an event union with these keys:

- `boundary_layer.maximum_update_interval_s`
- `convection.diagnostic_interval_s`
- `convection.update_mode`
- `radiation.diagnostic_interval_s`

Intervals are finite positive maxima. Each active process runs at its own deadline and
again at the dynamics-step end. Coincident deadlines are one event. At an event the
integration order is radiation, boundary layer/surface, dry adjustment, SBM,
saturation adjustment, then precipitation transfer. Saturation adjustment remains
eligible after every temperature or moisture update.

The old `moisture.maximum_physics_substep_s` and process-interval keys are mutually
exclusive. This prevents a config from silently selecting one of two clocks. The
`cached_relaxation` spelling is reserved but rejected until the SBM cache owns the
necessary reference state and invalidation metadata.

## Consequences

- A 1800 s step with BL 600 s, convection 900 s and radiation 1800 s has events at
  600, 900, 1200 and 1800 s.
- A 450 s dynamics step shortens all larger configured intervals to 450 s.
- Event construction is independent of accepted step count and is deterministic.
- The driver executes BL and SBM deadlines through the event union. A failed event
  restores the complete state and accepted diagnostics before bisecting all processes
  due at that event.
- In process-interval mode gray radiation is removed from every dynamics RHS stage.
  One column diagnosis at each radiation event is applied over that event's actual
  interval, before BL/surface. Its physical evaluation time is the event deadline.
- Radiation and local processes finish at every accepted dynamics-step boundary. No
  cache or unapplied increment crosses that boundary, so the moist checkpoint schema
  does not change. Restarts at step boundaries reproduce the uninterrupted state.
- A radiation update that exceeds its diagnosed explicit stability limit enters the
  same atomic event rollback and bisection path as local physics.
