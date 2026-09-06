# ADR 0015: Held--Suarez vertical thermodynamic transport and energy attribution

**Status:** accepted as a mitigation; exact dry-core energy conservation remains open

## Context

The completed 1200-day Phase 7 production integration was finite and conserved dry and
tracer mass, but its day-200--1200 total-energy budget did not close. The accumulated
Held--Suarez forcing was `-7.08599e22 J`, while the unattributed residual was
`+7.07235e22 J` (`+8.18558e14 W` on average). At the final checkpoint, a centered
directional derivative of the semi-discrete dry dynamics was `+5.51299e14 W`.

The final state was evaluated with every already-implemented horizontal and vertical
reconstruction combination. Keeping linear Barth--Jespersen horizontal transport and
changing only potential-temperature vertical reconstruction from minmod-limited to
unlimited linear reduced the dry-dynamics rate to `+1.22890e14 W`, a 77.7% reduction.
Changing momentum reconstruction as well produced no material additional benefit.
The minmod limiter was active in 14,923 of 17,280 final cell-level profiles, so its
nonlinear switches were not a rare-event correction.

## Decision

The three-dimensional dry driver uses separate vertical limiter selections:

- potential temperature uses unlimited linear reconstruction;
- momentum and passive tracer retain the configured limiter;
- donor-cell selection remains donor-cell because its limiter argument is inactive;
- SSPRK3 stage validation continues to reject non-positive temperature instead of
  clipping the state.

The generic column-coupling entry point keeps its previous shared-limiter behaviour.
An explicit overload carries the potential-temperature limiter, and the dry driver uses
that overload. This keeps the monotonic tracer contract while reducing the diagnosed
thermodynamic limiter source.

The physics CSV now also accumulates explicitly configured horizontal diffusion. Its
reported residual subtracts thermal forcing, Rayleigh-drag work, and diffusion, so a
deliberately dissipative term is not mislabeled as unexplained numerical error.

## Consequences

- The existing production input file remains unchanged and continues to state
  `vertical.limiter = minmod`; that setting applies to momentum and passive tracer in
  the dry driver.
- The old checkpoint can be inspected, but the 1200-day integration must be rerun to
  measure the end-to-end improvement because changing the operator changes the
  trajectory.
- This is not a claim of exact energy conservation. The remaining instantaneous defect
  is dominated by the collocated horizontal transport/least-squares pressure-force
  pairing. Eliminating it requires a jointly derived negative-adjoint or Hamiltonian
  horizontal operator, not a global energy correction.
- A Simmons--Burridge-style vertically compatible formulation and a compatible
  horizontal dry-core state remain the next structural steps if the rerun still misses
  the production energy gate.

## Validation

The coupling regression verifies that the split limiter changes potential-temperature
transport while leaving minmod-limited momentum and tracer tendencies unchanged. The
Phase 5 dry baseline, Phase 7 forced-core test, diffusion attribution test, allocation
gate, and Phase 7 CSV smoke test pass with this decision.
