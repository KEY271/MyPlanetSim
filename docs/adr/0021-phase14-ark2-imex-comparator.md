# ADR 0021: Phase 14 ARK2 IMEX comparator

## Status

Accepted as an independent comparison method; not adopted by a production preset.

## Context

P14.02 requires a known second-order additive Runge--Kutta comparator before a
predictor/IMEX candidate can be evaluated in the dry-hydrostatic driver. The split is

`y' = N(y,t) + L(y-y_ref)`, where `N = F_dyn - L` and `L` is the existing fixed
reference-linear fast operator. Calling `L` is not a replacement for evaluating `N`:
every `N` stage still requires a full nonlinear dynamics RHS.

The comparator must fix its coefficients, stage times, full-RHS count and implicit
solve count independently of any performance result. SUNDIALS is a coefficient source,
not a new runtime dependency.

## Decision

Use the SUNDIALS ARKODE v7.7.0 default second-order additive pair
`ARKODE_ARK2_ERK_3_1_2` / `ARKODE_ARK2_DIRK_3_1_2`. It has three stages at
`c = [0, 2-sqrt(2), 1]`, three explicit evaluations, three cheap linear-operator
evaluations and two non-trivial diagonally implicit solves. The shared diagonal is
`gamma = 1-1/sqrt(2)`.

`Ark2Imex` is a fixed-step, vector-valued independent comparator. Its caller supplies
the explicit RHS, the linear implicit RHS and the solve of
`[I - gamma*dt*L] z = r`. Unit tests pin the complete table, stage times, call counts,
split linear second-order convergence and failure atomicity.

The production dry-hydrostatic driver is unchanged by this decision. A later commit
must adapt the model state and fast operator, then compare wave accuracy, conservation,
explicit constraints and end-to-end cost. In particular, the ARK2 update is not the
last implicit stage, so a separate final full-RHS check would raise its nonlinear RHS
count from three to four. Omitting that check requires cheap state/CFL checks and
stage-weighted budgets first.

## Consequences

- The independent reference method and its accounting exist without SUNDIALS linkage.
- The initial model integration should retain all three full nonlinear stage
  evaluations; fast-operator calls must be reported separately.
- No preset, tolerance, diffusion or physical parameter changes are justified here.
- If the complete model comparison is not more accurate or faster in its natural form,
  ARK2 may be recorded as not adopted without coefficient fine tuning.
