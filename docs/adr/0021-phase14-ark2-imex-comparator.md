# ADR 0021: Phase 14 ARK2 IMEX comparator

## Status

Accepted as an independent comparison method; rejected for production use after
the integrated screening.

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

The dry-hydrostatic comparison path adapts the model state and complete fast operator,
solves every vertical mode, checks each stage, and uses the stage weights for budgets.
The ARK2 update is not the last implicit stage, so the retained final full-RHS check
raises its nonlinear RHS count from three to four. This preserves the existing final
state/CFL check and first-same-as-last reuse contract.

## Consequences

- The independent reference method and its accounting exist without SUNDIALS linkage.
- The model integration retains all three full nonlinear stage evaluations and reports
  fast-operator calls separately.
- No preset, tolerance, diffusion or physical parameter changes are justified here.
- At N=6/K=20, one moist day and dt=1800 s, a single Release screening run was 4.8%
  slower than ICI2 while both used four full RHS evaluations per step. It also increased
  linear iterations from 126 to 329 and added 144 fast-operator calls per day.
- This fails the pre-registered requirement to beat an accurate ICI2 path by at least
  10%. ARK2 is therefore not adopted and is not coefficient- or tolerance-tuned.
