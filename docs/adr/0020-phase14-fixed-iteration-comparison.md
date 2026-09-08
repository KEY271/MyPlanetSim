# ADR 0020: Phase 14 fixed ICI iteration comparison

Status: comparison contract accepted; default iteration count unchanged.

## Decision

Compare the existing ICI method at 2, 3, and 5 iterations before changing the
numerical method. Keep the diffusion, linear tolerances, divergence checks,
positivity constraints, and CFL checks unchanged. A speed target is not a reason
to tune these parameters. Keep the final full RHS evaluation and invalidate its
reuse after local physics, as in ADR 0016.

`tools/compare_phase14_ici.py` runs N=6/K=20 for one day at 450/900/1800 s,
from either the initial aquaplanet or a validated developed moist checkpoint.
Each condition has a warm-up, five ordinary runs, five profiled runs, and one
separate allocation run. Timings are compared using medians and ranges from the
same binary, config, and checkpoint; their hashes are checked across the matrix.
Selected physical diagnostics must agree exactly across instrumentation modes.

The legacy maximum physics interval remains 300 s. Its equal subdivision means
450 s dynamics actually uses two 225 s physics updates; 900/1800 s uses 300 s.
Therefore compare temperature to five iterations at the **same dt** to isolate
the iteration error. Also report the 450 s / five-iteration reference, explicitly
including its different physics subdivision. This does not replace the Phase 13
450 s / five-iteration / 75 s physics reference.

The spatial temperature RMS uses reference layer mass times cell area as its
weight. Integrated evaporation and precipitation exclude the imported cumulative
history. Near-zero flux comparisons use the manifest's absolute floor.
The reported one-day gates cover only the quantities actually measured; CAPE and
rain tails, normalized moist enthalpy closure, terrain/baroclinic comparisons,
and 30/1200-day climate validation remain separate requirements.

## Existing dry-wave constraint

The unchanged `phase10.long_step` test gives second-order convergence with two
iterations (observed orders 1.97693 and 2.04108). However, at 1800 s its normalized
iteration error is `1.85243e-9`, exceeding the temporal truncation comparison
`7.50432e-10`. Three iterations give `2.25569e-10`.
Thus second-order convergence alone does not justify changing every preset to
two iterations. Retain three for the registered dry-wave regime.

## Scope of this decision

This is the first part of P14.02, not completion of Phase 14 or adoption of a new
integrator. No predictor/IMEX experiment has been evaluated by this decision.
The fixed-iteration matrix and its limitations are recorded in the Phase 14
validation report before pursuing those separate numerical changes.
