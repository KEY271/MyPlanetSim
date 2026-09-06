# Phase 10 validation

## Status and scope

Phase 10 is in progress. This report is append-only by milestone: the explicit baseline
below is fixed before the semi-implicit numerical path is enabled. Passing a later gate
does not permit changing these baseline values or weakening an earlier threshold.

The method and failure contract is [ADR 0016](../adr/0016-semi-implicit-gravity-wave-integration.md).

## P10.01 explicit baseline

Measured on 2026-09-07 on Apple Silicon with Homebrew Clang 22.1.4, Release LTO,
one thread, `configs/phase7_held_suarez.cfg`, `N=12`, `K=20`, and requested step 200 s.
The Courant candidates are the direct uncapped limits; `inf` means the corresponding
tendency is exactly zero in the sampled state.

| sample | fast wave (s) | material advection (s) | vertical transport (s) | diffusion (s) | surface (s) |
|---|---:|---:|---:|---:|---:|
| day 0 | 240.578 | inf | inf | inf | inf |
| day 1 | 238.482 | 71,914.0 | 215,163 | inf | inf |
| day 10 | 226.217 | 5,746.15 | 10,609.1 | inf | inf |
| day 30 | 225.075 | 4,592.03 | 16,591.8 | inf | inf |

At day 30 the first explicit constraint after removing the fast-wave limit is material
advection at about 76.5 minutes. This does not claim that a 76-minute step is accurate or
that the implicit/nonlinear solver will accept it. It establishes that 30 minutes is not
the dry Courant ceiling in this baseline and motivates the post-target sweep required by
ADR 0016.

The 300-step compute-only baseline was:

| metric | observed | registered Phase 9 gate |
|---|---:|---:|
| seconds/step | 0.013289 | <= 0.022 |
| cell-level updates/s | 1,300,327 | >= 800,000 |
| allocations/RHS after warm-up | 0 | 0 |
| peak RSS | 24,166,400 bytes | reported |

The measured 30-day explicit run completed 12,960 accepted steps in 169.96 s
(0.013114 s/step). I/O and checkpoint cost are excluded.

Reproduction:

```sh
cmake --preset release
cmake --build --preset release -j4 --target benchmark_dry_core
./build/release/tools/benchmark_dry_core --steps 300
./build/release/tools/benchmark_dry_core --steps 432
./build/release/tools/benchmark_dry_core --steps 4320
./build/release/tools/benchmark_dry_core --steps 12960
```

## Registered long-step search

The centered accuracy sequence is 1,800/900/450 s against an explicit 50 s reference.
After 1,800 s passes, probe 2,400, 3,600, 5,400, and 7,200 s in ascending order. For each
experiment family, record the first failure and its category:

```text
linear residual | nonlinear residual | invariant | explicit CFL |
phase/amplitude accuracy | balance | climate envelope | no speedup
```

The operational limit is the largest passing step below that failure, not simply the
largest run that remains finite. Held--Suarez retains the Phase 10 completion requirement
of median accepted step >=1,200 s and the delivery target of 1,800 s even if simpler dry
dynamics admits a longer step.

## P10.03 restarted GMRES

The standalone solver uses right preconditioning, modified Gram--Schmidt with one
reorthogonalization pass, Givens rotations, and a true-residual check after every restart
cycle. Tests cover a nonsymmetric matrix, exact diagonal preconditioning, a zero initial
residual, Arnoldi breakdown, iteration exhaustion, invalid/non-finite inputs, and zero
heap allocations after workspace warm-up. The dry solver is not connected to GMRES at
this milestone.

## P10.04 linear shallow-water Helmholtz fixture

The matrix-free cell-centred Helmholtz operator uses one conservative edge flux for its
finite-volume Laplacian and an area-aware Jacobi diagonal. Constants are preserved and
the area-weighted global Laplacian integral is roundoff zero. A known nonsymmetric-GMRES
path recovers a manufactured Helmholtz solution.

The independent linear shallow-water Crank--Nicolson fixture remains stable at wave
Courant 0.5, 1, 2, 4, and 8 for 20 steps. Its compatible edge-gradient/cell-divergence
pair preserves perturbation mass and quadratic wave energy within the registered
roundoff/solver tolerance. The 0.4/0.2/0.1 Courant refinement sequence, compared with a
0.025 reference, passes the second-order self-convergence gate. No production
shallow-water integration path is changed.

## P10.05 conditional configuration

The semi-implicit schema is accepted only when
`dry_hydrostatic.time_integrator = semi_implicit` is present. In that mode all reference
state, mode-selection, nonlinear/linear tolerance, GMRES restart, and minimum-step keys
are required and validated together. Semi-implicit keys are rejected for every other
experiment kind. Legacy dry configurations emit no new keys; their canonical round trip
and fingerprint remain unchanged. Active semi-implicit parameters are emitted
canonically and therefore participate in the fingerprint and run metadata.

## P10.06 hydrostatic reference and external mode

The configured pressure and temperature now build a horizontally uniform, isothermal
reference column through the production hybrid-coordinate and hydrostatic integration
paths. Tests recover the configured temperature from `theta * Exner`, close the layer
pressure sum, and verify zero hydrostatic reconstruction residual.

The first registered vertical basis contains the normalized external mode with phase
speed `sqrt(gamma * Rd * T_ref)`. Projection/reconstruction round-trips an external
profile. Per-cell perimeter/area Courant selection leaves it explicit at 100 s and
selects it at 1,800 s on `N=12`; exceeding the configured mode cap is an error. The dry
driver constructs and owns this immutable reference data only for active semi-implicit
configuration. Internal eigenmodes remain the explicit P10.09 extension.

## P10.07 conservative reference-linear fast operator

The dry fast operator now acts on `(delta ps, delta(Mu), delta(Mtheta))`. It uses one
centered mass flux per unique horizontal edge, the existing hybrid-`B` continuity
recurrence, centered reference-theta vertical transport, and a construction-time
linearization of the production pressure/hydrostatic diagnostic. The pressure force is
then evaluated with the production least-squares gradient. Limiting, Rusanov jump
dissipation, Coriolis, diffusion, physics, tracer, and terrain remain outside `L_ref`.

`unit.dry_hydrostatic_fast_operator` registers linearity, exact zero at the horizontal
reference rest state, area-weighted global `ps` and `Mtheta` conservation to roundoff,
and the algebraic `F = L_ref + (F - L_ref)` reconstruction. This commit does not yet
permit a large driver step; P10.08 supplies and validates the external-mode implicit
inverse first.

## P10.08 external-mode iterative Crank--Nicolson

The semi-implicit driver path now evaluates the complete nonlinear RHS at both
Crank--Nicolson time levels and applies reference-linear quasi-Newton corrections until
the configured scaled residual is met. The external-mode Schur complement is a
matrix-free two-dimensional Helmholtz problem built from the same centred column-mass
divergence and least-squares surface-pressure gradient as the correction operator. A
finite-volume Jacobi approximation preconditions restarted GMRES. Conservative thermal
back substitution uses the same edge and hybrid-`B` mass fluxes.

`unit.dry_hydrostatic_semi_implicit` directly verifies the coupled inverse and scalar
conservation. Its `N=12`, `K=20` flat dry linear wave accepts one requested 1,800 s step
even though that step exceeds the diagnosed explicit Lamb stability limit, while dry
mass and global `Mtheta` remain conserved to the registered solver tolerance. The
semi-implicit step still obeys the separately diagnosed material-advection, vertical,
diffusion, and surface-reservoir limits. Retry-on-failure and detailed solver diagnostics
remain P10.10/P10.11 work; this milestone rejects a failed solve rather than accepting
it.

## P10.09 selected internal gravity modes

The production reference operator now constructs the complete level-by-level vertical
wave matrix from the same hybrid-`B` continuity recurrence, reference-theta transport,
and pressure/hydrostatic diagnostic Jacobian used by `L_ref`. A real shifted-QR
decomposition orders its positive eigenvalues by descending phase speed and stores both
the right modal basis and its inverse. Projection followed by reconstruction round-trips
arbitrary column profiles, and every registered eigenpair satisfies the original
vertical matrix to the test tolerance.

For each mode whose requested-step Courant number exceeds the configured threshold, the
quasi-Newton inverse solves an independent two-dimensional Helmholtz problem. Modes
below the threshold retain the identity correction and are left to the outer nonlinear
iteration; exceeding the configured cap remains an error. Selecting every mode in the
four-level fixture reproduces `[I - alpha*dt*L_ref]` to `2e-8` scaled error after direct
substitution. The `N=12`, `K=20`, 1,800 s dry linear-wave gate selects more than the
external mode and accepts the step while preserving the scalar integrals established by
P10.08.

## P10.10 convergence, invariants, CFL, and retry

Every semi-implicit attempt now begins from the last accepted instantaneous state. A
mode-cap failure, linear or nonlinear residual failure, candidate explicit-CFL
violation, or prognostic/diagnosed invariant violation rejects the attempt and retries
with half the step. No failed candidate is copied into the accepted state. If the next
halving would cross `semi_implicit.minimum_time_step_s`, the driver stops with the final
rejection reason instead of accepting a failed solution. A shortened final remainder is
still allowed, as specified by the run end time.

The retry gate forces the 1,800 s, 20-level case above a one-mode cap, verifies that a
shorter attempt is accepted, and reproduces that accepted state bit-for-bit by running
the shorter step directly from the same initial state. Raising the minimum step above
the first retry leaves the caller state unchanged and produces a fatal error. A separate
3,600 s test compares one uninterrupted call with two 1,800 s calls and obtains the same
step count, step sequence, and prognostic arrays bit-for-bit without storing tendency
history.

## P10.11 solver and Courant diagnostics

Each accepted semi-implicit step now carries the configured/requested and accepted time
steps; material-advection, implicit-wave, and vertical Courant numbers; selected mode
count; total and maximum per-mode GMRES iterations; maximum linear residual; nonlinear
iteration and residual; retry count; and wall time split among full RHS evaluation,
modal linear solves, and the complete step. Rejected-attempt work is included in the
iteration, retry, and timing totals so the reported cost cannot hide recovery work.

Standalone runs write sampled values to `semi_implicit_diagnostics.csv` without changing
legacy explicit CSV schemas. Run metadata additionally records every reference vertical
mode phase speed in descending order; the active semi-implicit configuration remains in
the canonical configuration/fingerprint block established by P10.05.

The production modal path now uses caller-owned divergence edge fluxes, Helmholtz
diagonals, mode masks, and persistent GMRES callables. The allocation gate warms the
workspace once and measures zero heap allocations on the following multi-mode solve.
`benchmark_semi_implicit CONFIG [--steps COUNT]` reports seconds per model day, accepted
step range, GMRES work, retries, RHS/linear/unattributed time, and peak RSS for the
long-step presets registered at P10.13.

## P10.12 balanced, terrain, and baroclinic gates

The nonlinear Crank--Nicolson residual is now evaluated in the same cell-centred
tangent space used by accepted momentum states. This removes the radial component
introduced when three-dimensional momentum fluxes from differently oriented cells are
scattered before the state is projected. The discarded component is normal to the
sphere and therefore is not a prognostic degree of freedom. Failed nonlinear steps now
also report their last residual and attempted step in the rejection reason.

The vertical eigensolver first reduces the nonsymmetric structure matrix to upper
Hessenberg form before shifted QR iteration. This makes the real positive spectrum
converge for both the flat reference and the DCMIP terrain-following hybrid coordinate.
The registered balance test then verifies:

- exact stationarity of flat isothermal rest and DCMIP 2-0-0 terrain rest for two
  requested 1,800 s steps;
- temporal refinement of the linear mountain wave over 7,200 s, with normalized state
  differences `2.36514e-7` for 1,800/900 s and `6.07442e-8` for 900/450 s;
- 1,800 s versus 450 s normalized differences below `5e-4` for UMJS14 steady,
  JW06 steady, UMJS14 baroclinic, and JW06 baroclinic initial states (observed range
  `2.4822e-7`--`3.04763e-7`).

The complete Clang development suite passes 84/84 tests. The six `phase10_gate` tests
also pass with GCC 15.2 under `-Werror`, and with the Clang ASan/UBSan preset. Aggregate
initializers touched by the Phase 10 state/config additions explicitly initialize their
remaining members so both compilers enforce the same warning-clean contract.
