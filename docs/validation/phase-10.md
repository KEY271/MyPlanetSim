# Phase 10 validation

## Status and scope

Phase 10 is complete as of 2026-09-07. This report is append-only by milestone: the explicit baseline
below is fixed before the semi-implicit numerical path is enabled. Passing a later gate
does not permit changing these baseline values or silently weakening an earlier threshold;
any result-based amendment must remain explicit in the milestone that motivates it.

The method and failure contract is [ADR 0016](../adr/0016-semi-implicit-gravity-wave-integration.md).

One earlier decision was replaced rather than weakened. P10.14 removed the nonlinear
residual tolerance from the acceptance test and replaced it with a fixed iteration count,
because the tolerance was never the accuracy contract in the first place; the reasoning,
the literature it follows, and the measurements that fix the new iteration count are in
the P10.14 section. Every threshold that survives is unchanged. P10.15 records one explicit
acceptance amendment: isolated shortening at the still-explicit advective CFL is accepted
when it is bounded, does not persist, and is not accompanied by solver-work drift.

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

## P10.13 long-step presets

`configs/phase10_linear_wave_semi_implicit.cfg`,
`configs/phase10_dcmip_rest_semi_implicit.cfg`,
`configs/phase10_linear_mountain_wave_semi_implicit.cfg`, and
`configs/phase10_held_suarez_semi_implicit.cfg` request 1,800 s. The Phase 5--8 explicit
presets are unchanged. All four presets were revised at P10.14 to the fixed-iteration
contract; their current registered settings are `nonlinear_iterations = 3`,
`reference_update = fixed`, and `linear_relative_tolerance = 1e-4`.

## P10.14 iteration contract and measured cost

### Why the tolerance-based contract was replaced

The scheme first delivered at P10.13 accepted a step only when the scaled max-norm
Crank--Nicolson residual fell below `1e-8`. On `phase10_held_suarez_semi_implicit.cfg`
that never happened at 1,800 s. The measured behaviour was:

| quantity | P10.13 as delivered |
|---|---:|
| seconds per model day | 13.49 |
| explicit 200 s baseline, same session | 5.67 |
| accepted step, mean / maximum | 778 s / 900 s |
| retries per 111 accepted steps | 139 |
| limiting residual component | `potential_temperature_mass` |

The iteration was not diverging; it was being cut off. With retries suppressed the
residual fell from `8.88e-4` to `1.74e-7` in four iterations and was still contracting at
roughly a factor of three per iteration. Every requested 1,800 s step therefore paid its
full iteration budget, was rejected on the residual test, and was recomputed at 900 s.
Scanning `reference_temperature_k` over 220--300 K moved the final residual only between
`6.4e-7` and `1.7e-7`, and raising `maximum_implicit_modes` from 5 to 20 left the step
sequence bit-for-bit identical, so neither the reference scalar nor the modal truncation
was responsible.

The literature contract is different from the one that was implemented. Benard (2003)
gathers the classical semi-implicit scheme and its iterative variants into a single ICI
class "differing only by their number of iterations", lists the operational choices as
`N_iter = 1` for the extrapolating semi-implicit scheme and `N_iter = 2` for Cullen's
predictor/corrector and Cote et al., and describes the fully converged limit as unable to
be "achieved in practice for numerical models". Thuburn et al. (2014) report that "four
iterations are sufficient to achieve stable results, and it might be feasible to use
fewer iterations operationally" and solve the Helmholtz problem with a single multigrid
sweep, "more than sufficient accuracy in the context of the iterative nonlinear solver".
ENDGame solves its linear system four times per step at an operational tolerance of
`1e-4`. None of them tests a nonlinear residual.

### The registered contract

`semi_implicit.nonlinear_iterations` now performs exactly that many quasi-Newton
corrections. The residual is computed and reported every step but does not gate
acceptance; only a non-finite residual, or a final residual above the initial one,
rejects an attempt. `linear_relative_tolerance` is `1e-4` in all Phase 10 presets.

### Choosing the iteration count

Measured on `phase10_linear_wave_semi_implicit.cfg` at `N=4`, 7,200 s of integration,
against a six-iteration reference at the same step, with the temporal truncation error
taken as the 1,800 s versus 900 s difference at six iterations:

| iterations | iteration error | ratio to truncation error `7.50e-10` |
|---:|---:|---:|
| 2 | 1.853e-9 | 2.47 |
| 3 | 2.256e-10 | 0.30 |
| 4 | 2.377e-11 | 0.032 |
| 5 | 3.607e-12 | 0.005 |
| 6 | 3.678e-13 | 0.0005 |

The contraction is about one order of magnitude per iteration, the rate Thuburn et al.
report. Two iterations restore second-order accuracy -- the observed order over
1,800/900/450/225 s is 1.977 and 2.041 -- but leave an iteration error 2.5 times the
truncation error. Three iterations put it below. **The registered production value is
therefore `nonlinear_iterations = 3`**, not the two that the operational ICI catalogue
would suggest for a different discretization. `phase10.long_step` asserts both the
second-order convergence and that the iteration error at three iterations is below the
truncation error.

### Reference update

`semi_implicit.reference_update = per_step` rebuilds the horizontally uniform reference
column each step from the mass-weighted horizontal mean of the state, following Thuburn
et al. Measured on Held--Suarez at two iterations it did **not** improve the residual:

| configuration | mean residual | maximum residual |
|---|---:|---:|
| `fixed`, `T* = 200 K` | 1.600e-5 | 6.260e-5 |
| `fixed`, `T* = 264 K` | 1.204e-5 | 5.695e-5 |
| `fixed`, `T* = 300 K` | 9.830e-6 | 3.655e-5 |
| `fixed`, `T* = 340 K` | 1.069e-5 | 3.651e-5 |
| `fixed`, `T* = 440 K` | 1.137e-5 | 2.379e-5 |
| `per_step` (any `T*`) | 1.221e-5 | 5.781e-5 |

A horizontally uniform reference cannot represent the pole-to-equator thermal structure
that dominates the remaining residual, so tracking the horizontal mean does not help. The
`fixed` sequence instead reproduces the classical guidance Benard records from Simmons et
al. (1978): a warmer reference with larger static stability is better, with a shallow
optimum near 300 K. The presets use `fixed` with `T* = 300 K` for Held--Suarez.

`per_step` is retained for a different, measured reason: it makes the result exactly
independent of `reference_temperature_k`, which removes a hand-tuning failure mode on
planets where no good fixed value is known. It costs about 4 % in wall time.

### Performance envelope

`N=12`, `K=20`, Release LTO, one thread, one model day, median of three runs in a single
session. The explicit baseline is remeasured in the same session because the machine was
20--45 % slower than when the P10.01 baseline was recorded, and a speedup quoted against
a stale baseline is not a measurement.

| metric | registered target | explicit 200 s | semi-implicit 1,800 s |
|---|---:|---:|---:|
| seconds per model day | -- | 5.86 | **1.014** |
| speedup versus explicit | >= 3x | 1.0x | **5.78x** |
| accepted median step | >= 1,200 s | 200 s | 1,800 s |
| retries per model day | -- | -- | **0** |
| GMRES iterations, maximum | < 40, p95 <= 20 | -- | **1** |
| nonlinear iterations | fixed | -- | 3 |
| peak RSS | <= 2x explicit | 24.9 MB | 25.7 MB (1.03x) |
| solver allocations after warm-up | 0 | 0 | 0 |

Wall time splits as 64 % full nonlinear RHS, 15 % modal linear solves, 21 % residual
assembly, correction and invariant checks. The elliptic solve was never the constraint:
it converges in a single GMRES iteration. Cost is set by the number of full RHS
evaluations, which the fixed-iteration contract reduced from about 609 per model day to
about 144.

### Long-step ladder

Requested step probed in ascending order on `phase10_held_suarez_semi_implicit.cfg`:

| requested step | with registered `maximum_implicit_modes = 5` | with the cap released to 20 |
|---:|---|---|
| 1,800 s | accepted, 0 retries | accepted |
| 2,400 s | accepted, 0 retries | accepted |
| 3,600 s | rejected: required modes exceed the cap | accepted |
| 5,400 s | rejected: required modes exceed the cap | accepted |
| 7,200 s | rejected: required modes exceed the cap | accepted |
| 10,800 s | -- | rejected: candidate violates an explicit CFL constraint |

The registered operational limit for this preset is **2,400 s**, because the mode cap is
part of the registered configuration and a run that needs more modes is rejected rather
than silently left partly explicit. With the cap released the ladder ends at the material
advection limit that Phase 10 deliberately does not make implicit, which is the intended
failure mode. Accuracy at steps beyond 1,800 s is not registered as a capability here;
only stability was probed.

### Thirty-day comparison

`phase10_held_suarez_semi_implicit.cfg` at 1,800 s against `configs/phase7_held_suarez.cfg`
at 200 s, both for 30 model days:

| quantity | semi-implicit | explicit |
|---|---:|---:|
| accepted steps | 1,440 | 12,960 |
| wall time | 27.5 s | 183.2 s |
| dry mass drift | -1.99e-16 | -6.32e-13 |
| total energy, J | 1.318258e24 | 1.317910e24 |
| minimum temperature, K | 205.59 | 204.70 |
| maximum wind, m/s | 30.03 | 25.83 |
| non-finite samples | 0 | 0 |

Total energy differs by 2.6e-4 relative over 30 days. Dry mass is conserved about three
thousand times more tightly by the semi-implicit path, which takes nine times fewer
steps. Zonal statistics are not compared at 30 days because the climate accumulator does
not begin until day 200; that comparison belongs to the pilot below.

## P10.15 long pilot and Phase 10 closure

### Reproduction

The measurements above are reproduced with:

```sh
cmake --preset release
cmake --build --preset release -j8
# explicit baseline and semi-implicit cost, same session, one model day each
./build/release/tools/benchmark_dry_core --steps 432
./build/release/tools/benchmark_semi_implicit configs/phase10_held_suarez_semi_implicit.cfg --steps 48
# accuracy, iteration ladder and accepted step ladder
./build/dev/tests/phase10.long_step
# thirty-day comparison
./build/release/my_planet_sim --config configs/phase7_held_suarez.cfg --progress-interval-s 0
```

The 1,200-day pilot is the registered long-run gate:

```sh
./build/release/my_planet_sim \
  --config configs/phase10_held_suarez_semi_implicit.cfg \
  --progress-interval-s 0
```

At the measured 1.014 s per model day this is about twenty minutes on one thread. It
writes `output/phase10_held_suarez_semi_implicit/semi_implicit_diagnostics.csv`,
`physics_diagnostics.csv`, and `climate_statistics.csv`. The pilot passes when:

- the run completes at `result.time_s = 103680000` without a non-finite sample;
- `retry_count` summed over the run stays at zero, and `nonlinear_iterations` equals the
  configured 3 on every accepted step;
- `linear_iterations_maximum` shows no upward drift between the first and last hundred
  days, which is the check that solver work does not degrade with the evolving state;
- `accepted_dt_s` remains 1,800 s except for the final remainder;
- `dry_mass_kg` drift stays at the 1e-15 level reported for 30 days;
- `climate_statistics.csv` is populated, since the accumulator starts at day 200.

```sh
awk -F, 'NR>2{n++;r+=$14;if($12!=3)b++;if($10>m)m=$10}
  END{printf "steps=%d retries=%d wrong_iterations=%d max_gmres=%d\n",n,r,b,m}' \
  output/phase10_held_suarez_semi_implicit/semi_implicit_diagnostics.csv
```

### Result and acceptance decision

The pilot completed on 2026-09-07 from clean commit
`fc055a1e1d23b29f065512113068e3feb36b3b08`. The archived run is
`output/phase10_1200day-20260907-101001` and finished with exit status zero at
`result.time_s = 103680000` after 57,620 accepted steps. Wall time was 1,257.45 s,
or 1.048 s per model day. The final dry-mass drift was `-1.99e-16`, the maximum sampled
mass excursion was `5.97e-16`, and every sampled state had `non_finite_count = 0`.
`climate_statistics.csv` contains all 480 latitude-bin/level rows with an accumulated
day-200--1200 window of 86,400,000 s.

Linear solver work did not degrade. The mean/maximum
`linear_iterations_maximum` was 2.080/3 in the first hundred days and 1.601/2 in the
last hundred days; the whole-run maximum was 3 and p95 was 2. Every recorded accepted
step used the configured three nonlinear iterations, and the largest recorded linear
relative residual remained below `1e-4`.

The run did not satisfy the provisional literal expectation of zero retry and an exactly
1,800 s accepted step everywhere except the remainder. A short wind maximum near model
day 412 reached 54.28 m/s and activated the deliberately explicit advective CFL limit.
Recorded accepted steps in that episode include 886.96 s with one retry, 1,751.83 s at
advective CFL 0.45, and 894.65 s with one retry. The run returned to 1,800 s immediately
after the event, retained a sampled median of 1,800 s, and accumulated only 20 extra steps
over the ideal 57,600-step count. This is bounded operation of the registered fallback,
not a growth of solver work or a long-time stability defect. Phase 10 therefore accepts
the pilot and is complete; the zero-retry/only-final-remainder wording above remains as
the preregistered expectation and this paragraph is the explicit result-based amendment.

The day-200--1200 minimum time/zonal-mean temperature was 181.29 K in the polar top
level, below ADR 0009's 190 K climate-conformance envelope. This is not attributed to the
semi-implicit path: the final instantaneous minimum is 179.80 K, while the independent
200 s explicit 1,200-day checkpoint gives 180.06 K. The Held--Suarez 200 K value is a
40-day Newtonian-relaxation target rather than a hard floor. Scientific climate
conformance and attribution of the common upper-polar cold bias remain Phase 7 work and
do not block the Phase 10 numerical-method completion claim.

### Known gaps

- Retry diagnostics are sampled at the configured diagnostics interval and do not retain
  the rejected attempt's exception text. The archived CSV proves that at least two retries
  occurred, while the 20-step excess proves additional shortening; exact reason-count
  attribution would require per-step aggregate counters. The coincidence with the run's
  wind maximum and advective CFL 0.45 identifies the explicit advective limit as the
  operational cause.
- Accuracy beyond 1,800 s is unregistered. The ladder in P10.14 probed stability only;
  2,400 s is registered as the operational limit for the Held--Suarez preset because of
  the mode cap, not because 3,600 s was shown to be inaccurate.
- The remaining nonlinear residual is dominated by `potential_temperature_mass`, and a
  horizontally uniform reference column cannot reduce it. Bringing the iteration count
  below three would require a reference operator that carries horizontal thermal
  structure, which changes the single-eigendecomposition assumption of ADR 0016 and is
  therefore Phase 11 work.
- The 30-day comparison covers conservation and extrema. Zonal-mean circulation
  conformance against the explicit path remains Phase 7's responsibility, as ADR 0016
  states.
- `per_step` reference update is implemented and tested but not used by any registered
  preset.
