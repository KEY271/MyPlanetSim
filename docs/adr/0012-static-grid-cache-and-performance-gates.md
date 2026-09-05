# ADR 0012: Static grid cache, RHS workspace, and performance gates

**Status:** accepted for bounded Phase 9 implementation

## Context

`configs/phase7_held_suarez.cfg` (`N=12`, `K=20`, 864 cells x 20 levels) runs at
0.220 s/step, or `7.85e4` cell-level updates per second, with clang `-O2` on one thread.
A sampling profile attributes more than 55% of self time to quantities that do not
change in time:

| share | symbols |
|---:|---|
| 43% | `atan2`, `acos`, `normalize`, `tangent_basis`, `inverse_map`, `tan`, `sin`, `map_to_unit_sphere` |
| 12% | `neighbor_across`, `cell_index`, `cell_id`, `cell_edges`, `cell` |
| 9% | `reconstruct_dry_hydrostatic_face_states` body |
| 9% | `least_squares_gradient` and `least_squares_vector_gradient` |

The cause is structural. `reconstruct_tangent_vector` re-derives `inverse_map`,
`tangent_basis`, and `logarithmic_displacement` on every call, though all three are fixed
per `(cell, edge)` pair. One RHS evaluation makes `4K` scalar plus `K` vector
least-squares passes for reconstruction and `2K` scalar passes for sources, which is 140
grid passes per RHS and 420 per step at `K=20`. `couple_dry_hydrostatic_columns`
allocates roughly 30 `std::vector` objects per cell, and `DryHydrostaticReconstruction`
allocates about 93 MB per RHS at `N=48`/`K=30`.

With the per-cell CFL of ADR 0011 the Phase 7 candidate step falls to about 264 s, so one
1200-day run costs about 24 hours and the ADR 0009 six-run matrix costs 150--250
core-hours. That, and not any missing physics, is why the matrix has never been run.

Parallelism is not the answer yet: more than half of the work is redundant, and
parallelising redundant work only parallelises waste.

## Decision

### Precompute the static grid quantities

`CubedSphereGrid` owns a derived cache built once at construction. Every entry is
time-invariant:

```text
per cell          tangent basis          (one inverse_map + tangent_basis)
per cell          1 / area_m2
per (cell, edge)  logarithmic displacement vector
per (cell, edge)  displacement components in the cell tangent basis
per (cell, edge)  neighbour flat index and edge sign
per (cell, edge)  outward unit normal and edge tangent basis
per cell          sum over faces of ( sign * L_f * n_f )
per cell          least-squares neighbour weights (wx, wy) x 4
```

The bounds-checked `neighbor_across`, `edge_sign_for_cell`, `cell_index`, `cell_id`, and
`cell_edges` remain part of the public API and stay in use outside hot loops; hot loops
read flat arrays.

### Separate memoization from reformulation

Performance commits are split by what they may do to the numbers, and the two kinds never
share a commit:

| commit | kind | gate |
|---|---|---|
| static geometry cache | memoization: same expression, same order, evaluated once | **all presets byte-exact** |
| least-squares stencils | reformulation: the `a00/a01/a11` solve folded into neighbour weights | registered tolerance, relative `1e-12` initially |
| RHS workspace | allocation removal, arithmetic order unchanged | **byte-exact** |
| observer interval | diagnostic call frequency only | **state and checkpoint byte-exact** |

Exactly one commit is allowed to move a fingerprint. That is what makes a moved
fingerprint diagnostic rather than ambiguous.

### Allocation contract

RHS evaluation performs no heap allocation in steady state.

- `DryHydrostaticDriver` owns the level slices, column slices, reconstruction buffers, and
  tendency buffers, and passes spans into the RHS.
- `couple_dry_hydrostatic_columns` and `vertical_scalar_rhs` write into caller buffers.
  `interface_coordinate` and `center_coordinate` are computed once per column and shared.
- `AtmosphericHybridCoordinate::geometry` is evaluated once per diagnose and handed to the
  physics kernels, rather than recomputed per cell in both `diagnose` and the forcing.
- The allocation count itself is pinned by a test that counts global `operator new`.

### Observer sampling

`DryHydrostaticDriver::advance` currently re-runs `diagnose` and a full physics kernel on
every accepted step whenever an observer is attached, while
`PhysicsDiagnosticsAccumulator` keeps only rows where
`step % diagnostics.interval_steps == 0` (default 144). The observer signature carries
whether derived fields and physics rates are needed for this call. The per-step
stage-weighted accumulators stay per step, since the RHS already produces them.
`ClimateStatisticsAccumulator` is a time-weighted integral, so interval sampling preserves
its definition; that this is true is checked by comparing two intervals on one run.

### Performance gates, separate from correctness gates

```text
tools/benchmark_dry_core   fixed preset, fixed step count
reports: cell-level-update/s, s/step, peak RSS per cell, allocations per RHS
```

Registered target at `N=12`/`K=20`, one thread, clang `-O2`, `phase7_held_suarez`:

| metric | current | Phase 9 target | stretch |
|---|---:|---:|---:|
| cell-level-update/s | 7.85e4 | **>= 8.0e5** (10x) | 1.6e6 (20x) |
| s/step | 0.220 | <= 0.022 | 0.011 |
| heap allocations per RHS | ~1e5 | **0** (steady state) | 0 |

Removing the 43% transcendental share and the 12% topology share, plus the allocation and
observer duplication, accounts for just over 60%; a factor of about two on the remaining
arithmetic gives the 10x commitment and the 20x stretch. The cost consequence at 1200 days
and `dt ~ 264 s` (392,700 steps):

```text
today       1 run ~24 h        six-run matrix 150--250 core-hours
after 10x   1 run ~2.4 h       six-run matrix 15--25 core-hours
```

The six runs are independent, so process parallelism alone brings the matrix to about
three hours of wall time. **Making the production matrix runnable without MPI or OpenMP is
the point of the target**, not speed for its own sake.

CI measures a small preset and fails if cell-update/s falls more than 30% below its
registered value. Absolute rates are machine-dependent, so the threshold is expressed
against a reference case measured in the same job.

Phase 9 adds no threading, but its refactor is required to leave the RHS free of hidden
global state and to keep per-cell work independent, so that it does not obstruct later
parallelisation.

The Phase 9 pilot on Apple Silicon, one thread, Release with interprocedural optimization,
measured 833,382 cell-level-update/s, 0.020735 s/step, 1,254.4 peak RSS bytes per
cell-level, and zero post-warm-up allocations per RHS over 300 steps. At the amended 200 s
production step this extrapolates to about 2.99 compute-hours for one 1200-day run.

## Consequences

- The ADR 0009 six-run matrix becomes affordable on one machine.
- A moved fingerprint has exactly one candidate cause.
- Memory per cell becomes a reported quantity, so the `N=48`/`K=30` estimate is checked
  rather than assumed.
- The public bounds-checked grid API is unchanged, so existing callers keep their
  diagnostics.

## Rejected alternatives

- **Add OpenMP or MPI first.** Half the current work is redundant; parallelising it wastes
  the same cycles on more cores and hides the redundancy behind a scaling number.
- **Cache lazily inside the operators.** Lazy per-call memoization reintroduces hidden
  mutable state in the RHS and blocks later parallelisation.
- **Combine the cache and the stencil rewrite in one commit.** Then no commit is
  byte-exact and no fingerprint change is attributable.
- **Gate performance inside the correctness tests.** A slow machine would then report a
  physics failure ([docs/plan.md](../plan.md) 4.1 keeps the two separate).
- **Switch to binary checkpoints for speed.** Not yet justified; revisited only if
  checkpoint I/O exceeds 5% of wall time in the production matrix.

## Validation obligations

- Every preset is byte-exact across the cache, workspace, and observer commits.
- The stencil commit's deviation is measured and stays inside its registered tolerance.
- A steady-state RHS makes zero heap allocations, pinned by an allocation-counting test.
- Two observer intervals give identical state and identical climate statistics.
- `cell-level-update/s` and memory per cell are recorded with their reproduction command.

Commit order is recorded in [the Phase 9 implementation plan](../phase-9-plan.md).
