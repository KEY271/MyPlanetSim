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
