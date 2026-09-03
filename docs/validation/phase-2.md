# Phase 2 validation report

## Scope and reproducibility

Phase 2 implements the rotating one-layer shallow-water equations on the
Phase 1 cubed sphere: a versioned depth/Cartesian-momentum state, the
cubed-sphere dual topology and signed incidence, tangent-vector operators,
a Rusanov unique-edge flux with depth-positive linear reconstruction, a
gravity-wave CFL driver on SSP-RK3, budgeted Laplacian/biharmonic diffusion,
the Williamson 2, Williamson 6 and Galewsky benchmarks, and a compatible
C-grid candidate with its own diagnostic operators.

The executable writes the run summary and the conservation drift to standard
output, a cell snapshot to `output.directory/shallow_water.csv` and the
interval invariants to `output.directory/diagnostics.csv`.

The local matrix used AppleClang 21.0.0 and GCC 15.2.0. All 40 CTest cases
passed in Debug with both compilers, in the AppleClang Release preset and in
the AppleClang ASan/UBSan preset. `format-check` and `git diff --check` passed.

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
cmake --build build/dev --target format-check

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan
```

The published tables come from the Release build:

```sh
cmake --preset release
cmake --build --preset release
./build/release/my_planet_sim --config configs/phase2_williamson2.cfg
./build/release/my_planet_sim --config configs/phase2_williamson2_compatible.cfg
./build/release/my_planet_sim --config configs/phase2_comparison_rusanov.cfg
./build/release/my_planet_sim --config configs/phase2_comparison_compatible.cfg
./build/release/my_planet_sim --config configs/phase2_williamson6.cfg
```

The Galewsky table below uses a reduced resolution; its command is given with
the table.

## State, topology and operators

The flatten/unflatten round trip is bitwise, and a layout or shape mismatch is
rejected before a restart. Every stage boundary checks finiteness, the depth
floor and the tangency `|dot(m,k)|`; a violation aborts with the stage, the
CFL and the minimum depth.

The dual topology gate builds `N=1..64`, checks the Euler counts, the cyclic
vertex rings including the eight three-valent cube corners, positive dual areas
and lengths, and the exact identity `D2 D1 = 0` at every cell-vertex pair.

A uniform resting state has a rounding-scale right-hand side. An arbitrary
unique-edge mass flux cancels globally, and the Coriolis tendency has the
documented sign and stays tangent. Linear reconstruction preserves a constant
state exactly, keeps near-dry face depths above the floor and does not break
the unique-edge mass cancellation.

## Linear waves, rest and balance

| gate | result |
|---|---|
| non-rotating inertia-gravity wave, `N=8,16` | depth L2 error falls by more than 20 per cent per refinement |
| Coriolis-only inertial oscillation | shows the SSP-RK3 design time order |
| geostrophic adjustment | stays finite and depth-positive |
| uniform rest, 1000 steps | mass drift below `5e-14`, no new depth extrema |

## Williamson test 2

The five-day Release runs of `configs/phase2_williamson2.cfg` (`N=12`) and
`configs/phase2_comparison_rusanov.cfg` (`N=24`) give:

| measure | N=12 | N=24 | observed order |
|---|---:|---:|---:|
| depth L1 [m] | 21.99 | 3.29 | 2.74 |
| depth L2 [m] | 26.42 | 3.68 | 2.84 |
| depth Linf [m] | 79.10 | 11.26 | 2.81 |
| velocity L2 [m/s] | 1.549 | 0.192 | 3.01 |
| mass drift | -3.87e-14 | -7.90e-14 | — |
| energy drift | -8.73e-3 | -1.02e-3 | — |
| potential enstrophy drift | -2.86e-2 | -3.55e-3 | — |
| axial angular momentum drift | -7.99e-3 | -9.16e-4 | — |
| minimum depth [m] | 1188.1 | 1108.2 | — |

The observed orders exceed the planned `1.7` for linear reconstruction, and
mass drift stays below the planned `5e-13`. The convergence gate additionally
runs the oblique flow axis `(1,1,1)/sqrt(3)` and requires its errors to stay
within twice the axis-aligned errors at the same resolution.

## Williamson test 6

The Rossby-Haurwitz wave is checked as a regression rather than against an
independent high-resolution reference. The unit gate verifies that the initial
state is positive, finite and tangent, that its depth spectrum is dominated by
wave four by more than two orders of magnitude over wave three, that the field
checksum is deterministic and sensitive to a one-ulp change, and that a short
integration keeps the mass drift below `5e-13` and at least 90 per cent of the
initial wave-four amplitude.

The published fourteen-day run of `configs/phase2_williamson6.cfg` at `N=24`:

| measure | value |
|---|---:|
| mass drift | -4.55e-13 |
| energy drift | -1.65e-2 |
| potential enstrophy drift | -1.15e-1 |
| axial angular momentum drift | -4.28e-3 |
| minimum depth [m] | 8019.8 |
| maximum depth [m] | 10331.4 |
| wall time [s] | 235.7 |

## Galewsky jet

The balanced jet is generated from the published constants with a trapezoidal
gradient-wind integral over 4096 latitude intervals, and the depth is shifted so
that the global mean matches `shallow_water.mean_depth_m`. The gate checks that
the jet is confined to the northern hemisphere between `pi/7` and `5 pi/14`,
that the mean depth is reproduced, that the height perturbation is local and
leaves the velocity untouched, and that the initial gradient-wind residual
converges: the momentum-tendency RMS of the balanced state, relative to its own
Coriolis term, falls by a factor of about four between `N=24` and `N=48`.

The published six-day run reduces `configs/phase2_galewsky.cfg` to `N=24`, which
is the same run that ADR 0003 uses for its robustness row:

```sh
sed 's/cells_per_panel = 48/cells_per_panel = 24/' configs/phase2_galewsky.cfg \
  > build/galewsky_n24.cfg
./build/release/my_planet_sim --config build/galewsky_n24.cfg
```

| measure | value |
|---|---:|
| mass drift | -1.76e-13 |
| energy drift | -2.57e-3 |
| potential enstrophy drift | -8.02e-2 |
| minimum depth [m] | 8984.9 |
| maximum depth [m] | 10152.6 |
| wall time [s] | 96.5 |

## Diffusion

Laplacian and biharmonic diffusion are globally mass neutral, damp a smooth
depth mode at the expected sign and rate, and are bitwise zero when disabled.
The coefficient helper derives the coefficient from a resolution and an
e-folding time, so no resolution-dependent constant is written into a config.
Diffusion tendencies enter the budget separately from the Rusanov dissipation,
and the budget components sum to the total.

## Scheme comparison

[ADR 0003](../adr/0003-shallow-water-horizontal-discretization.md) records the
measured comparison and the decision to keep the reference Rusanov scheme for
Phase 3 and Phase 4. In summary, both schemes conserve mass to rounding, the
compatible candidate is about five times cheaper per cell-step but one to two
orders of magnitude less accurate, its potential enstrophy drift does not fall
with resolution, and it fails the Galewsky jet outright without dissipation.
The cause is the non-orthogonality of the equiangular cubed sphere together
with the three-valent cube corners, where neither the dual circulation nor the
Bernoulli gradient converges.

## Known gaps

These items of the Phase 2 plan were reduced in scope and are not claimed:

- No seam/corner error-ratio diagnostic exists for shallow water; the Phase 1
  transport diagnostic was not extended. Grid imprinting is instead measured by
  the wave-four amplitude of a zonally symmetric state and by the operator unit
  tests.
- Williamson 6 is compared against its own regression checksums, not against an
  independently generated high-resolution reference, and no day 1/7/14
  vorticity and PV norm table is produced.
- The Galewsky jet is verified at a single resolution against its balance and
  conservation properties. The planned `N=48,96,192` self-convergence study and
  the comparison against a published spectral reference were not run, and the
  `N=48` default of `configs/phase2_galewsky.cfg` has no published table yet.
- The compatible candidate reuses the collocated state instead of a separate
  edge-staggered checkpoint layout, as recorded in ADR 0003.
