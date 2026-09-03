# Phase 1 validation report

## Scope and reproducibility

Phase 1 implements the equiangular gnomonic cubed sphere, global topology,
panel fields and halos, Cartesian tangent vectors, finite-volume operators,
prescribed winds, conservative upwind/MUSCL transport, SSP-RK3, CFL selection,
checkpoint/restart, and machine-readable diagnostics. The executable writes a
small `tracer.csv` snapshot beneath `output.directory`.

The local matrix used AppleClang 21.0.0 and GCC 15.2.0. All 23 CTest cases
passed in Debug with both compilers, in the AppleClang Release preset, and in
the AppleClang ASan/UBSan preset. `format-check` and `git diff --check` passed.

Reproduce the local gates with:

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
cmake --build build/dev --target format-check

cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan
```

Run the published cases with:

```sh
./build/dev/my_planet_sim --config configs/phase1_solid_body.cfg
./build/dev/my_planet_sim --config configs/phase1_williamson1.cfg
./build/dev/my_planet_sim --config configs/phase1_deformational.cfg
./build/dev/my_planet_sim --config configs/phase1_divergent.cfg
```

## Geometry and operators

The property gate uses `N=1,2,4,16,64`. It verifies unit positions, positive
areas and lengths, two cells per edge, outward normals, and the exact topology
counts `6N^2`, `12N^2`, and `6N^2+2`. The global area agrees with
`4 pi R^2` within the planned relative tolerance `5e-13` at every resolution.
All 24 directed panel connections round-trip, including index reversals.

Arbitrary unique-edge fluxes cancel globally to rounding error. A constant
field has bitwise-zero least-squares gradient and Laplacian. For the manufactured
field `q=x`, gradient RMS error decreases on `N=8,16,32`.

## Transport gates

A broad, globally smooth Gaussian `exp(-(1-x))` was transported for one full
rotation at CFL 0.35. The observed orders from the two finest resolutions were:

| scheme | resolutions | norm | observed order | gate |
|---|---:|---:|---:|---:|
| piecewise-constant upwind | 12, 24, 48 | L1 | 0.879 | >= 0.8 |
| unlimited linear reconstruction | 6, 12, 24 | L2 | 2.010 | >= 1.7 |

The same gate checks relative mass drift below `5e-13`. Axis-aligned and oblique
rotation cases are finite, conservative, and keep seam RMS below three times
global RMS. A stopped run resumed from its flattened checkpoint is bitwise
identical to an uninterrupted run. The Barth-Jespersen limiter keeps slotted
cylinder values inside the initial `[0,1]` range to `1e-12`.

For the checked-in `N=12` oblique-axis sample, relative mass drift was
`-7.94e-15`; seam and corner RMS errors were below global RMS. The reversible
nondivergent deformation sample had drift `-7.07e-15`, no overshoot, and a
threshold-averaged filament preservation score of `0.804`.

Diagnostics include normalized `L1/L2/Linf`, min/max, seam/corner RMS,
threshold-area filament preservation, and one-tracer unmixing, overshooting,
and in-range mixing measures. Williamson test 1 uses the cosine-bell preset;
the deformational and divergent presets use a time-reversing streamfunction
with a solid-body background flow.

The benchmark selection and diagnostic categories follow the open-access
[Lauritzen et al. (2012) test-suite paper](https://doi.org/10.5194/gmd-5-887-2012).
