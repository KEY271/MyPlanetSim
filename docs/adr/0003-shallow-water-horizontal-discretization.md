# ADR 0003: horizontal shallow-water discretization

**Status:** accepted for Phase 3 and Phase 4

## Context

ADR 0002 fixed a collocated cell-centred finite-volume reference scheme with a
Rusanov unique-edge flux, and a compatible C-grid candidate that keeps `h` on
primal cells, the normal velocity on primal edges and the circulation on the
dual cells around primal vertices.  The choice between them was deliberately
left to a measurement.

Both schemes are now implemented and share the grid, the state and checkpoint
layout, the time control, the positivity checks, the explicit diffusion and the
diagnostics.  The compatible candidate deviates from ADR 0002 in one respect:
it derives the edge normal velocity from the same cell-centred state instead of
using a separate edge-staggered layout, so no second checkpoint layout was
needed.  The prognostic degrees of freedom are therefore the reference ones,
while the operators, the Bernoulli gradient and the potential vorticity flux
follow the compatible construction.

The equiangular cubed sphere is not orthogonal: the dual edge that joins two
cell centres leans away from the primal edge normal by up to about 25 degrees.
Both the dual circulation and the Bernoulli gradient therefore carry an
explicit metric correction with a Perot-reconstructed tangential velocity.
Without it the Williamson 2 balance residual does not converge at all.  With
it, the smooth part of the grid converges at second order, but the eight
three-valent cube corners do not converge in either operator.

## Measurements

All runs use the same grid, initial state, CFL target of 0.45 and no explicit
diffusion. The reference scheme uses linear reconstruction with the
Barth--Jespersen limiter; the compatible candidate is centred and reconstructs
no face states.  Williamson test 2 is integrated for five days, and the errors are
taken against its steady analytic solution.

```sh
cmake --preset release && cmake --build --preset release
./build/release/my_planet_sim --config configs/phase2_williamson2.cfg
./build/release/my_planet_sim --config configs/phase2_williamson2_compatible.cfg
./build/release/my_planet_sim --config configs/phase2_comparison_rusanov.cfg
./build/release/my_planet_sim --config configs/phase2_comparison_compatible.cfg
```

Accuracy, Williamson 2 at day 5:

| measure | N=12 Rusanov | N=12 compatible | N=24 Rusanov | N=24 compatible |
|---|---|---|---|---|
| depth L2 [m] | 26.42 | 102.48 | 3.682 | 45.05 |
| depth Linf [m] | 79.10 | 347.5 | 11.26 | 181.0 |
| velocity L2 [m/s] | 1.549 | 22.08 | 0.1920 | 12.57 |
| observed depth order | — | — | 2.84 | 1.19 |
| observed velocity order | — | — | 3.01 | 0.81 |

Conservation, relative drift over the same five days:

| measure | N=12 Rusanov | N=12 compatible | N=24 Rusanov | N=24 compatible |
|---|---|---|---|---|
| mass | -3.9e-14 | -4.0e-14 | -7.9e-14 | -8.1e-14 |
| energy | -8.7e-3 | +1.9e-2 | -1.0e-3 | +4.8e-3 |
| potential enstrophy | -2.9e-2 | +1.4e-1 | -3.6e-3 | +1.2e-1 |
| axial angular momentum | -8.0e-3 | +1.1e-3 | -9.2e-4 | -1.3e-3 |

Robustness, Galewsky jet for six days at `N=24`:

```sh
sed -e 's/cells_per_panel = 48/cells_per_panel = 24/' \
    -e 's/scheme = rusanov/scheme = compatible/' \
    configs/phase2_galewsky.cfg > build/galewsky_compatible.cfg
./build/release/my_planet_sim --config configs/phase2_galewsky.cfg
./build/release/my_planet_sim --config build/galewsky_compatible.cfg
```

| measure | Rusanov | compatible |
|---|---|---|
| six-day run | completes | aborts in the first stage |
| minimum depth [m] | 8985 | -581 |
| mass drift | -1.8e-13 | — |
| potential enstrophy drift | -8.0e-2 | — |

Cost, Williamson 2 at day 5 on the same machine and build:

| measure | N=12 Rusanov | N=12 compatible | N=24 Rusanov | N=24 compatible |
|---|---|---|---|---|
| wall time [s] | 5.91 | 1.19 | 45.2 | 9.56 |
| cell-steps per second | 1.48e5 | 7.53e5 | 1.58e5 | 7.64e5 |

Grid artifacts are measured by the operator unit tests rather than by a
spectrum, because the accuracy result already decides the question.  For solid
rotation the dual circulation error has a median that falls by a factor of six
between `N=8` and `N=16`, while its maximum stays at about 28 per cent of the
signal at every resolution and is located at the three-valent cube corners.
The compatible Williamson 2 residual shows the same pattern: the median and
the ninetieth percentile converge at second order while the maximum grows.

## Decision

Phase 3 and Phase 4 use the reference collocated finite-volume scheme with the
Rusanov unique-edge flux, linear reconstruction and the Barth--Jespersen
limiter.

The compatible operators stay in the code base as diagnostics and as the basis
for a future candidate.  They are not offered as the production scheme, and
`shallow_water.scheme = compatible` remains available for comparison runs only.

## Consequences

- Both schemes conserve mass to rounding at the same level, so mass is not a
  discriminator.  The decision rests on accuracy, conservation of the derived
  invariants and robustness.
- The compatible candidate is roughly five times cheaper per cell-step because
  it needs neither reconstruction nor a Riemann solver, but it is one to two
  orders of magnitude less accurate at the same resolution and does not reach
  its design order, so the cost advantage does not buy anything.
- The candidate is not energy or enstrophy conserving as built.  The
  non-orthogonality correction breaks the antisymmetry that the conserving form
  requires, and the measured enstrophy drift is positive and does not fall with
  resolution.  This confirms that the conservation properties of the published
  C-grid construction may not be claimed on a non-orthogonal grid.
- Without numerical dissipation the candidate fails the unstable jet outright.
  Adding enough diffusion to survive it would change the comparison conditions,
  and the amount needed is itself a negative result.
- Re-evaluate the compatible family in Phase 4 if any of these change: the
  horizontal grid becomes orthogonal or nearly so, the corner treatment gains a
  dedicated construction, or the full set of TRiSK tangential-velocity weights
  is implemented so that the antisymmetry and the conservation proof survive on
  a skewed dual mesh.
- The reference scheme's own weakness is recorded here: its energy and
  enstrophy both decay, so long climate integrations in later phases must
  budget that dissipation explicitly rather than assume it is negligible.
