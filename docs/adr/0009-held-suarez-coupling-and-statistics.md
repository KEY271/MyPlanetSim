# ADR 0009: Held--Suarez coupling and climate statistics

**Status:** accepted for Phase 7 implementation

## Context

Phase 7 needs one deliberately small dry-physics benchmark. It is not a general physics
framework: the only supported process is the standard Held--Suarez (1994) Newtonian
temperature relaxation and lower-atmosphere Rayleigh drag, coupled to the existing dry
hydrostatic state.

The Phase 5/6 production entry gate is still open. On 2026-09-05 it was explicitly deferred
because completing the long production matrix first would delay the bounded Phase 7
implementation substantially. This deferral does not turn CI-scale Phase 5/6 evidence into
production validation, does not accept ADR 0008, and must remain visible in the Phase 7
validation report.

## Decision

### Selection and pairing

Add one optional configuration key:

```text
physics.kind = held_suarez
```

Omitting it means `none`; in that case canonical configuration text and its fingerprint are
byte-for-byte unchanged. `held_suarez` is valid only with
`experiment.kind=dry_hydrostatic`, `dry_hydrostatic.test_case=held_suarez`, flat
orography, and the registered Earth constants used by the Phase 5/6 presets. No forcing
constant is configurable.

### Forcing

For latitude `phi`, full-level pressure `p`, surface pressure `ps`,
`sigma=p/ps`, and `kappa=Rd/cp`, use

```text
T_eq = max(200 K,
           (315 K - 60 K sin(phi)^2
                  - 10 K log(p/p0) cos(phi)^2) (p/p0)^kappa)

k_T = 1/(40 day) + (1/(4 day) - 1/(40 day))
                     max(0, (sigma - 0.7)/0.3) cos(phi)^4
k_v = 1/(1 day) max(0, (sigma - 0.7)/0.3)

dT/dt = -k_T (T - T_eq)
du/dt = -k_v u
```

Here `day` is exactly 86400 SI seconds. Temperature forcing is converted with the current
full-level Exner function to a tendency of `M*theta`; drag is converted to a tendency of
Cartesian tangent `M*u`. Neither source changes surface pressure, layer air mass, tracer
mass, or terrain. The kernel returns tendencies and diagnostics and never mutates state.

The complete dynamics-plus-physics RHS is evaluated from the same state at every SSP-RK3
stage. There is no splitting, subcycling, exact relaxation update, or configurable coupling
interval. `none` bypasses physics entirely and must preserve the Phase 6 RHS, accepted
states, checkpoint bytes, canonical config, and output.

### Initial state

The selected benchmark starts from a flat, hydrostatic, zero-wind, uniform 264 K atmosphere.
Add a deterministic 0.1 K peak-to-peak temperature perturbation. For each `(cell, level)`,
mix the unsigned 64-bit `run.random_seed`, cell index, and level index with SplitMix64; map
the high 53 bits directly to `[0,1)`, then to `[-0.05 K,0.05 K)`. Subtract the cell-area
weighted mean independently at each level. This fixes the bit construction and avoids
standard-library random distributions. The subtraction may change the strict amplitude
bound by the magnitude of the removed mean; zero mean and toolchain reproducibility are the
contract.

### Source budgets

At every configured diagnostic sample, write `physics_diagnostics.csv` with instantaneous
global rates for applied potential-temperature mass, eastward and northward momentum,
thermal energy, drag kinetic-energy work, and their sum. Also write cumulative thermal,
drag, and total physics energy contributions. Cumulative production budgets are reset at
day 200. Model-energy residual is the measured discrete energy change minus the
time-integrated physics contribution; it deliberately includes all dynamics/numerical
energy change. Mass, tracer mass, minimum temperature, maximum wind, and non-finite count
are health columns.

The energy rates are diagnosed from the tendencies actually added to `M*theta` and `M*u`:
thermal power is the directional derivative of the current discrete internal plus
hydrostatic potential energy, including the geopotential response above a heated layer,
and drag work is `area * dot(u, d(Mu)/dt)`. Thus Rayleigh work must be non-positive;
Newtonian relaxation may have either sign.

### Online climate statistics

Use `2*N` bins uniformly spaced in `mu=sin(phi)`. A cell belongs to exactly one bin according
to its centre, including seams and corners. At each accepted state after the spin-up
boundary, first compute cell-area-weighted bin/level moments, then weight that sample by the
elapsed accepted `dt`. Accumulate only means of `p,T,u,v,u^2,v^2,T^2,uv,vT`, layer mass,
and northward mass transport. East and north are defined against the planet rotation axis.

`climate_statistics.csv` has one row per `(latitude bin, native model level)`, ordered from
south to north and top to bottom. Eddy quantities use deviations from each instantaneous
zonal mean. Eddy kinetic energy is one half of the sum of the zonal and meridional velocity
variances. Hadley mass streamfunction is the northward layer mass transport accumulated
downward from the model top; positive values denote clockwise circulation in a
latitude--pressure section. No pressure interpolation or analysis checkpoint is added.

The standard window is day 200 through day 1200. A restart used for equality testing must
occur before day 200, allowing statistics and production-budget accumulation to restart
from zero without a second checkpoint schema.

## Preregistered validation

CI gates exercise formulas, signs, tangency, deterministic initialization, exact `none`
identity, source-only third-order time refinement, bounded forced evolution, CSV shape, and
restart before a shortened spin-up boundary. They do not claim climate conformance.

The production candidate is `N=12`, 20 uniformly spaced sigma levels between 1000 Pa and
100000 Pa, a requested 600 s step, linear/Barth--Jespersen horizontal reconstruction,
linear/minmod vertical transport, and no horizontal diffusion. A pilot must establish wall
time and sustained eddy activity before this candidate is frozen. The planned comparison is
seeds 0/1/2, a 300 s step, `N=16` with 30 levels, and one measured horizontal-diffusion
setting. Because the shared dry core currently parses but does not apply its diffusion
selection, the diffusion member and therefore the six-run matrix stay blocked rather than
silently comparing identical runs.

Before seeing MyPlanetSim's 1200-day result, register these deliberately broad sanity
envelopes from the published plots:

- a pair of tropospheric eastward midlatitude jets, one in each hemisphere;
- maximum time/zonal-mean eastward wind between 20 and 45 m/s, located between 20 and
  60 degrees absolute latitude and between sigma 0.15 and 0.55;
- time/zonal-mean temperature between 190 and 320 K everywhere, with a warmer tropical
  lower atmosphere than the polar lower atmosphere;
- nonzero midlatitude eddy momentum flux, eddy heat flux, eddy kinetic energy, and
  temperature variance in both hemispheres;
- opposite-signed hemispheric Hadley cells reaching at least sigma 0.5, with equatorial
  antisymmetry no worse than 20% of their larger peak;
- each seed's second 500-day mean differs from its first 500-day mean by no more than 10 K
  in temperature and 10 m/s in zonal wind at every reported bin/level.

These are structural acceptance bounds, not digitized reference fields. A later change to
them requires a new ADR or an explicit amendment explaining the independent evidence; a
failed run is not grounds for loosening them.

## Amendment (Phase 9, 2026-09-05)

Two of the assumptions above have been measured and are wrong. This amendment changes the
production candidate and unblocks the six-run matrix. It changes nothing about the forcing
constants, the 200-to-1200-day window, or the registered sanity envelopes.

- **The candidate's time step becomes about 264 s, not 600 s.** The dry core limited its
  horizontal step per edge rather than per cell, which is roughly four times too weak on a
  quadrilateral cell. The candidate reports CFL 0.45 but runs at an effective cell CFL of
  about 1.02. [ADR 0011](0011-vector-operators-and-time-step-normalization.md) normalizes
  all three solvers to the per-cell condition; the candidate's step follows from it rather
  than from a request. No production run existed at the old step, so nothing is
  invalidated.
- **The diffusion member is no longer blocked.** The dry core now applies
  `dry_hydrostatic.diffusion_*` level by level, as ADR 0006 already permitted, using the
  corrected spherical vector Laplacian. `diffusion_kind = none` remains byte-exact with
  every earlier run, so the five existing members are unchanged.

The consequence is cost, not meaning: at the corrected step a 1200-day run is about
392,700 steps. [ADR 0012](0012-static-grid-cache-and-performance-gates.md) registers the
performance work that makes the six-run matrix affordable without MPI or OpenMP. The
matrix itself is still not run by Phase 9, which delivers one pilot run and its measured
cost.

## Consequences

- Phase 7 adds no generic process registry and no new prognostic or checkpoint field.
- Restart equality inside the statistics window is intentionally unsupported.
- Terrain plus Held--Suarez is rejected rather than assigned an unvalidated meaning.
- The complete production claim remains unavailable while the deferred Phase 5/6 gate
  remains open. The dry horizontal-diffusion gap is closed by the amendment above.

## References

- I. M. Held and M. J. Suarez (1994), [A Proposal for the Intercomparison of the
  Dynamical Cores of Atmospheric General Circulation Models](https://www.gfdl.noaa.gov/bibliography/related_files/ih9401.pdf).
- H. Wan, M. Giorgetta, and L. Bonaventura (2008), [Ensemble Held--Suarez Test with a
  Spectral Transform Model: Variability, Sensitivity, and Convergence](https://doi.org/10.1175/2007MWR2044.1).
- K. Heng, K. Menou, and P. J. Phillipps (2011), [Atmospheric circulation of tidally
  locked exoplanets: a suite of benchmark tests for dynamical solvers](https://doi.org/10.1111/j.1365-2966.2011.18315.x).
- N. J. Mayne et al. (2014), [Using the UM dynamical cores to reproduce idealised 3-D
  flows](https://doi.org/10.5194/gmd-7-3059-2014).
- D. E. Sergeev et al. (2023), [Simulations of idealised 3D atmospheric flows on
  terrestrial planets using LFRic-Atmosphere](https://doi.org/10.5194/gmd-16-5601-2023).
