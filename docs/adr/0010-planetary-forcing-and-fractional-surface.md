# ADR 0010: Planetary forcing and fractional land--ocean surface

**Status:** accepted for bounded Phase 8 implementation

## Context

The dry hydrostatic core already reads planetary radius, signed rotation rate, gravity,
gas constants, reference pressure, and fixed orography. Phase 7 deliberately restricts
Held--Suarez forcing to its registered Earth-like constants and has no prognostic surface.
Changing that meaning in place would invalidate the Phase 7 regression and its published
comparison.

The Phase 5--7 production gates and acceptance of ADR 0008 remain open. On 2026-09-05 the
user explicitly requested the planned Phase 8 implementation; this accepts the contracts in
this ADR for bounded implementation, but does not turn the earlier CI evidence into a
production-validation claim. Those gaps remain visible in the Phase 8 validation report.

Phase 8 also needs a minimal representation of land and ocean. An ocean dynamics model is
out of scope. A binary land mask would make coastlines resolution-dependent and would
discard sub-grid islands and coast cells, so it is not an acceptable intermediate state.

## Decision

### Separate benchmark forcing from surface thermodynamics

Keep `physics.kind=held_suarez` unchanged. Add two distinct selections:

```text
physics.kind = planetary_newtonian
physics.kind = surface_energy_balance
```

`planetary_newtonian` implements only the preregistered axisymmetric and substellar analytic
benchmark profiles. It has no prognostic surface and accepts flat, uniform geography only.
`surface_energy_balance` owns a surface-temperature tendency and is used to test orbital
insolation, terrain, coastlines, and land--ocean thermal inertia. Its results are not claimed
to conform to the analytic Newtonian benchmarks.

All time-dependent sources are evaluated from the same stage state and stage time as the
complete dry-dynamics RHS at every SSP-RK3 stage.

For both analytic profiles, set `sigma=p/p_s`, `kappa=Rd/cp`, and

```text
T_eq = max(200 K,
           (315 K + 60 K chi
                  - 10 K log(p/p0) cos(phi)^2) (p/p0)^kappa)

k_T = 1/(40 day) + (1/(4 day)-1/(40 day))
                     max(0,(sigma-0.7)/0.3) cos(phi)^4
k_v = 1/(1 day) max(0,(sigma-0.7)/0.3)

dT/dt = -k_T (T-T_eq)
du/dt = -k_v u
```

where `day=86400 s`. The axisymmetric case uses `chi=-sin(phi)^2`, making it
algebraically the existing Held--Suarez profile while permitting a different rotation rate.
The substellar case uses the signed dot product of the cell direction and the planet-fixed
star direction. It is not clipped on the night side. This is the Merlis--Schneider/Heng
longitudinal dipole generalized to the shared orbit geometry.

### Fractional surface

Every horizontal cell has one immutable fraction:

```text
0 <= f_land <= 1
C_surface = f_land*C_land + (1-f_land)*C_ocean
```

No threshold converts `f_land` to a Boolean. A coast cell therefore retains
`0<f_land<1`. Pure ocean and pure land are the endpoints. A mixed cell has one surface
temperature and one area-weighted heat capacity; it does not have separate tile
temperatures or fluxes.

`C_land` and `C_ocean` are configured directly in `J m-2 K-1`. In particular, there is no
ocean depth parameter and no conversion from depth, density, or volumetric heat capacity.
Albedo, emissivity, sensible-exchange coefficient, internal heat flux, and low-level drag
are cell-independent in Phase 8. Thus heat capacity is the only land--ocean physical
difference.

The initial registered surface values are `C_land=2.0e6 J m-2 K-1`,
`C_ocean=4.0e7 J m-2 K-1`, `T_s=288 K`, `alpha=0.3`, `epsilon=1`,
`lambda_air=10 W m-2 K-1`, and `Q_internal=0 W m-2`. They define an idealized heat-capacity
ratio of 20, not an implied ocean depth. They remain explicit config values rather than
hidden code defaults.

### Surface energy equation

For upward sensible flux `H`:

```text
Q_sw = (1-alpha) S(t) max(0, r_hat dot s_hat)
H = lambda_air (T_s - T_air,bottom)

C_surface dT_s/dt = Q_sw + Q_internal - epsilon sigma_SB T_s^4 - H
d(M theta)_bottom/dt += H/(cp Pi_bottom)
```

The surface loses exactly the stage-weighted `H` used to construct the atmospheric source.
Diagnostics report absorbed stellar, emitted longwave, internal, sensible, surface storage,
the discrete atmospheric thermal contribution, and their residual separately. Temperature
clipping is not part of the method; a source step that would violate positivity is retried
with a smaller time step.

This is an idealized dry surface energy balance. It has no atmospheric shortwave absorption,
downwelling longwave, greenhouse gas, latent heat, ocean transport, or subsurface diffusion.

### Geography and Earth data

The run-owned `SurfaceBoundary` contains fixed surface geopotential and land fraction plus
one provenance record. `surface.geography=uniform` supplies a constant fraction and may use
the existing analytic/imported orography. `surface.geography=earth` supplies both fields
from one checked derived product and rejects a second non-flat orography selection.

The Earth product is generated offline from fixed versions of NOAA ETOPO 2022 Ice Surface
60 arc-second elevation and NOAA GSHHG 2.3.7 shoreline polygons. Water elevations, including
bathymetry, are set to zero. Land elevation is retained even where it is below mean sea
level. Source-cell and target-cell values are area averages. The repository stores the
bounded 1 degree derived product, generator, versions, licenses, cryptographic source hashes,
and a runtime fingerprint; the simulator never downloads data.

### Orbit and synchronization

Orbital parameters are SI values in the expanded config. A deterministic Kepler solve and
signed planet rotation produce a planet-fixed unit star vector. Synchronous rotation is an
exact contract for the zero-obliquity/eccentricity benchmark: `Omega=2*pi/P_orbit` and the
substellar vector remains fixed. Its stellar day is reported as undefined, not as an
arbitrary large finite number.

Longitude of periapsis is measured prograde from northern vernal equinox. For true anomaly
`nu`, the planet-to-star solar longitude is `nu+longitude_of_periapsis+pi`; obliquity rotates
that direction into the equatorial inertial frame. The rotation phase is anchored by the
configured substellar longitude at `t=0`. This convention, rather than only a period formula,
is part of the restart and reproducibility contract.

### Compatibility and state layout

Absent Phase 8 keys leave existing canonical configs unchanged. Existing `none` and
`held_suarez` runs keep their atmospheric payload and FrameV2 byte-for-byte. Surface runs add
an optional `surface_temperature_k[C]` tail and use checkpoint layout
`dry_hydrostatic_surface_v1`; fixed geography is reconstructed from config and fingerprint
and is not checkpointed. Phase 8 adds standalone surface CSV output but no FrameV3 or Web
protocol change.

## Consequences

- Coastline behavior converges with the fractional source and is not controlled by a
  resolution-dependent land/sea cutoff.
- Ocean thermal inertia can be changed without pretending to resolve an ocean depth or
  circulation.
- Earth terrain, land fraction, and provenance stay mutually consistent.
- The surface model is small enough for unit energy budgets, but it is not an Earth climate
  or habitability model.
- Published analytic benchmarks remain comparable because the Held--Suarez path is not
  silently generalized.
- Surface-enabled restart needs a new explicit layout, while all earlier checkpoints remain
  readable and unchanged.

## Rejected alternatives

- **Boolean land mask:** loses mixed coast cells and changes discontinuously with resolution.
- **Separate land and ocean temperatures in one cell:** introduces tile flux partitioning
  and more prognostic state before it is needed.
- **Ocean depth times volumetric heat capacity:** exposes a physical depth while omitting all
  other depth-dependent ocean processes; direct areal heat capacity states the approximation
  more honestly.
- **Negative Earth elevation as ocean:** misclassifies below-sea-level land. Shoreline
  polygons decide land/water; elevation does not.
- **Generalize `held_suarez` in place:** changes a completed Phase 7 contract and prevents an
  exact Earth regression.
- **Use the surface model as the tidally locked reference:** conflates a new radiation/surface
  approximation with a published analytic-forcing benchmark.
- **Runtime online/GDAL ingestion:** expands dependencies and weakens reproducibility.

## Validation obligations

- Exact old-config/state/checkpoint/FrameV2 identity.
- Orbital analytic cases, global insolation integral, and synchronous fixed-star property.
- Capacity interpolation at endpoints and multiple strict mixed fractions.
- Surface-only third-order time refinement and equal-and-opposite sensible exchange budget.
- Earth land area, mixed coast presence, water-zero height, seam/corner, and provenance gates.
- Earth-like regression, nondimensional similarity, zero-rotation safety, and preregistered
  slow/rapid/tidally locked comparisons.

Implementation detail, commit order, preset matrix, and reference sources are recorded in
[the Phase 8 implementation plan](../phase-8-plan.md).
