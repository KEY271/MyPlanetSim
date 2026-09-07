# ADR 0019: Dilute water and moist-physics coupling

**Status:** accepted for staged implementation

## Context

Phase 13 adds a closed, liquid-only water cycle to the dry hydrostatic model.  The
existing state, fluxes, diffusion, time integrators, diagnostics, and checkpoint
layout all assume one passive tracer.  Phase 12 also established an accepted-step
boundary for implicit boundary-layer mixing and dry convective adjustment.  Moist
physics must preserve that boundary: condensation and convection are finite state
changes and must not be evaluated inside Runge--Kutta stages or nonlinear iterations.

The reference implementation used to define the Simple Betts--Miller lineage is
ExeClim/Isca master commit `d1321ac372698bcbe054da1cc5b4eef9f998a096`, in particular
`qe_moist_convection.F90`.  It is a behavioural reference only.  MyPlanetSim remains
an independent C++ implementation with no source or runtime dependency on Isca.

The pre-moisture regression point is MyPlanetSim commit
`c8d5f4c71559164ea107bf96478a2c4889c6ca27` (P12.09).  At that commit the registered
Phase 5 baseline test has SHA-256
`9d9f11b32816f1bd44167965f99cefb427d5b8e74cbebdfd9cab8ff17f6593d8` and the
`phase5_umjs14_steady.cfg` preset has SHA-256
`0da750db7610c851058c311636aa107c9125ef1f5c3bbf9f7c162b36d1e2d0ea`.  Old
configuration serialization, its fingerprint, its trajectory, and the V1 checkpoint
round trip remain regression contracts.

## Decision

### State and tracer identity

Every atmospheric tracer stores dry-carrier mass times mixing ratio:

```text
Q_i = M q_i
M   = Delta p / g
q_i = kg tracer / kg dry air.
```

The dry carrier controls surface pressure, the equation of state, hydrostatic
diagnosis, heat capacities, and all fast modes.  Tracer mixing ratios are neither
renormalized nor included in `M`.  A registry supplies stable names, roles, and
nonnegativity/diffusion requirements.  Moist configurations contain exactly one
`water_vapor` role; physics receives its resolved index and never assumes index zero.

The contiguous representation is tracer-major, with each component retaining the
old cell-column ordering:

```text
Q[(tracer * cells + cell) * levels + level].
```

An omitted registry denotes the legacy one-component passive tracer.  That path keeps
the existing public field view, canonical configuration text, FrameV2 meaning, and V1
checkpoint layout.

### Dilute thermodynamics

Moisture does not alter `R_d`, `cp`, `cv`, density, Exner pressure, hydrostatic
geopotential, Richardson number, or parcel buoyancy through virtual temperature.
Those omissions are one explicit dilute-water approximation, not a consequence of
the storage convention.  Moisture affects temperature only through constant liquid
vaporization enthalpy:

```text
e_s(T) = e0 exp[(L_v/R_v)(1/T0 - 1/T)]
q_s(T,p) = (R_d/R_v) e_s(T) / p
H_m = sum_k M[k] (cp T[k] + L_v q_v[k]).
```

The constants are `T0=273.16 K`, `e0=611.657 Pa`, `L_v=2.5e6 J/kg`, and
`R_v=461.5 J/(kg K)`.  The exact specific-humidity denominator `p-e_s`, virtual
temperature, ice saturation, temperature-dependent latent heat, condensate heat
capacity, and precipitation potential/kinetic energy are deliberately excluded.
`H_m` is therefore a process-local liquid-reference enthalpy, not the existing global
dry-energy diagnostic.

Grid-scale saturation adjustment holds pressure, layer mass, and Exner pressure fixed.
For a supersaturated layer it solves

```text
q-c = q_s(T + L_v c/cp, p),  0 <= c <= q,
```

with a bracketed scalar solve.  This conserves `cp*T + L_v*q`; temperature-fixed
clipping is forbidden.  Condensate is transferred immediately to the surface of the
same cell and receives no second latent-heat increment.

### Process order and rollback

One accepted moist step has this order:

```text
saved state and accepted ledgers
  -> dynamics + radiation candidate, transporting every tracer
  -> repeat local physics substeps until their durations sum to dynamics dt
       hydrostatic and boundary-layer coefficient diagnosis
       implicit boundary layer + sensible/latent surface exchange + bucket
       dry convective preprocessing when SBM is selected
       finite-time Simple Betts--Miller adjustment
       instantaneous saturation adjustment
       route convective and grid-scale rain to land/ocean ledgers
  -> final diagnosis and invariant/budget checks
  -> atomically accept state and ledgers.
```

The default maximum physics substep is 300 s and must also not exceed one quarter of
the convection relaxation time.  A local solve failure first subdivides the physical
substep.  Failure at the supported lower bound rejects the entire dynamics candidate.
Rejected precipitation, evaporation, runoff, ocean exchange, bucket changes, and
accepted counters are rolled back with all tracers and surface temperature.  Work
counters may record rejected calls separately.  Any accepted post-RHS mutation
invalidates the same-time RHS cache.

Initialization and restart loading diagnose but do not silently remove initial
supersaturation.  The first accepted physical step owns that change.

### Surface water and ledger signs

Evaporation `E` is positive from surface to atmosphere; precipitation `P` is positive
from atmosphere to surface; sensible heat `H` is positive into the atmosphere.
Land bucket water `W` is per unit land area.  Ocean water is an unlimited reservoir
represented only by its cumulative change from the initial condition.  Runoff is the
land-capacity excess and is transferred immediately to the ocean ledger when any
ocean area exists.  In an all-land experiment it is recorded as external outflow.

The exact exchanged water flux is shared by the lowest atmospheric layer, surface
latent heat, bucket, and ocean ledgers.  Limiting evaporation without limiting latent
cooling, or repairing water afterwards with an unrecorded clamp, is forbidden.  The
accepted global identities are

```text
Delta W_atm   = integral(E_cell - P_conv - P_lsc) dA dt
Delta W_land  = integral f_land(P_conv + P_lsc - E_land - runoff) dA dt
Delta W_ocean = integral[(1-f_land)(P_conv + P_lsc - E_ocean)
                         + f_land runoff] dA dt
Delta(W_atm + W_land + W_ocean) = 0.
```

For all-land experiments the last reservoir is external outflow.  No bucket is
allocated or divided by land fraction in all-ocean cells.

### Simple Betts--Miller contract

The fixture in `tests/fixtures/phase13_sbm_columns.csv` fixes top-to-bottom column
inputs before the implementation exists.  It is interpreted with the constants and
dilute equations in this ADR.  It records branch expectations rather than copying
upstream floating-point output.

The adopted behaviour and intentional Isca differences are:

| Behaviour | MyPlanetSim Phase 13 | Isca reference |
|---|---|---|
| stored humidity | dry-air mixing ratio | input specific humidity, internally mixing ratio |
| saturation ratio | dilute `epsilon e_s/p` | exact `epsilon e_s/(p-e_s)` |
| buoyancy | temperature | virtual temperature |
| parcel integration | RK2 or better in `ln(p)` | second-order pressure-level integration |
| trigger | positive CAPE; CIN diagnostic only | positive CAPE |
| deep reference | parcel temperature and `RH_ref*q_s` | parcel temperature and RH reference |
| finite relaxation | `dt/(tau+dt)` for T and q | explicit `dt/tau`, with branch-dependent rates |
| deep conservation | uniform temperature-reference offset | humidity relaxation-time correction |
| shallow branch | fractional top giving zero integrated water change | shallower level search |
| condensate | liquid rain only | rain output; snow path inactive in this scheme |
| virtual effects | deliberately omitted | included |

Before SBM, the existing weighted dry adjustment runs once to remove only dry
superadiabatic instability.  It has a separate budget.  The driver does not invoke a
second dry adjustment.  SBM then diagnoses its parcel and reference field from the
preprocessed state.  Deep adjustment conserves `H_m` while exporting nonnegative
rain.  Shallow adjustment conserves both column water and `H_m` and exports no rain.

### Checkpoint contract

The legacy one-passive-tracer configuration continues to read and write
`dry_hydrostatic_cell_column_v1` or `dry_hydrostatic_surface_v1` byte-for-byte.
Moist/multi-component state uses a new named layout.  Its payload contains registry
count, ordered names and roles, every `Q_i`, surface temperature, land bucket, and all
cumulative water ledgers required to continue a run.  The registry metadata is also
represented in the configuration fingerprint.

Readers reject unknown layout versions, duplicate tracer names, missing or duplicate
water-vapor roles, mismatched count/name/order/role metadata, and payload-size
mismatches.  A legacy dry checkpoint is never upgraded by inventing a water field or
bucket.  Restart equivalence means the same accepted-step partition produces the same
state and cumulative ledgers as uninterrupted execution.

## Consequences

- Existing dry runs remain an executable baseline rather than merely a dry limit of a
  new code path.
- Process-local water and enthalpy closure can be tested independently of known global
  dry-energy residuals.
- Every water transfer has one owner and one sign convention.
- Strong latent heating cannot leak into RK stages or semi-implicit iterations.
- The intentional differences from Isca are visible in tests and validation reports;
  agreement with an upstream climate is not implied.
- The first model applies only while water remains dilute.  Exceeding configured
  limits is reported as an invalid input/domain condition, not hidden by clipping.

## Validation

P13.02--03 retain the frozen dry baseline and validate component permutation,
constant fields, conservative/nonnegative transport, and both checkpoint families.
P13.04 validates analytic saturation derivatives and adjustment enthalpy.  P13.05--06
validate each surface and convective transfer against the identities above and the
frozen column fixture.  P13.07 validates end-of-step placement, rollback, cache
invalidation, and restart.  P13.08--10 add column, global, timestep, sensitivity,
performance, and long-pilot evidence without weakening these kernel contracts.
