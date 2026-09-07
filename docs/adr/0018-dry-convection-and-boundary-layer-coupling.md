# ADR 0018: Dry convection and boundary-layer coupling

**Status:** accepted for staged implementation

## Context

Phase 11 deliberately leaves dry convective instability and vertical turbulent
mixing unresolved.  Its gray-radiation adapter also owns a legacy linear sensible
heat exchange and a lower-atmosphere Rayleigh drag.  Adding a boundary layer without
changing those ownership rules would count surface heat and momentum exchange twice.

The dry state stores layer mass times potential temperature and Cartesian horizontal
momentum.  Hydrostatic geopotential is diagnosed from that state, while the time
integrators evaluate right-hand sides more than once per accepted step.  Instantaneous
adjustment and implicit vertical diffusion therefore cannot be hidden inside a right-
hand-side evaluation.

## Decision

### Process ownership and ordering

Convection and boundary-layer schemes are configuration axes independent of
`physics.kind`.  Their default is `none`, preserving the canonical serialization,
fingerprint, and trajectory of every existing configuration.  The first implementation
supports either scheme only with the dry hydrostatic gray-radiation path.

An accepted coupled step is ordered as follows:

```text
saved state n
  -> dynamics + radiation candidate
  -> implicit boundary-layer/surface step, when enabled
  -> dry convective adjustment, when enabled
  -> hydrostatic diagnosis and invariant checks
  -> accept state and all process budgets together
```

A failure after the dynamics candidate rejects the whole candidate.  Retry starts
from the saved state and accepted budgets exclude rejected work.  Neither process is
run during initialization, restart loading, Runge--Kutta stages, or semi-implicit
iterations.  Any accepted post-RHS state change invalidates same-time RHS caches.

When the boundary layer is disabled, gray radiation retains its legacy sensible heat
and Rayleigh drag.  When the boundary layer is enabled, it owns all surface sensible
heat, surface stress, and lower-atmosphere friction; gray radiation owns only
radiative convergence, radiative surface storage, and internal surface heat.  A
nonzero legacy `surface.air_exchange_coefficient_w_m2_k` is invalid in this mode.
Held--Suarez benchmark drag is not reinterpreted.

### Dry convective adjustment

Vertical indices increase downward.  During adjustment, pressure, layer mass `M`,
and Exner functions are fixed.  A column is stable when
`theta[k] >= theta[k+1]`.  Adjacent violating blocks are merged with a weighted PAVA,
using

```text
w[k] = cp M[k] Pi_full[k]
theta_block = sum(w[k] theta[k]) / sum(w[k]).
```

Every member of a merged block receives `theta_block`.  Thus the kernel is `O(K)`,
stable to the roundoff-scale comparison tolerance, and idempotent.  It preserves the
fixed-mass column enthalpy exactly to floating-point reduction error:

```text
E_h = sum_k cp M[k] T[k] = sum_k cp M[k] Pi_full[k] theta[k].
```

It does not alter mass, pressure, momentum, tracers, or surface temperature.  The
operation is a finite increment, not a tendency whose result depends on the caller's
time step.

This invariant is distinct from the diagnosed dry energy.  For
`Delta theta[k] = theta_after[k]-theta_before[k]`, hydrostatic re-diagnosis gives the
exact finite attribution because the implemented geopotential diagnosis is linear in
theta:

```text
Delta Phi_half[K] = 0
Delta Phi_full[k] = Delta Phi_half[k+1]
  + cp Delta theta[k] (Pi_half[k+1]-Pi_full[k])
Delta Phi_half[k] = Delta Phi_half[k+1]
  + cp Delta theta[k] (Pi_half[k+1]-Pi_half[k])

Delta E_dry / A = sum_k M[k]
  [cv Pi_full[k] Delta theta[k] + Delta Phi_full[k]].
```

Tests compare this expression with a before/after diagnosis.  The fixed lower
boundary has `Delta Phi_half[K]=0`, so it contributes no moving-boundary work.  The
finite-pressure model top and the model's full-level quadrature are already represented
by the recurrence; no synthetic top or surface work and no uniform energy fixer are
added.  `Delta E_h`, `Delta E_dry`, and their attribution residual remain separate
diagnostics.

### Boundary-layer fluxes and implicit solve

Height is measured from local surface geopotential.  Internal interface heat,
momentum, and tracer fluxes are defined once, upward positive, by

```text
F_H[j] = -rho[j] cp Pi_half[j] K_h[j]
         (theta_upper-theta_lower)/dz[j]
F_m[j] = -rho[j] K_m[j] (u_upper-u_lower)/dz[j]
F_q[j] = -rho[j] K_q[j] (q_upper-q_lower)/dz[j].
```

The top flux is zero.  Surface fluxes are the bulk sensible heat `H`, the stress on
the atmosphere `tau_s`, and zero tracer flux.  The same interface value is applied
with opposite signs to its two reservoirs:

```text
cp Pi[k] d(M theta)[k]/dt = F[k+1]-F[k] + D[k]
d(M u)[k]/dt              = F_m[k+1]-F_m[k]
d(M q)[k]/dt              = F_q[k+1]-F_q[k]
C_s dTs/dt                = -H.
```

All geometry, density, stability, and transfer coefficients are diagnosed at the
start of the boundary-layer substep and held fixed.  Backward Euler then solves the
atmospheric heat reservoirs and `Ts` as one tridiagonal system and solves each
momentum or tracer component as a closed atmospheric system.  Positive capacities
and nonnegative conductances produce an M-matrix.  Consequently a passive tracer
obeys the discrete maximum principle, while shared fluxes give column conservation:

```text
sum(cp M Delta T) + C_s Delta Ts = sum(D dt)
sum(M Delta q) = 0
sum(M Delta u) = dt tau_s.
```

Uniform potential temperature, paired with `Ts/Pi_s` at the surface, is the zero heat-
flux state.  This is an enthalpy-flux discretization; it does not additionally claim
unweighted potential-temperature conservation.

### Stress dissipation

For fixed interface momentum conductance `G`, fixed surface drag coefficient `beta`,
and the backward-Euler result `u+`, the kinetic-energy decrease is partitioned by the
exact discrete identity

```text
KE(n) - KE(+)
 = dt sum_interfaces G |u_upper+ - u_lower+|^2
 + dt beta |u_bottom+|^2
 + sum_layers M |u+ - u(n)|^2 / 2.
```

The first two terms are physical shear and surface-drag dissipation for the frozen
operator.  The last is backward-Euler time-discretization dissipation.  All three are
reported separately and returned as atmospheric heat: interface physical terms are
split equally between adjacent layers, surface physical dissipation goes to the
bottom layer, and numerical dissipation stays in the layer where it is incurred.
The resulting `D[k]` is included in the heat solve.  This closes the local mixing
budget only; it never compensates dynamics or radiation defects.

The initial coupled scheme is first-order splitting with backward Euler.  Stability
of the implicit matrix is not used as an accuracy claim; time-step refinement measures
coefficient freezing and splitting errors.

## Consequences

- Standalone column kernels remain independent of the grid, driver, and I/O.
- Convection can be validated before boundary-layer parameterization is introduced.
- Surface exchange has one owner for every supported configuration.
- Restart needs no new prognostic field, but configuration fingerprints and accepted
  accumulated budgets include the selected schemes and coefficients.
- Process-local enthalpy closure and the existing whole-system dry-energy residual are
  both visible and cannot be mistaken for one another.

## Validation

P12.02 validates analytic and multi-block adjustment, stability, enthalpy conservation,
idempotence, and dry-energy attribution.  P12.03 validates end-of-step placement,
rollback, cache invalidation, restart, and the disabled-scheme baseline.  P12.04--05
validate implicit diffusion, the maximum principle, shared-flux budgets, dissipation,
bulk exchange, stability limits, and fractional surfaces.  Later stages validate the
fully coupled column and global matrices described in the Phase 12 plan.
