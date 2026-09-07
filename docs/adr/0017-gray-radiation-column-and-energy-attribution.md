# ADR 0017: Gray-radiation column fluxes and energy attribution

**Status:** accepted for staged implementation

## Context

The Phase 8 surface scheme absorbs stellar radiation and emits longwave radiation
directly to space. It cannot represent atmospheric absorption, emission, or downward
longwave flux. Phase 11 adds those processes without changing the dry prognostic
layout or hiding the known dry-core energy defect documented by ADR 0015.

The transport kernel must also remain testable without a grid, dynamics state, or I/O.
Its boundary conditions and signs therefore have to be fixed before it is connected to
the three-dimensional driver.

## Decision

### Column coordinates, opacity, and boundaries

Half levels are ordered from model top to surface, `j=0..K`; layer `k` lies between
`j=k` and `j=k+1`. Downward and upward flux arrays are both stored as non-negative
directional magnitudes. The net flux used for budgets is upward positive:

```text
F[j] = LW_up[j] - LW_down[j] + SW_up[j] - SW_down[j]
Q_rad[k] = F[k+1] - F[k]
```

`Q_rad>0` heats the layer. The model top is a finite-pressure radiative boundary with
no incident longwave radiation above it. Incoming direct stellar flux is
`S max(0, r_hat dot s_hat)`. No unmodelled atmosphere is placed above `p_top`.

Layer optical depths are physical vertical optical depths, separate from the angular
diffusivity factors:

```text
Delta tau_SW = kappa_SW (p_bottom-p_top)/g
Delta tau_LW = kappa_ref p_ref/(g n)
               [(p_bottom/p_ref)^n-(p_top/p_ref)^n]
```

They use absolute pressure rather than `p/ps`; terrain and gravity may therefore alter
total column opacity. The implementation evaluates thin pressure-power differences in
a form that does not subtract nearly equal powers.

The shortwave surface boundary is diffuse reflection,
`SW_up[K]=alpha SW_down[K]`. The longwave surface boundary obeys Kirchhoff's law,

```text
LW_up[K] = epsilon sigma Ts^4 + (1-epsilon) LW_down[K].
```

Night and the exact terminator have zero shortwave flux and never divide by `mu0`.
There is no artificial zenith-angle floor, temperature clipping, or energy fixer.

### Source conversion and local conservation

For layer mass per area `M=(p_bottom-p_top)/g`, the column convergence enters the
existing potential-temperature-mass prognostic through the already diagnosed full-
level Exner function:

```text
d(M theta)/dt = Q_rad/(cp Pi_full).
```

At fixed mass this gives `dT/dt=Q_rad/(cp M)`. Radiation does not directly change
surface pressure, dry mass, tracer mass, or momentum. Sensible heat `H`, upward
positive, enters the bottom layer as `H/(cp Pi_full)` and the surface as `-H`.
The surface reservoir uses the same bottom-interface flux as the atmosphere:

```text
C_s dTs/dt = -F[K] + Q_internal - H.
```

Consequently the transport kernel must satisfy, separately for shortwave and
longwave and to roundoff,

```text
sum_k Q_rad[k] = F[K] - F[0].
```

No independently recomputed surface or TOA flux may be used to update a reservoir.

### Relation to the diagnosed dry total energy

The diagnostic dry energy is

```text
E_dry = sum_k A M_k (|u_k|^2/2 + cv T_k + Phi_k).
```

For a radiation-only direction, pressure, mass, and velocity are fixed, so Exner is
fixed. Write `X_k=d(M theta)_k/dt` and `theta_dot_k=X_k/M_k`. The exact directional
derivative of the implemented diagnostic is

```text
dE_dry/dt = A sum_k [cv Pi_full,k X_k + M_k Phi_dot_full,k],

Phi_dot_half,K = 0,
Phi_dot_full,k = Phi_dot_half,k+1
                 + cp theta_dot_k (Pi_half,k+1-Pi_full,k),
Phi_dot_half,k = Phi_dot_half,k+1
                 + cp theta_dot_k (Pi_half,k+1-Pi_half,k).
```

The first term is `(cv/cp) sum Q_rad`, not `sum Q_rad`; hydrostatic re-diagnosis
supplies the second term. Because the top pressure is finite and the model uses its
existing full-level quadrature, their sum is not assumed algebraically equal to the
radiative boundary flux. Phase 11 reports these quantities separately:

1. interface-flux convergence and its telescoping residual;
2. the temperature/internal-energy direction;
3. the hydrostatic-geopotential direction, including the finite-top contribution;
4. the centered finite-difference derivative of the implemented total-energy
   diagnostic; and
5. the remaining dynamics and time-integration residuals.

The analytic direction above must match a centered perturbation ladder of the actual
state diagnosis before the radiation source is accepted. Sensitivity to a lower
`p_top` is measured separately; it is not corrected by adding boundary work whose
discrete counterpart is absent from the state.

### Compatibility baseline

`physics.kind = gray_radiation` is a new, exclusive physics path. It owns atmospheric
radiation, the existing surface-style sensible exchange, and lower-layer drag. It
does not call `surface_energy_balance` in addition. The meanings and canonical
serialization of `none`, `held_suarez`, `planetary_newtonian`, and
`surface_energy_balance` remain unchanged.

Before the first implementation commit, the Phase 8/config compatibility subset was
run at commit `77ba092`:

```text
unit.experiment_config
unit.run_metadata
unit.checkpoint
unit.frame
unit.surface_energy_balance
phase8.planetary_core
```

All six tests passed. Later stages compare transparent gray radiation against the
Phase 8 tendency using the same explicit time step; a new stability bound is not
allowed to disguise a tendency difference.

## Consequences

- A standalone `O(K)` column kernel is the sole owner of radiative interface fluxes.
- Atmospheric and surface tendencies are derived from those stored fluxes.
- Local radiative conservation can pass even while the dry-core total-energy
  diagnostic retains a separately attributable defect.
- A finite model top, gray absorption, and the diffusivity factors are explicit model
  choices whose sensitivity must be reported rather than tuned silently.
- Restart payloads need no radiative prognostic state; accumulated diagnostics still
  need an explicit continuity contract when they are added.

## Validation

Phase 11 gates test optical-depth limits, analytic longwave and Beer--Lambert fluxes,
shortwave and longwave telescoping identities, non-black surface reflection, the
source directional derivative, source-only time convergence, driver coupling,
restart, and the registered global comparison matrix. Quantitative results are kept
in `docs/validation/phase-11.md` as stages are completed.
