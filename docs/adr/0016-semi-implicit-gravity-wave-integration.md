# ADR 0016: Reference-linear semi-implicit gravity-wave integration

**Status:** accepted for staged implementation

## Context

The Phase 9 dry-hydrostatic core advances every tendency with SSP-RK3. At `N=12`,
`K=20`, the initial Held--Suarez state has a fast horizontal stability limit of
240.578 s, so the production preset uses 200 s. The same state has no material or
vertical Courant restriction at rest. After 30 model days the measured material-
advection and vertical-transport limits are still 4,592 s and 16,592 s, while the
fast-wave limit is 225.075 s. The immediate stiffness is therefore the horizontal
Lamb/gravity-wave block, not transport.

Increasing only `run.time_step_s` is invalid: the explicit method would cross its
registered wave CFL. Conversely, changing the spatial operator or adding a mass fixer
would make it difficult to attribute differences from the explicit reference.

## Decision

### Equation split

Keep the prognostic layout

```text
U = (ps, M*u, M*theta, M*q, optional Ts)
```

and the existing full nonlinear tendency `F(U)`. Construct a time-independent resting,
horizontally uniform hydrostatic reference `U_ref` and a tangent-linear fast operator
`L_ref`. Define the slow residual by subtraction using the same discrete paths:

```text
N(U) = F(U) - L_ref(U - U_ref)
```

The fast operator contains surface-pressure continuity, horizontal divergence,
temperature/hydrostatic diagnosis, and pressure acceleration. It excludes material
upwinding, limiters, tracers, Coriolis, diffusion, vertical transport, surface storage,
and physical parameterizations. Tests must establish both `L_ref(0)=0` and the
component-wise add--subtract identity before the implicit driver is enabled.

### Time discretization

Use one-step iterative Crank--Nicolson with implicit weight `alpha`:

```text
R(U_l) = U_l - U_n
         - dt * ((1 - alpha) * F(U_n) + alpha * F(U_l))

(I - alpha * dt * L_ref) deltaU = -R(U_l)
U_(l+1) = U_l + deltaU
```

The accuracy reference uses `alpha = 0.5`. Production configurations may use values up
to 0.55 only when their damping and phase effects are recorded. The method carries no
past tendency, so the existing instantaneous checkpoint state remains sufficient for a
deterministic restart.

### Vertical modes and horizontal solve

Build the hydrostatic reference column once at driver construction. Diagonalize its
vertical structure matrix in deterministic descending phase-speed order. A mode is
implicit when its per-cell wave Courant exceeds the configured threshold. If more modes
are needed than the configured cap, reject the run rather than silently leaving an
unstable mode explicit.

Eliminate the thermodynamic and pressure variables per mode to obtain matrix-free 2-D
Helmholtz systems. The existing least-squares gradient and finite-volume divergence are
not assumed symmetric, so use restarted GMRES rather than CG. The first preconditioner
is area-scaled Jacobi; later coarse-grid work is driven by measured iteration counts.
Both relative and absolute residuals must pass. A non-finite residual or iteration-limit
exit is a failed solve.

Implicit `ps` and `M*theta` corrections are applied through shared oriented edge-flux
scatter. No global-mean repair, pressure clipping, temperature clipping, or tracer
clipping is allowed.

### Step selection and failure

`run.time_step_s` is a requested maximum, not a promise. A semi-implicit step must obey
the remaining material-advection, vertical-transport, diffusion, surface-reservoir, and
invariant constraints. A failed linear solve, nonlinear residual, stage constraint, or
state invariant causes a retry from `U_n` at half the step. Crossing
`semi_implicit.minimum_time_step_s` is a diagnosed hard failure.

The registered delivery target remains 1,800 s, with median accepted Held--Suarez step
at least 1,200 s. This is not an artificial ceiling. After the 1,800/900/450 s accuracy
gate passes, each validation family also probes 2,400, 3,600, 5,400, and 7,200 s in order,
stopping after the first failed stability, accuracy, invariant, or cost gate. The last
passing value is reported separately for:

- linear flat waves;
- balanced terrain and mountain waves;
- dry baroclinic dynamics;
- Held--Suarez and each later physical configuration.

Thus a dynamics-only experiment may register a step longer than 30 minutes while a
physics-coupled preset retains a smaller operational maximum. Stability alone is not
sufficient: phase and amplitude error, climate statistics, retry frequency, and solver
cost are part of the limit.

### Configuration contract

The accepted conditional keys are:

```ini
dry_hydrostatic.time_integrator = semi_implicit
dry_hydrostatic.advective_cfl = 0.45
semi_implicit.reference_surface_pressure_pa = 100000
semi_implicit.reference_temperature_k = 300
semi_implicit.implicit_weight = 0.5
semi_implicit.wave_cfl_threshold = 0.45
semi_implicit.maximum_implicit_modes = 5
semi_implicit.nonlinear_relative_tolerance = 1e-8
semi_implicit.nonlinear_maximum_iterations = 4
semi_implicit.linear_relative_tolerance = 1e-8
semi_implicit.linear_absolute_tolerance = 1e-12
semi_implicit.linear_maximum_iterations = 40
semi_implicit.gmres_restart = 20
semi_implicit.minimum_time_step_s = 60
```

Absent keys preserve the current explicit SSP-RK3 path and its canonical text,
fingerprint, checkpoint, Frame, and numerical result. Semi-implicit keys are rejected for
other experiment kinds and participate in the fingerprint when active.

## Consequences

- The explicit solver remains the regression oracle and emergency operational path.
- Fast-wave Courant is reported but does not limit a semi-implicit step; material,
  vertical, diffusion, surface, and physics constraints remain authoritative.
- Thirty minutes is a required capability, not a hard-coded maximum.
- Different physical packages may have different registered operational limits.
- The first implementation is optimized for the current small serial cubed-sphere
  problem. MPI/GPU solvers and semi-Lagrangian transport remain later work.
- Terrain failure cannot be waived as a flat-only success.

## Validation

The staged gates are recorded in [Phase 10 validation](../validation/phase-10.md). They
cover algebraic GMRES failures, zero-allocation reuse, linear shallow-water waves,
reference-column modes, add--subtract identity and conservation, flat and terrain dry
dynamics, restart/retry determinism, forced climate, and measured model-day cost.
