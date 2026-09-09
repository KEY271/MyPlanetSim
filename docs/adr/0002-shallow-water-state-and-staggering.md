# ADR 0002: shallow-water state and staggering

**Status:** retired in Phase 15

The implementation, dedicated configuration, and tests described here were removed in
P15.07. This ADR is retained only as the Phase 2 design record; it does not describe a
currently selectable model.

## Context

Phase 2 adds the rotating, one-layer shallow-water equations on the existing
equiangular gnomonic cubed sphere.  Panel-local vector components are not a
single-valued representation at seams, so all stored and exchanged vectors use
the planet-fixed Cartesian frame.

## Decision

The reference scheme is a collocated, cell-centred finite-volume method.  Cell
averages store layer depth `h` in metres and Cartesian momentum `m = h u` in
`m^2 s^-1`.  At a cell centre with outward unit normal `k`, momentum satisfies
`k dot m = 0`; velocity is diagnosed as `u = m / h`.  Depth must be finite and
strictly greater than the configured `depth_floor_m` at every Runge--Kutta
stage boundary.

The continuous vector-invariant equations used to define signs are

```text
partial_t h + div(h u) = 0
partial_t u + (zeta + f) k cross u
  + grad(g h + 0.5 |u|^2) = D_u
zeta = k dot curl(u),       f = 2 Omega dot k
Omega = (0, 0, planet.rotation_rate_rad_s)
```

The reference implementation evaluates the equivalent conservative Cartesian
momentum flux.  For an oriented unique edge, `n` points out of the left cell,
`t = k_edge cross n`, and primitive states are expressed as `(h,u_n,u_t)`.
The physical normal flux is

```text
F_h       = h u_n
F_m       = h u_n u + 0.5 g h^2 n
lambda    = |u_n| + sqrt(g h)
F_Rusanov = 0.5 (F(q_L) + F(q_R)) - 0.5 max(lambda_L,lambda_R) (q_R-q_L)
```

The vector jump is compared in the shared Cartesian tangent plane at the edge.
Each edge flux is evaluated once, multiplied by edge length, subtracted from
the left cell and added to the right cell.  This makes global mass cancellation
independent of summation order.  The Cartesian tensor-flux divergence contains
the spherical curvature contribution: after division by cell area its momentum
tendency is projected with `P_k(v)=v-(v dot k)k`.  Pressure therefore remains in
the edge flux and is not added a second time as a panel-coordinate source.

Coriolis is a cell source,

```text
S_m = -2 h P_k(Omega cross u) = -2 P_k(Omega cross m).
```

Flux, Coriolis and diffusion tendencies are retained separately for budgets.
After they are summed, each SSP-RK3 stage projects updated momentum onto the
cell-centre tangent plane.  Radial momentum immediately before and after this
projection is diagnosed.

The comparison scheme has C-grid staggering: `h` on primal cells, oriented
normal velocity or mass flux on primal edges, and circulation, relative
vorticity and potential vorticity on dual cells around primal vertices.  It
shares experiment setup, time control and diagnostics with the reference
scheme, but has a distinct state/checkpoint layout.  Its PV flux and Hodge
operators are separate from the reference Rusanov flux; their conservation
claims will be decided from the P2.18 comparison rather than assumed here.

## Consequences

- Seam exchange never copies panel-local vector components.
- Mass conservation follows from unique-edge scatter; momentum is tangent but
  is not claimed to be a globally conserved Cartesian vector on a sphere.
- The reference scheme is deliberately diffusive and is the independent
  comparator for the compatible scheme.
- Orography and its pressure/topographic source are outside Phase 2.
