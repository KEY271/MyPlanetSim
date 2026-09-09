# ADR 0027: Tile horizontal edge work without changing the stage state

## Decision

Each cubed-sphere panel is partitioned into logical square tiles. The production width
is 16 cells; widths 8, 16, and 32 are available through the benchmark constructor.
The halo radius is derived from the widest current stencil and is two cells because
explicit diffusion applies a two-stage Laplacian.

Tile halos are built by graph distance over the grid edge cache, so panel seams and
cube corners use the same global cell and edge identities as panel interiors. Every
edge is assigned once to the lower numbered incident tile. Edge reconstruction and
flux therefore run once per edge, while each tile gathers only its disjoint interior
cells in edge-ID order.

## Consequences

- All kernels read the same immutable stage state. No tile can observe a neighboring
  tile's partially updated tendency or prognostic value.
- The global prepared-gradient and compact edge-flux buffers remain. Removing them
  would require a boundary exchange protocol whose cost and benefit have not yet been
  established; this commit does not force that expansion.
- Tile metadata records interior, halo, and owned-edge indices once per driver and is
  reported by RHS profiling.
- `benchmark_tiles` compares widths and rejects any RHS bit-pattern difference before
  reporting timing. High-resolution/repeated measurements remain foreground work.
