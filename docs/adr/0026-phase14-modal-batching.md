# ADR 0026: Batch vertical modal transforms across cells

## Decision

The semi-implicit momentum projection and reconstruction process eight independent
cell lanes together. A cell block is packed as `[level][lane]` for x, y, and z, then
the existing mode-major coefficient matrices are applied as multiple right-hand
sides. The level or mode accumulation order within each cell is unchanged.

Only selected implicit modes are projected and accumulated. Unselected modes remain
the identity contribution handled by the nonlinear iteration, as before. The
eigendecomposition, mode ordering, Helmholtz solve, tolerances, and global fast
operator are not changed.

## Consequences

- Tail blocks are masked by their logical width; no vertical-level padding is added.
- Pack storage is retained in the semi-implicit workspace and reused.
- The compiler may vectorize the independent lane loop; OpenMP builds make that
  independence explicit with `omp simd` without enabling global fast-math.
- Scalar projection/reconstruction remain as reference fixtures. Performance is
  evaluated with the existing foreground semi-implicit benchmark, including solver
  iterations rather than timing the matrix multiply alone.
