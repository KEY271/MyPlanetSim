# ADR 0025: Keep the prognostic cell-column layout behind field views

## Decision

The canonical prognostic and checkpoint order remains
`[component][cell][level]`. `Field3DView`, `ColumnView`, and the component axis used
as a tracer view now make bounds and strides explicit without changing stored values.

A blocked `[component][cell_block][level][lane]` view is available for temporary
kernel packing with lane widths 4, 8, and 16. It is not adopted as the prognostic
layout: horizontal traversal can benefit, while column physics becomes strided and
packing, padding, restart, and named-tracer conversion costs must be included in the
decision.

## Consequences

- Existing restart files and visualization order remain unchanged.
- Production tracer preparation uses the view abstraction while retaining the same
  arithmetic and results.
- `benchmark_field_layout` provides a foreground comparison of horizontal traversal,
  column traversal, padding, and one-time packing. A full model-day and RSS comparison
  is still required before expanding the blocked layout beyond temporary workspaces.
