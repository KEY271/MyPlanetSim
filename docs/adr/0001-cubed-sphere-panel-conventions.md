# ADR 0001: cubed-sphere panel conventions

Status: accepted for Phase 1.

## Mapping and panel order

The global panel order is `PX, PY, NX, NY, PZ, NZ`.  Every panel uses
`-pi/4 <= alpha,beta <= pi/4` and

```text
r = normalize(n + tan(alpha) e_alpha + tan(beta) e_beta).
```

The bases are:

| panel | n | e_alpha | e_beta |
|---|---|---|---|
| PX | (1,0,0) | (0,1,0) | (0,0,1) |
| PY | (0,1,0) | (-1,0,0) | (0,0,1) |
| NX | (-1,0,0) | (0,-1,0) | (0,0,1) |
| NY | (0,-1,0) | (1,0,0) | (0,0,1) |
| PZ | (0,0,1) | (1,0,0) | (0,1,0) |
| NZ | (0,0,-1) | (1,0,0) | (0,-1,0) |

All bases satisfy `e_alpha cross e_beta = n`. Cells are flattened panel-major,
then `j`, then `i`. Cell corners follow `SW, SE, NE, NW`; this traversal is
positive on the outward-oriented sphere.

The inverse map selects the largest absolute Cartesian component. Exact ties use
Cartesian-axis order `x, y, z`, then the sign of that component. This canonical
tie-break is for coordinate queries only; connectivity is assembled from shared
vertices and never depends on the inverse-map tie-break.

## Edge names and connectivity

West/east edges have increasing-beta natural index; south/north edges have
increasing-alpha natural index. `reversed` means index `k` maps to `N-1-k`.

| panel edge | neighbor edge | index |
|---|---|---|
| PX W | NY E | same |
| PX E | PY W | same |
| PX S | NZ E | reversed |
| PX N | PZ E | same |
| PY W | PX E | same |
| PY E | NX W | same |
| PY S | NZ S | reversed |
| PY N | PZ N | reversed |
| NX W | PY E | same |
| NX E | NY W | same |
| NX S | NZ W | same |
| NX N | PZ W | reversed |
| NY W | NX E | same |
| NY E | PX W | same |
| NY S | NZ N | same |
| NY N | PZ S | same |
| PZ W | NX N | reversed |
| PZ E | PX N | same |
| PZ S | NY N | same |
| PZ N | PY N | reversed |
| NZ W | NX S | same |
| NZ E | PX S | reversed |
| NZ S | PY S | reversed |
| NZ N | NY S | same |

A unique edge stores the first incident cell as `left`. Its vertices follow that
cell's positive boundary traversal. Positive flux is outward from `left` and
inward to `right`; the accumulator applies exactly the same number with opposite
signs. The tangent-plane normal is `edge_tangent cross edge_center`.
