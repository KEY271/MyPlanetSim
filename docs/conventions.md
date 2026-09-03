# Scientific and numerical conventions

This document fixes conventions shared by model code, configurations, diagnostics,
and validation cases. A change to these rules requires an ADR and updated tests.

## Units

All internal calculations use SI units. Configuration keys and public C++ members carry
a unit suffix unless the quantity is dimensionless.

| Quantity | Symbol | Unit | Name suffix |
|---|---:|---:|---|
| length | `x`, `R` | m | `_m` |
| time | `t` | s | `_s` |
| velocity | `u` | m s⁻¹ | `_m_s` |
| acceleration | `g` | m s⁻² | `_m_s2` |
| angular velocity | `Omega` | rad s⁻¹ | `_rad_s` |
| pressure | `p` | Pa | `_pa` |
| temperature | `T` | K | `_k` |
| specific energy/heat capacity | — | J kg⁻¹ / J kg⁻¹ K⁻¹ | `_j_kg` / `_j_kg_k` |

Angles passed through APIs are radians. Degrees are allowed only in human-facing input
or output and must carry `_deg`. Radian values are dimensionless in C++, so naming is
the required guard against accidental mixing.

## Floating-point representation

- `mps::Real` is IEEE 754 binary64 (`double`).
- Prognostic calculations do not silently replace NaN or infinity with defaults.
- Exact equality is used for identities that must be preserved structurally, such as
  copying a constant field or checkpoint round-trip.
- Analytical comparisons state an absolute and/or relative tolerance at the call site.
- Signed zero is not physically distinguished, but diagnostics may normalize it when
  producing canonical text.

For scalar comparison, the intended predicate is

```text
abs(actual - expected) <= absolute_tolerance
  + relative_tolerance * abs(expected)
```

Tolerance is selected from truncation error and conditioning, not increased merely to
make a regression pass.

## Planet-fixed Cartesian frame

MyPlanetSim uses a right-handed, planet-fixed Cartesian frame:

- origin at the planet center;
- `+z` along the declared north rotation axis;
- `+x` at longitude 0 on the equator;
- `+y` at longitude `+pi/2` on the equator;
- positive longitude increases eastward;
- latitude is geocentric and positive northward.

For positive `rotation_rate_rad_s`, the rotation vector is `Omega = (0, 0, +Omega)`.
Negative values represent a planet rotating in the opposite direction. The Coriolis
acceleration in the rotating frame is `-2 Omega cross u`. The traditional horizontal
momentum equations therefore use `f = 2 Omega sin(latitude)` with `f > 0` in the
northern hemisphere for positive rotation.

Local horizontal basis vectors are east then north. A physical vector crossing a
cubed-sphere panel boundary will be transformed through the Cartesian frame; panel-local
components are never copied directly.

## Array and grid indices

- Horizontal interior indices are `i = 0..nx-1`, `j = 0..ny-1`.
- `i` is the fastest-varying index in row-major storage.
- A halo of width `h` extends logical indices to `-h..nx+h-1` and
  `-h..ny+h-1`.
- Global cubed-sphere panel order and edge orientation will be fixed by the Phase 3
  geometry ADR; code must not infer it from enum integer values before then.

`Field2D` owns storage but does not own coordinates or cell metrics. Algorithms that
integrate a field receive weights explicitly.

## Vertical indices

Vertical index `k` increases downward. With `nz` layers:

```text
half levels: k+1/2 = 0 .. nz
full levels: k     = 0 .. nz-1

top    half 0 ─────────────
             full 0
         half 1 ────────────
             full 1
                 ...
surface half nz ────────────
```

Pressure increases downward. Hybrid interface pressure is

```text
p_half[k] = A_pa[k] + B[k] * surface_pressure_pa
```

where `A_pa` is in Pa and `B` is dimensionless. Interfaces must be monotone and every
layer mass `delta_p/g` must be positive. Full-level horizontal wind, potential
temperature, and tracers use `k`; pressure and vertical mass flux use interfaces.

## Time and tendency signs

- Simulation time and time step are seconds and increase forward.
- A positive flux is in the positive coordinate direction of its oriented face.
- A cell's conservative tendency is inflow minus outflow divided by its measure.
- Heating is positive when it increases temperature or potential temperature.
- Surface geopotential and topographic height are positive outward from the reference
  sphere.

The final time step may be shortened to land exactly on `end_time_s`; the solver must
never advance past it and then interpolate backward.

## Error norms and conservation drift

For cell value error `e_i` and positive measure `w_i`, normalized norms are

```text
L1   = sum(w_i * abs(e_i)) / sum(w_i)
L2   = sqrt(sum(w_i * e_i^2) / sum(w_i))
Linf = max(abs(e_i))
```

Relative drift of a conserved scalar `Q` is

```text
drift(Q, t) = (Q(t) - Q(0)) / max(abs(Q(0)), Q_scale)
```

The validation case defines `Q_scale`; it prevents division by zero but must not hide a
large absolute residual. Reports include both absolute change and relative drift.

Observed order between resolutions or time steps `h` and `h/2` is

```text
p_observed = log(error(h) / error(h/2)) / log(2)
```

Norms use positive physical measures. Unweighted norms are allowed only when the grid is
uniform and the test states this explicitly.

## Naming map

| Mathematical quantity | C++/configuration name |
|---|---|
| planet radius `R` | `radius_m` |
| rotation rate `Omega` | `rotation_rate_rad_s` |
| gravity `g` | `gravity_m_s2` |
| dry gas constant `R_d` | `gas_constant_j_kg_k` |
| heat capacity `c_p` | `heat_capacity_cp_j_kg_k` |
| reference pressure `p_0` | `reference_pressure_pa` |
| surface pressure `p_s` | `surface_pressure_pa` |
| potential temperature `theta` | `potential_temperature_k` |
| layer pressure thickness `Delta p` | `pressure_thickness_pa` |

ASCII identifiers are used in source and data files. Unicode mathematical notation is
reserved for prose and figures.

## Open decisions

The following are intentionally not fixed in Phase 0:

- cubed-sphere panel numbering and local coordinate mapping;
- prognostic momentum representation and horizontal staggering;
- exact hybrid `A/B` profile and model top pressure;
- large field file format and parallel I/O library;
- deterministic reduction strategy across multiple threads or MPI ranks.

Each decision is made before its first implementation and recorded in an ADR.
