# ADR 0008: fixed orography and terrain lower boundary

**Status:** proposed for Phase 6

## Context

Phase 5 fixes a flat, impermeable lower boundary and diagnoses each hydrostatic
column from a scalar surface pressure. Phase 6 must admit horizontally varying
surface geopotential without creating a second atmospheric state, changing the
hybrid-coordinate continuity recurrence, or making terrain part of the browser
protocol.

The existing dry pressure-gradient form is already the continuous generalized
pressure-coordinate identity

```text
G[k] = grad_eta(Phi_full[k]) + alpha[k] * grad_eta(p_full[k])
alpha = Rd*T/p
```

and hydrostatic integration already accepts a surface-geopotential lower boundary.
The smallest terrain extension is therefore to supply a fixed cell field to that
integration and test the existing identity on sloping coordinate surfaces. A second,
terrain-special pressure solver would obscure whether the current discretization
converges and would reopen ADR 0006 before measurements require it.

Phase 6 also reuses Williamson shallow-water test 5. Its terrain is fixed data used by
the momentum source, not a new shallow-water prognostic variable.

## Decision

### One immutable horizontal field

For `C = 6*N*N` horizontal cells, construct one run-owned field

```text
SurfaceOrography
  surface_geopotential_m2_s2[C]
  source_fingerprint
```

in the existing panel-major cell order. `surface_height_m` is derived as `Phi_s/g`
for input, output, and diagnostics; the solver passes geopotential so gravity is not
applied in more than one dynamics component.

The field is immutable during a run. It is finite and has exactly `C` values. Phase 6
benchmark and imported fields are nonnegative, but the reusable value contract only
requires finiteness so that a later reference datum may lie above or below zero.
Surface pressure and every hybrid layer must still satisfy the existing pressure
bounds and positivity checks after the terrain-dependent initial state is built.

Orography is fixed configuration, not prognostic state. It is not added to
`DryHydrostaticState`, `ShallowWaterState`, or either checkpoint payload. Restarts
reconstruct it from the canonical configuration and reject a changed configuration or
input fingerprint. The Phase 4 scalar
`vertical.surface_geopotential_m2_s2` remains a column-only lower-boundary fixture; it
is not overloaded as a global terrain amplitude.

### Configuration and named fields

`OrographyParameters` is additive and defaults to `flat`. Existing Phase 0--5
configuration files therefore remain valid and their canonical text remains
unchanged. Non-flat configurations write only these keys:

```text
orography.kind = dcmip_2_0_0 | williamson5 | linear_bell | jw06 | latlon_csv
orography.input_file = <config-relative path>              # latlon_csv only
orography.input_fingerprint_fnv1a64 = <16 lowercase hex>   # latlon_csv only
orography.smoothing_passes = <nonnegative integer>         # latlon_csv only
```

Published benchmark constants remain inside their named initializer. They are not
expanded into independent mountain-height, width, centre, and taper knobs. A named
test case validates its required orography kind, and published analytic profiles
require zero smoothing passes.

The file fingerprint uses the same deterministic FNV-1a 64-bit byte hash as the
current configuration fingerprint. It detects accidental source changes but is not a
security primitive. The loader verifies it before interpolation; the expected value
is part of canonical configuration and therefore of checkpoint and frame identity.
Relative paths are resolved against the configuration file, while canonical text
retains the portable relative spelling.

### Dry hydrostatic lower boundary

`DryHydrostaticDriver` owns the orography field. Diagnosis passes `Phi_s[c]` to the
existing hydrostatic integration for column `c`:

```text
Phi_half[c,K] = Phi_s[c]
Phi_full, Phi_half = existing upward hydrostatic integral
```

No pressure, temperature, geopotential, or terrain value is added to the prognostic
state. Interface pressure remains

```text
p_half[c,k] = A[k] + B[k]*ps[c]
```

and the Phase 4/5 horizontal-convergence recurrence remains the only source of
`ps_dot` and coordinate-relative interface mass flux. Top and bottom relative mass
fluxes remain zero. The geometric atmosphere follows the fixed surface because a
parcel at the lower boundary may have nonzero geometric vertical velocity while its
velocity relative to the fixed terrain-following coordinate is zero.

The dry momentum source continues to use the ADR 0006 identity and the same
least-squares spherical gradient for `Phi_full` and `p_full`. Phase 6 does not add an
exact-rest perturbation formulation or a second well-balanced solver. DCMIP 2-0-0
measures the resulting cancellation error and its horizontal and vertical refinement.
If that error does not converge at the pre-registered rate, the pressure-gradient
form is reconsidered in a separate ADR rather than patched with clipping, a hidden
filter, or a benchmark-specific cancellation.

For `Phi_s[c] = 0`, diagnosis, RHS values, time steps, state updates, and outputs must
be exactly the Phase 5 path. A spatially uniform nonzero `Phi_s` must shift diagnosed
geopotential and total potential energy by the analytic constant while leaving every
dynamics tendency unchanged to rounding.

### Shallow-water topographic source

Williamson test 5 uses the existing Rusanov shallow-water path and adds one separately
reported momentum tendency

```text
S_orography = -h * grad(Phi_s)
```

where `h` is fluid depth and `Phi_s = g*h_s`. The same cell-centred least-squares
gradient used elsewhere supplies `grad(Phi_s)`. It changes neither mass flux nor the
checkpoint state. The total momentum tendency is the existing flux, pressure,
Coriolis, and diffusion sum plus `S_orography`.

Phase 6 does not retrofit the unselected compatible candidate or introduce
hydrostatic reconstruction solely for this regression. A non-flat shallow-water
configuration with the compatible scheme is rejected explicitly. Flat configurations
of both schemes remain unchanged.

### Minimal imported-data path and filter

Phase 6 accepts one external format: a strict rectilinear latitude--longitude CSV with
the exact columns

```text
longitude_deg,latitude_deg,height_m
```

The tensor-product grid must cover one periodic longitude cycle and the requested
latitude range, contain every pair exactly once, and contain only finite values.
Bilinear interpolation is periodic in longitude. No NetCDF, GeoTIFF, missing-value
repair, projection library, or online download is introduced.

Interpolation samples cell centres. Optional smoothing then applies a fixed,
synchronous, unique-edge diffusion pass. Each pass exchanges equal and opposite
area-integrated increments across an edge. Its globally chosen stable coefficient
makes the update convex, so the pass preserves constants, the area-weighted mean, and
the input extrema while treating panel seams and corners like ordinary neighbours.
Only the integer pass count is configurable. Published analytic benchmarks bypass
this import/filter path.

### Diagnostics, not duplicate forces

Terrain diagnostics add:

- minimum and maximum surface height and maximum resolved slope;
- pressure-gradient cancellation norms and maximum spurious wind;
- dry-core pressure-gradient axial torque;
- boundary mountain-torque estimate and their residual;
- terrain pressure-work conversion.

For physical horizontal position vector `r_c = R*cell_center[c]`, cell area `A_c`,
layer pressure source `S_pg`,
and surface gradient `grad(Phi_s)`, diagnose

```text
tau_model = sum_c A_c * r_c cross sum_k(S_pg[c,k])
tau_boundary = sum_c A_c * r_c cross (-(ps[c]/g)*grad(Phi_s[c]))
W_terrain = sum_c A_c * sum_k u[c,k] dot (-M[c,k]*grad(Phi_s[c]))
```

The axial components and `tau_model - tau_boundary` are reported. These are
diagnostics of the pressure force already present in the momentum equation; they are
never applied as additional tendencies. Fixed, impermeable terrain performs no
external moving-boundary work. `W_terrain` is labelled as conversion between resolved
kinetic and potential/internal energy, not as an external energy source.

For the bounded hydrostatic mountain-wave case, absolute interface pressure velocity
may be diagnosed from existing quantities as

```text
omega_abs = B*ps_dot + u_half dot grad_eta(p_half) + g*F_relative
```

and the pressure-coordinate momentum flux is reported. Phase 6 does not add vertical
momentum or claim a nonhydrostatic `w` solution.

### Publication boundary

Phase 6 numerical acceptance is C++-only. `FrameV2` has a fixed payload and remains
unchanged; the Web UI, gateway, controls, and checkpoint layouts are not extended to
display terrain. A future measured visualization requirement may define a new frame
schema separately.

## Consequences

- Terrain changes one immutable driver input and derived geopotential, not the
  atmospheric state or continuity coupling.
- The flat case remains an exact regression of Phase 5 rather than a parallel code
  path.
- DCMIP 2-0-0 exposes the actual pressure-gradient truncation error; Phase 6 does not
  make that benchmark pass by construction.
- Williamson test 5 validates the selected shallow-water reference method without
  reviving the compatible candidate.
- External terrain support is deliberately narrow but reproducible and seam-safe.
- Moving terrain, surface stress/heat exchange, nonhydrostatic mountain waves, sponge
  layers, multiple raster formats, and terrain visualization remain outside Phase 6.
