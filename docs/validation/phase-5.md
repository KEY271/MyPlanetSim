# Phase 5 validation

## C++ core gate

Validated with:

```sh
cmake --preset dev
cmake --build --preset dev -j4
ctest --preset dev --output-on-failure
cmake --build build/dev --target format-check
build/dev/my_planet_sim --config configs/phase5_isothermal_rest.cfg
```

The unit gate covers cell-major flattening, strict configuration, derived hybrid-pressure
columns, constant-scalar Rusanov consistency, uniform pressure-gradient cancellation,
the Phase 4 continuity recurrence, tangent projection, checkpoint layout, and deterministic
budget reductions. The coupled rest preset completes 60 steps to 600 s with unchanged
state and finite dry-mass and energy diagnostics.

The DCMIP and UMJS presets in this phase are bounded initialization/regression cases. They
do not claim full published-test conformance or a reference envelope; quantitative long-run
convergence at production resolution remains a known validation gap. The reference core is
explicit and intentionally has no implicit solver, terrain, physics, multi-tracer registry,
or vertical mixing.

## ADR 0006 implementation audit

The Phase 7 entry gate requires the dry driver to match the
[ADR 0006](../adr/0006-dry-hydrostatic-state-and-coupling.md) contract before any forcing
is added. Three mismatches were found and closed on 2026-09-04.

| ADR 0006 clause | Previous implementation | Now |
|---|---|---|
| Every SSP-RK3 stage performs the complete sequence | one forward Euler update per step | three stages, each re-evaluating the complete RHS from its own stage state, with tangent projection and invariant validation at every stage boundary |
| One step size satisfies the horizontal CFL, the Phase 4 vertical transport CFL, surface-pressure bounds, and the requested maximum | horizontal CFL and requested maximum only | `DryHydrostaticRhs` also reports the per-column vertical stable step; the step is limited by the surface-pressure headroom and halved and retried when any stage leaves the configured pressure range or exceeds the stage CFL |
| Each model level uses linear reconstruction and the Barth--Jespersen limiter | piecewise-constant face values | `dry_hydrostatic_reconstruction.cpp` reconstructs `M`, tangent `u`, `theta`, `q`, and `T` per level and honours `dry_hydrostatic.reconstruction` and `dry_hydrostatic.limiter` |

`test_dry_hydrostatic_driver.cpp` fixes the first two: one test reproduces all three stages
independently and requires `advance` to match it, and one requires a CFL-limited run to take
more than one step, to land exactly on the requested end time, and to stay inside the
pressure bounds.

Turning the reconstruction on had two consequences worth recording.

The Barth--Jespersen limiter now also applies to the tangent velocity, with the same face
normal/tangent bounds and speed factor the shallow-water path uses. The first version left
the least-squares vector gradient unlimited, which is not what ADR 0006 selects.

Second-order faces also let a mixing ratio that is physically zero pick up a value of
order 1e-75 from flux cancellation, which the exact `q >= 0` invariant rejected. The
invariant is now checked to a 1e-12 roundoff tolerance, as tangency already was. This
surfaced as `configs/phase5_solid_body_transport.cfg` failing while every registered test
passed, because the shipped presets are not otherwise exercised by the suite;
`smoke.phase5_solid_body_transport` now runs that preset as part of `phase5_gate`.

## Unforced flat baseline

`test_dry_hydrostatic_baseline.cpp` (`phase5.dry_hydrostatic_baseline`, label
`phase5_gate`) registers the reference values Phase 7 forced runs are compared against. The
case is the flat UMJS14 steady preset at `N = 4`, `K = 8`, integrated to 600 s with a 10 s
requested step.

| Quantity | Registered bound on relative drift | Observed |
|---|---:|---:|
| Dry mass | 1e-13 | 5.0e-15 |
| Potential-temperature mass | 1e-13 | 1.0e-15 |
| Tracer mass | 1e-13 | 0 |
| Total energy | 2e-6 | 1.1e-6 |
| Absolute axial angular momentum | 5e-5 | 3.1e-5 |

Energy and angular momentum are bounded rather than conserved, as ADR 0006 states for the
Rusanov/least-squares reference scheme. Restart is exact: a run interrupted at 300 s,
written and read back through the `dry_hydrostatic_cell_column_v1` checkpoint payload, and
continued to 600 s reproduces the uninterrupted run bit for bit.

These are CI-scale registrations. The quantitative 3D convergence study and the UMJS14
published reference envelope remain the open Phase 5 gap above.

## Visualizer gate

Passing this section alone is not a Phase 5 numerical validation claim; it gates only the
browser-facing path from a validated C++ state to the screen.

Validated with:

```sh
cd web
npm ci
npm run typecheck && npm run lint && npm run test && npm run build
npm run test:live --workspace @myplanetsim/gateway
```

### Vertical resolution override

[ADR 0007](../adr/0007-interactive-vertical-resolution-override.md) lets a control request
carry `control.levels`. Absent or `0` keeps the configured coordinate, so every request
written before it is unchanged. A nonzero value regenerates the uniform sigma ramp from the
preset's own model top, and is accepted only when the preset already uses that ramp.
`test_vertical_column.cpp` gates the generated endpoints, the exact interior values, and
rejection of a stretched or slightly nudged coordinate. `test_dry_hydrostatic_io.cpp` gates
that no override leaves the coordinate byte-identical, that requesting the preset's own `K`
is the identity down to the configuration fingerprint, that a refined `K` changes the
fingerprint, and that a stretched preset is refused rather than flattened.

### Native publication

`experiment.kind = dry_hydrostatic` accepts `--control-request` and `--event-stream ndjson`.
The driver's existing observer and cancel hooks publish a FrameV2 file at the requested
frame interval and again for the final step, each announced by a `frame.ready` event that
carries `frameSchemaVersion: 2`. Initial-condition edits are rejected because the dry core
has no depth field. `--describe-control` reports the supported model kinds, frame schema
versions, and the interactive dry limits (`N <= 24`, `K <= 30`).

### Decoding and slicing

`web/packages/protocol/src/frame-v2.test.ts` builds the fixed byte layout directly and
checks magic dispatch against FrameV1, per-field offsets in `cell-major, then k` order,
agreement between a level slice, a selected column, and a single sample, wind speed as a
derived magnitude, and the adapted single-field dataset for every field. Rejection is
gated for a wrong schema version, `N > 24`, `K > 30`, a mismatched declared cell count,
truncation, trailing bytes, a non-finite value, a fingerprint mismatch, and an
out-of-range level or cell.

### Viewer

`web/apps/ui/src/visualization.test.ts` gates that the map, globe, profile, and inspector
report the same decoded number for every cell and level; that time, step, and fingerprint
follow each adapted level; that the top level, the bottom level, and panel-corner and
seam cells are ordinary members of the slice; that a surface field, a constant field, and
a signed field display correctly; that the memoized cell and edge geometry is identical
across field, level, and frame changes; and that pressure is plotted downward on a
logarithmic axis with a constant column centred rather than divided by a zero span. The
narrow-viewport and reduced-motion rules are gated against the stylesheet.
`controller.test.ts` gates that a FrameV2 run is stored as published, with no
shallow-water derivation applied, alongside the unchanged FrameV1 mock regression.

### Gateway

`web/apps/gateway/test/gateway.test.js` gates the preset descriptors returned by
`/api/v1/capabilities` (without leaking a configuration path), the shallow-water fallback
for an undescribed preset, the per-preset N limit and edit rejection, the exact FrameV2
byte size `8 * C * (1 + 7K)`, the pre-flight and streaming byte budget, and the startup
compatibility check between `--describe-control` and the configured preset descriptors.
The existing loopback auth, shell-free argv, path containment, reconnect, cancel, and
bundle regressions are unchanged.

### Live end to end

`web/apps/gateway/test/live.test.js` runs the native binary against
`configs/phase5_visualizer_rest_n4.cfg` at `N = 4`, `K = 8`, and again at `K = 20` through
the ADR 0007 override. It checks the advertised
preset descriptor, rejection of a `gaussian_depth` edit, `run.completed` with every frame
announced as schema 2, the exact frame byte length, frame 0 matching the initialized state
before integration (uniform 100000 Pa surface pressure, pressure increasing downward,
288 K, zero wind), a later frame advancing time and step, and a cancelled run whose
announced frames are all still servable at their announced length. The refined run is
checked for `K = 20` in the frame header, the matching byte length, a downward-increasing
generated coordinate, a configuration fingerprint distinct from the preset's, and rejection
of `K = 0` and `K = 31` before the solver is started. The shallow-water
FrameV1 live run and cancel in the same file are unchanged.

### Browser gate

Every other gate drives the viewer through node, which accepts API shapes a browser
rejects and never renders. Two defects reached the viewer that way.

The HTTP client called `fetch` as a method of the client. A browser's `fetch` is a WebIDL
operation on Window and rejects a foreign receiver with `Illegal invocation`, so every
request from the browser failed. `frame-v2.test.ts` now asserts that the client never
invokes `fetch` with itself as the receiver, and drives the real client over a real socket
through a `fetch` that enforces the browser rule; both fail with the browser error if the
binding is removed.

Separately, opening the plain dev-server URL instead of the tokenised one printed by
`just dev` silently selected the offline demo, which serves only the shallow-water preset.
That is indistinguishable from the dry hydrostatic model being absent, and `--open` cannot
repair it because reusing an already-open tab is a same-document navigation, so the session
is never re-read. `just dev` now hands the session to the dev server, which injects it into
the client, so the URL Vite prints is the live one and no fragment is involved. A dev server
started without a session still falls back to the demo and now says so in a banner.

The third defect was that the only dry preset on offer was `isothermal_rest`, which is
exactly steady by construction: every published frame was identical, so the viewer rendered
a still image. `configs/interactive_dry_presets.txt` is now the single list of presets the
viewer offers, read by both the Justfile and the live gate.

```sh
npx playwright install chromium
npm run test:browser --workspace @myplanetsim/ui
```

`web/apps/ui/test/browser.test.js` starts the native gateway and two dev servers, one with
an injected session and one without, then opens the real pages in Chromium. It checks that
the session-less server announces the offline demo and offers only the shallow-water preset;
that the plain URL of the session-carrying server reports the native engine and offers only
the steady dry preset in addition to it; that selecting the dry preset announces the level
count and shows the level and `K` controls before any run; and that a `K = 12` run yields a
level slider spanning `0..11`, a column profile, an inspected value that remains steady
between frames, and an inspector that follows the level slider. It skips when the binary,
the preset, or the Chromium download is absent.

`live.test.js` additionally requires `configs/interactive_dry_presets.txt` to contain only
`phase5_visualizer_rest_n4` and verifies that its surface pressure does not change. The
steady case proves the coupled core leaves a resting atmosphere alone.

### Known gaps

- The offline/mock client publishes FrameV1 only, so the dry hydrostatic viewer requires
  the native gateway; there is no browser-side dry model to regress against.
- The browser gate runs Chromium only. Firefox and WebKit coverage, and a WebGL-failure
  path for the profile panel, remain the Phase 3 automation gap.
- The interactive allowlist is bounded at `N <= 24`, `K <= 30` and a 256 MiB run budget.
  Partial frame retrieval is deliberately not designed until a measurement needs it.
