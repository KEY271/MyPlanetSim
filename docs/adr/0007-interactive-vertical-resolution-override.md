# ADR 0007: interactive vertical resolution override

**Status:** superseded by ADR 0028 in Phase 15

P15.06 removed interactive runs and their `N`/`K` overrides. Dataset dimensions now come
from the immutable saved manifest, so this control-path decision is retained only as a
record of the former Phase 5 interface.

## Context

[ADR 0005](0005-hybrid-vertical-coordinate-and-column-state.md) fixed the hybrid
sigma-pressure coordinate and made the `A/B` arrays explicit, immutable preset data. The
Phase 5 control protocol already lets an interactive run override the horizontal
resolution `N`, because a cubed-sphere grid is fully determined by one integer.

The vertical resolution `K` is not. Configuration validation requires
`A_pa.size() == B.size() == K + 1`, so a request cannot change `K` without also supplying
a coordinate. A reader who wants to see how the dry hydrostatic column behaves at a finer
vertical resolution therefore has to edit a configuration file and restart the gateway,
which the horizontal resolution does not require.

Two ways to close that asymmetry were considered.

1. Ship one configuration file per level count and select between them. This keeps ADR
   0005 untouched, but the available `K` values become a fixed, discrete list and each new
   value is a new file to maintain.
2. Let the control request carry `K` and derive the coordinate from it. This makes `K`
   continuous over the interactive range, at the cost of a coordinate that is computed
   rather than written down.

The second was chosen. The risk it introduces is that a preset with a deliberately
designed coordinate — a stretched boundary layer, a sponge near the model top — could be
silently replaced by a generated one, and the run would still look plausible.

## Decision

### The override

A control request may carry an optional `control.levels`. Absent or `0` means the
configured `A/B` are used unchanged, so every control request written before this decision
keeps its exact meaning. A nonzero value is bounded by `1 <= K <= 30`, matching the
interactive FrameV2 allowlist in the Phase 5 plan.

When `control.levels` is nonzero, the coordinate is regenerated as a uniform sigma ramp
that preserves the model top pressure the preset already declared:

```text
p_top   = A_pa[0]   from the configured preset
A_pa[k] = p_top * (1 - k/K)
B[k]    = k/K,                 k = 0 .. K
```

This satisfies the ADR 0005 endpoint requirements by construction: `A_pa[0] = p_top`,
`B[0] = 0`, `A_pa[K] = 0`, `B[K] = 1`. The generated arrays are then validated by the same
`HybridPressureCoefficients::validate` as any configured coordinate, so a `K` whose layers
would violate the minimum pressure thickness over the configured surface-pressure interval
is rejected rather than run.

### The guard

The override is accepted only when the preset's own `A/B` already equal the uniform sigma
ramp at its configured `K`, within a scale-aware tolerance. Otherwise the request is
rejected with an error naming the preset coordinate.

This is the part that keeps the decision honest. Refining a uniform preset stays inside
the coordinate family the preset already chose, and it reproduces the preset exactly when
the requested `K` equals the configured one. A preset with a designed non-uniform
coordinate cannot be flattened by a browser control; it needs either its own explicit
configuration at the desired resolution or a later ADR that defines how to refine it.

### Provenance

The resolved `A/B` and `K` are written into the canonical configuration text before the
run starts, so they enter the configuration fingerprint, the run metadata, the checkpoint
header, and the published FrameV2 metadata. A run at `K = 16` is a distinguishable
configuration, not the preset's `K = 8` configuration with a hidden substitution.

## Consequences

- `K` becomes an interactive control alongside `N` for uniform sigma presets.
- Nothing about a standalone configuration-file run changes; the override exists only in
  the machine control path.
- Vertical convergence can be examined from the viewer, but a comparison across `K` is a
  comparison across configurations and must be reported with the fingerprint of each.
- Non-uniform coordinates remain configuration-only. Supporting a refinement rule for them
  requires a new ADR that states the rule and its validation, not a change here.
