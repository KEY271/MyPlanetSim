import { ProtocolError } from "./errors";

export type VisualFieldKind = "scalar" | "derived";

export interface VisualFieldV1 {
  readonly id: string;
  readonly label: string;
  readonly unit: string;
  readonly kind: VisualFieldKind;
  readonly provenance: string;
  readonly values: Float64Array;
}

export interface VisualDatasetV1 {
  readonly schemaVersion: 1;
  readonly sourceKind: "visual_dataset";
  readonly grid: {
    readonly topology: "cubed_sphere";
    readonly mapping: "equiangular_gnomonic_v1";
    readonly cellsPerPanel: number;
    readonly flattenOrder: "panel_major_then_j_then_i";
  };
  readonly frame: {
    readonly timeSeconds: number | null;
    readonly step: number | null;
    readonly configFingerprint: string;
  };
  readonly fields: readonly VisualFieldV1[];
}

export function fieldRange(field: VisualFieldV1): readonly [number, number] {
  if (field.values.length === 0) {
    throw new ProtocolError("invalid_visual_dataset", `field ${field.id} is empty`);
  }
  let minimum = Infinity;
  let maximum = -Infinity;
  for (const value of field.values) {
    if (!Number.isFinite(value)) {
      throw new ProtocolError("invalid_visual_dataset", `field ${field.id} contains a non-finite value`);
    }
    minimum = Math.min(minimum, value);
    maximum = Math.max(maximum, value);
  }
  return [minimum, maximum];
}
