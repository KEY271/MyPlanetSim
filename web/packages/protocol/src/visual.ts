import { ProtocolError, protocolVersion } from "./index";

export type VisualSourceKind = "offline_csv" | "live_run";
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
  readonly schemaVersion: typeof protocolVersion;
  readonly sourceKind: VisualSourceKind;
  readonly grid: {
    readonly topology: "cubed_sphere";
    readonly mapping: "equiangular_gnomonic_v1";
    readonly cellsPerPanel: number;
    readonly flattenOrder: "panel_major_then_j_then_i";
  };
  readonly frame: {
    readonly timeSeconds: number | null;
    readonly step: number | null;
    readonly configFingerprint: string | null;
  };
  readonly fields: readonly VisualFieldV1[];
}

export interface VisualFrameV2 {
  readonly schemaVersion: 2; readonly cellsPerPanel: number; readonly levels: number;
  readonly timeSeconds: number; readonly step: number; readonly configFingerprint: string;
  readonly surfacePressurePa: Float64Array; readonly pressurePa: Float64Array;
  readonly potentialTemperatureK: Float64Array; readonly temperatureK: Float64Array;
  readonly windXMetersPerSecond: Float64Array; readonly windYMetersPerSecond: Float64Array;
  readonly windZMetersPerSecond: Float64Array; readonly tracerMixingRatio: Float64Array;
}

export type FrameV2Field = "pressure" | "potential_temperature" | "temperature" | "tracer" | "wind_speed";
export type FrameV2FieldId = FrameV2Field | "surface_pressure";

export interface FrameV2FieldDescriptor {
  readonly id: FrameV2FieldId;
  readonly label: string;
  readonly unit: string;
  readonly kind: VisualFieldKind;
  // A surface field has no model level, so the viewer disables the level selector and
  // the column profile instead of inventing a vertical copy of a single value.
  readonly volume: boolean;
}

export const frameV2Fields: readonly FrameV2FieldDescriptor[] = Object.freeze([
  { id: "surface_pressure", label: "Surface pressure", unit: "Pa", kind: "scalar", volume: false },
  { id: "pressure", label: "Pressure", unit: "Pa", kind: "scalar", volume: true },
  { id: "potential_temperature", label: "Potential temperature", unit: "K", kind: "scalar", volume: true },
  { id: "temperature", label: "Temperature", unit: "K", kind: "scalar", volume: true },
  { id: "tracer", label: "Tracer", unit: "1", kind: "scalar", volume: true },
  { id: "wind_speed", label: "Wind speed", unit: "m/s", kind: "derived", volume: true },
].map((field) => Object.freeze(field)) as readonly FrameV2FieldDescriptor[]);

export function frameV2Field(id: string): FrameV2FieldDescriptor {
  const field = frameV2Fields.find((candidate) => candidate.id === id);
  if (!field) fail(`unknown FrameV2 field ${id}`);
  return field;
}

export function decodeFrameV2(bytes: ArrayBuffer, expectedFingerprint?: string): VisualFrameV2 {
  const view = new DataView(bytes);
  if (view.byteLength < 60 || String.fromCharCode(...new Uint8Array(bytes.slice(0, 8))) !== "MPSFRAM2" || view.getUint32(8, true) !== 2 || view.getUint32(12, true) !== 0) fail("invalid FrameV2 header");
  const n=Number(view.getBigUint64(16,true)),levels=Number(view.getBigUint64(24,true)),timeSeconds=finite(view.getFloat64(32,true)),step=Number(view.getBigUint64(40,true)),cells=Number(view.getBigUint64(48,true)),fingerprintLength=view.getUint32(56,true),payload=60+fingerprintLength,volume=cells*levels;
  if(!Number.isInteger(n)||n<1||n>24||!Number.isInteger(levels)||levels<1||levels>30||cells!==6*n*n||fingerprintLength<1||payload+8*(cells+7*volume)!==view.byteLength)fail("invalid FrameV2 shape");
  const configFingerprint=new TextDecoder().decode(new Uint8Array(bytes,60,fingerprintLength));if(expectedFingerprint!==undefined&&configFingerprint!==expectedFingerprint)fail("frame fingerprint mismatch");let cursor=payload;
  const read=(length:number)=>{const values=new Float64Array(length);for(let i=0;i<length;i+=1){values[i]=finite(view.getFloat64(cursor,true));cursor+=8;}return values;};
  return Object.freeze({schemaVersion:2 as const,cellsPerPanel:n,levels,timeSeconds,step,configFingerprint,surfacePressurePa:read(cells),pressurePa:read(volume),potentialTemperatureK:read(volume),temperatureK:read(volume),windXMetersPerSecond:read(volume),windYMetersPerSecond:read(volume),windZMetersPerSecond:read(volume),tracerMixingRatio:read(volume)});
}
export function decodeFrame(bytes:ArrayBuffer,expectedFingerprint?:string):VisualDatasetV1|VisualFrameV2{const magic=String.fromCharCode(...new Uint8Array(bytes.slice(0,8)));if(magic==="MPSFRAM1")return decodeFrameV1(bytes,expectedFingerprint);if(magic==="MPSFRAM2")return decodeFrameV2(bytes,expectedFingerprint);fail("unsupported frame magic");}
export function frameV2CellCount(frame: VisualFrameV2): number {
  return 6 * frame.cellsPerPanel * frame.cellsPerPanel;
}

function volumeValue(frame: VisualFrameV2, field: FrameV2Field, offset: number): number {
  switch (field) {
    case "pressure": return frame.pressurePa[offset];
    case "potential_temperature": return frame.potentialTemperatureK[offset];
    case "temperature": return frame.temperatureK[offset];
    case "tracer": return frame.tracerMixingRatio[offset];
    case "wind_speed": return Math.hypot(frame.windXMetersPerSecond[offset],
      frame.windYMetersPerSecond[offset], frame.windZMetersPerSecond[offset]);
  }
}

export function levelSlice(frame: VisualFrameV2, field: FrameV2Field, level: number): Float64Array {
  if (!Number.isInteger(level) || level < 0 || level >= frame.levels) fail("level is out of range");
  const cells = frameV2CellCount(frame);
  const values = new Float64Array(cells);
  for (let cell = 0; cell < cells; cell += 1) {
    values[cell] = volumeValue(frame, field, cell * frame.levels + level);
  }
  return values;
}

export interface FrameV2Column {
  readonly pressurePa: Float64Array;
  readonly values: Float64Array;
}

export function selectedColumn(frame: VisualFrameV2, field: FrameV2Field, cell: number): FrameV2Column {
  if (!Number.isInteger(cell) || cell < 0 || cell >= frameV2CellCount(frame)) fail("cell is out of range");
  const values = new Float64Array(frame.levels);
  const pressurePa = new Float64Array(frame.levels);
  for (let level = 0; level < frame.levels; level += 1) {
    const offset = cell * frame.levels + level;
    values[level] = volumeValue(frame, field, offset);
    pressurePa[level] = frame.pressurePa[offset];
  }
  return Object.freeze({ pressurePa, values });
}

// The renderers only understand VisualDatasetV1, so a V2 level is adapted into a
// single-field dataset. The grid identity is unchanged, which keeps the memoized geometry
// caches in geometry.ts alive across field, level, and frame changes.
export function frameV2LevelDataset(
  frame: VisualFrameV2,
  fieldId: FrameV2FieldId,
  level: number,
): VisualDatasetV1 {
  const field = frameV2Field(fieldId);
  const values = field.volume ? levelSlice(frame, fieldId as FrameV2Field, level) : frame.surfacePressurePa;
  return dataset("live_run", frame.cellsPerPanel, frame.timeSeconds, frame.step,
    frame.configFingerprint, [{ id: field.id, label: field.label, unit: field.unit,
      kind: field.kind, provenance: "FrameV2", values }]);
}

export interface FrameV2Sample {
  readonly value: number;
  readonly pressurePa: number | null;
}

// The inspector reports the decoded C++ number for one cell without rescaling it, so a
// reader can compare the view against the solver output directly.
export function frameV2Sample(
  frame: VisualFrameV2,
  fieldId: FrameV2FieldId,
  cell: number,
  level: number,
): FrameV2Sample {
  const field = frameV2Field(fieldId);
  if (!Number.isInteger(cell) || cell < 0 || cell >= frameV2CellCount(frame)) fail("cell is out of range");
  if (!field.volume) return Object.freeze({ value: frame.surfacePressurePa[cell], pressurePa: null });
  if (!Number.isInteger(level) || level < 0 || level >= frame.levels) fail("level is out of range");
  const offset = cell * frame.levels + level;
  return Object.freeze({ value: volumeValue(frame, fieldId as FrameV2Field, offset),
    pressurePa: frame.pressurePa[offset] });
}

const panels = ["PX", "PY", "NX", "NY", "PZ", "NZ"] as const;
const panelIndex = new Map<string, number>(panels.map((panel, index) => [panel, index]));

function fail(message: string): never {
  throw new ProtocolError("invalid_dataset", message);
}

function finite(value: number): number {
  if (!Number.isFinite(value)) fail("dataset contains a non-finite value");
  return value;
}

function dataset(
  sourceKind: VisualSourceKind,
  cellsPerPanel: number,
  timeSeconds: number | null,
  step: number | null,
  configFingerprint: string | null,
  fields: readonly VisualFieldV1[],
): VisualDatasetV1 {
  if (!Number.isInteger(cellsPerPanel) || cellsPerPanel < 1 || cellsPerPanel > 1024) {
    fail("cellsPerPanel is outside the supported range");
  }
  const cellCount = 6 * cellsPerPanel * cellsPerPanel;
  for (const field of fields) {
    if (field.values.length !== cellCount) fail(`field ${field.id} has the wrong length`);
    for (const value of field.values) finite(value);
  }
  return Object.freeze({
    schemaVersion: protocolVersion,
    sourceKind,
    grid: Object.freeze({ topology: "cubed_sphere", mapping: "equiangular_gnomonic_v1",
      cellsPerPanel, flattenOrder: "panel_major_then_j_then_i" }),
    frame: Object.freeze({ timeSeconds, step, configFingerprint }),
    fields: Object.freeze(fields.map((field) => Object.freeze({ ...field,
      values: new Float64Array(field.values) }))),
  });
}

export function decodeFrameV1(bytes: ArrayBuffer, expectedFingerprint?: string): VisualDatasetV1 {
  const view = new DataView(bytes);
  const magic = "MPSFRAM1";
  if (view.byteLength < 52 || String.fromCharCode(...new Uint8Array(bytes.slice(0, 8))) !== magic) {
    fail("invalid FrameV1 magic");
  }
  if (view.getUint32(8, true) !== 1 || view.getUint32(12, true) !== 0) fail("invalid FrameV1 header");
  const cellsPerPanel = Number(view.getBigUint64(16, true));
  const timeSeconds = view.getFloat64(24, true);
  const step = Number(view.getBigUint64(32, true));
  const cellCount = Number(view.getBigUint64(40, true));
  const fingerprintLength = view.getUint32(48, true);
  const payloadOffset = 52 + fingerprintLength;
  if (!Number.isSafeInteger(cellsPerPanel) || cellCount !== 6 * cellsPerPanel * cellsPerPanel ||
      fingerprintLength < 1 || payloadOffset + 32 * cellCount !== view.byteLength) {
    fail("FrameV1 shape or byte length is invalid");
  }
  const fingerprint = new TextDecoder().decode(new Uint8Array(bytes, 52, fingerprintLength));
  if (expectedFingerprint !== undefined && fingerprint !== expectedFingerprint) fail("frame fingerprint mismatch");
  const values = ["depth", "momentum_x", "momentum_y", "momentum_z"].map((id, fieldIndex) => {
    const fieldValues = new Float64Array(cellCount);
    for (let cell = 0; cell < cellCount; cell += 1) {
      fieldValues[cell] = finite(view.getFloat64(payloadOffset + (fieldIndex * cellCount + cell) * 8, true));
    }
    return { id, label: id, unit: id === "depth" ? "m" : "m²/s", kind: "scalar" as const,
      provenance: "FrameV1", values: fieldValues };
  });
  return dataset("live_run", cellsPerPanel, finite(timeSeconds), step, fingerprint, values);
}

function csvRows(text: string, expectedHeader: readonly string[]): string[][] {
  const lines = text.split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  if (lines.length === 0 || lines[0].split(",").map((value) => value.trim()).join(",") !== expectedHeader.join(",")) {
    fail("CSV header is invalid");
  }
  return lines.slice(1).map((line) => line.split(",").map((value) => value.trim()));
}

function cellIndex(panel: string, i: number, j: number, cellsPerPanel: number): number {
  const panelNumber = panelIndex.get(panel);
  if (panelNumber === undefined || !Number.isInteger(i) || !Number.isInteger(j) ||
      i < 0 || i >= cellsPerPanel || j < 0 || j >= cellsPerPanel) fail("CSV cell identifier is invalid");
  return panelNumber * cellsPerPanel * cellsPerPanel + j * cellsPerPanel + i;
}

function csvDataset(
  text: string,
  cellsPerPanel: number,
  sourceKind: VisualSourceKind,
  header: readonly string[],
  fields: readonly { id: string; label: string; unit: string; column: number }[],
): VisualDatasetV1 {
  const rows = csvRows(text, header);
  const cellCount = 6 * cellsPerPanel * cellsPerPanel;
  if (rows.length !== cellCount) fail("CSV does not contain exactly one row per cell");
  const values = fields.map((field) => ({ ...field, kind: "scalar" as const,
    provenance: "offline CSV", values: new Float64Array(cellCount) }));
  const seen = new Set<number>();
  for (const row of rows) {
    if (row.length !== header.length) fail("CSV row has the wrong number of columns");
    const index = cellIndex(row[0], Number(row[1]), Number(row[2]), cellsPerPanel);
    if (seen.has(index)) fail("CSV contains a duplicate cell");
    seen.add(index);
    for (const field of values) field.values[index] = finite(Number(row[field.column]));
  }
  if (seen.size !== cellCount) fail("CSV is missing a cell");
  return dataset(sourceKind, cellsPerPanel, null, null, null, values);
}

export function parseShallowWaterCsv(text: string, cellsPerPanel: number): VisualDatasetV1 {
  return csvDataset(text, cellsPerPanel, "offline_csv",
    ["panel", "i", "j", "depth_m", "momentum_x", "momentum_y", "momentum_z"],
    [{ id: "depth", label: "Depth", unit: "m", column: 3 },
      { id: "momentum_x", label: "Momentum X", unit: "m²/s", column: 4 },
      { id: "momentum_y", label: "Momentum Y", unit: "m²/s", column: 5 },
      { id: "momentum_z", label: "Momentum Z", unit: "m²/s", column: 6 }]);
}

export function parseTransportCsv(text: string, cellsPerPanel: number): VisualDatasetV1 {
  return csvDataset(text, cellsPerPanel, "offline_csv", ["panel", "i", "j", "tracer"],
    [{ id: "tracer", label: "Tracer", unit: "1", column: 3 }]);
}

export function addShallowWaterDerivedFields(
  value: VisualDatasetV1,
  initial?: VisualDatasetV1,
): VisualDatasetV1 {
  const depth = value.fields.find((field) => field.id === "depth");
  const x = value.fields.find((field) => field.id === "momentum_x");
  const y = value.fields.find((field) => field.id === "momentum_y");
  const z = value.fields.find((field) => field.id === "momentum_z");
  if (!depth || !x || !y || !z) fail("shallow-water fields are incomplete");
  const speed = new Float64Array(depth.values.length);
  const anomaly = new Float64Array(depth.values.length);
  const initialDepth = initial?.fields.find((field) => field.id === "depth")?.values;
  if (initial && !initialDepth) fail("initial dataset has no depth field");
  for (let index = 0; index < speed.length; index += 1) {
    if (!(depth.values[index] > 0)) fail("depth must be positive for speed");
    speed[index] = Math.hypot(x.values[index], y.values[index], z.values[index]) / depth.values[index];
    anomaly[index] = initialDepth ? depth.values[index] - initialDepth[index] : depth.values[index];
  }
  return dataset(value.sourceKind, value.grid.cellsPerPanel, value.frame.timeSeconds,
    value.frame.step, value.frame.configFingerprint, [...value.fields,
      { id: "speed_m_s", label: "Speed", unit: "m/s", kind: "derived", provenance: "derived from momentum/depth", values: speed },
      { id: "depth_anomaly_m", label: "Depth anomaly", unit: "m", kind: "derived", provenance: "depth minus initial depth", values: anomaly }]);
}

export function fieldRange(field: VisualFieldV1): readonly [number, number] {
  if (field.values.length === 0) fail("field is empty");
  let minimum = field.values[0];
  let maximum = field.values[0];
  for (const value of field.values) { minimum = Math.min(minimum, value); maximum = Math.max(maximum, value); }
  return [minimum, maximum];
}

export function inspectFieldCell(datasetValue: VisualDatasetV1, fieldId: string, cell: number): number {
  const field = datasetValue.fields.find((candidate) => candidate.id === fieldId);
  if (!field || !Number.isInteger(cell) || cell < 0 || cell >= field.values.length) fail("field cell is invalid");
  return field.values[cell];
}
