import { ProtocolError } from "./errors";

export type DatasetFieldLocation = "surface" | "atmosphere";

export interface DatasetField {
  readonly id: string;
  readonly unit: string;
  readonly location: DatasetFieldLocation;
  readonly valueOffset: number;
  readonly valueCount: number;
}

export interface DatasetPeriod {
  readonly index: number;
  readonly scheduledStartS: number;
  readonly scheduledEndS: number;
  readonly actualStartS: number;
  readonly actualEndS: number;
  readonly weightS: number;
  readonly coverage: number;
  readonly complete: boolean;
  readonly path: string;
  readonly byteLength: number;
}

export interface VisualDatasetManifest {
  readonly schemaVersion: 1;
  readonly modelKind: "dry_hydrostatic";
  readonly cellsPerPanel: number;
  readonly cellCount: number;
  readonly levelCount: number;
  readonly arrayOrder: "field-cell-level";
  readonly planet: { readonly radiusM: number; readonly gravityMS2: number };
  readonly hybridAHalfPa: readonly number[];
  readonly hybridBHalf: readonly number[];
  readonly configFingerprint: string;
  readonly terrainSourceFingerprint: string;
  readonly averagingRule: "accepted_step_end_time_weighted_rectangle";
  readonly terrain: {
    readonly path: string;
    readonly byteLength: number;
    readonly fields: readonly { readonly id: string; readonly unit: string }[];
  };
  readonly fields: readonly DatasetField[];
  readonly periods: readonly DatasetPeriod[];
}

export interface TerrainDataset {
  readonly fields: ReadonlyMap<string, Float64Array>;
  readonly centers: readonly (readonly [number, number, number])[];
}

export interface PeriodDataset {
  readonly descriptor: DatasetPeriod;
  readonly values: Float64Array;
}

export interface DatasetFiles {
  readText(relativePath: string): Promise<string>;
  readBytes(relativePath: string): Promise<ArrayBuffer>;
}

const headerBytes = 24;
const maximumDatasetBytes = 512 * 1024 * 1024;

function fail(message: string): never {
  throw new ProtocolError("invalid_visual_dataset", message);
}

function record(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

function finite(value: unknown, name: string): number {
  if (typeof value !== "number" || !Number.isFinite(value)) fail(`${name} must be finite`);
  return value;
}

function integer(value: unknown, name: string, minimum = 0): number {
  if (!Number.isSafeInteger(value) || (value as number) < minimum) fail(`${name} is invalid`);
  return value as number;
}

function text(value: unknown, name: string): string {
  if (typeof value !== "string" || value.length === 0) fail(`${name} is invalid`);
  return value;
}

export function validateDatasetPath(value: unknown, name = "dataset path"): string {
  const path = text(value, name);
  if (path.startsWith("/") || path.startsWith("\\") || /^[A-Za-z]:/.test(path) ||
      path.includes("\\") || path.split("/").some((part) => part === "" || part === "." || part === "..")) {
    fail(`${name} must be a safe relative path`);
  }
  return path;
}

function numberList(value: unknown, expectedLength: number, name: string): readonly number[] {
  if (!Array.isArray(value) || value.length !== expectedLength) fail(`${name} has the wrong length`);
  return Object.freeze(value.map((item, index) => finite(item, `${name}[${index}]`)));
}

export function parseVisualDatasetManifest(json: string): VisualDatasetManifest {
  let parsed: unknown;
  try { parsed = JSON.parse(json); } catch { fail("manifest is not valid JSON"); }
  if (!record(parsed)) fail("manifest must be an object");
  if (parsed.schemaVersion !== 1 || parsed.modelKind !== "dry_hydrostatic" ||
      parsed.arrayOrder !== "field-cell-level" ||
      parsed.averagingRule !== "accepted_step_end_time_weighted_rectangle") {
    fail("manifest contract is unsupported");
  }
  const cellsPerPanel = integer(parsed.cellsPerPanel, "cellsPerPanel", 1);
  const cellCount = integer(parsed.cellCount, "cellCount", 1);
  const levelCount = integer(parsed.levelCount, "levelCount", 1);
  if (cellCount !== 6 * cellsPerPanel * cellsPerPanel) fail("cellCount does not match the cubed sphere");
  if (!record(parsed.planet)) fail("planet is invalid");
  const planet = Object.freeze({ radiusM: finite(parsed.planet.radiusM, "radiusM"),
    gravityMS2: finite(parsed.planet.gravityMS2, "gravityMS2") });
  if (planet.radiusM <= 0 || planet.gravityMS2 <= 0) fail("planet values must be positive");
  const hybridAHalfPa = numberList(parsed.hybridAHalfPa, levelCount + 1, "hybridAHalfPa");
  const hybridBHalf = numberList(parsed.hybridBHalf, levelCount + 1, "hybridBHalf");
  if (!record(parsed.terrain) || !Array.isArray(parsed.terrain.fields)) fail("terrain descriptor is invalid");
  const terrainFields = parsed.terrain.fields.map((value, index) => {
    if (!record(value)) fail(`terrain field ${index} is invalid`);
    return Object.freeze({ id: text(value.id, "terrain field id"), unit: text(value.unit, "terrain field unit") });
  });
  const terrain = Object.freeze({ path: validateDatasetPath(parsed.terrain.path, "terrain path"),
    byteLength: integer(parsed.terrain.byteLength, "terrain byteLength", headerBytes),
    fields: Object.freeze(terrainFields) });
  if (terrain.byteLength > maximumDatasetBytes) fail("terrain file is too large");
  if (!Array.isArray(parsed.fields) || parsed.fields.length === 0) fail("fields are invalid");
  let nextOffset = 0;
  const identifiers = new Set<string>();
  const fields = parsed.fields.map((value, index) => {
    if (!record(value)) fail(`field ${index} is invalid`);
    const id = text(value.id, "field id");
    if (identifiers.has(id)) fail("field identifiers must be unique");
    identifiers.add(id);
    const location = value.location;
    if (location !== "surface" && location !== "atmosphere") fail(`field ${id} has invalid location`);
    const valueOffset = integer(value.valueOffset, `field ${id} offset`);
    const valueCount = integer(value.valueCount, `field ${id} count`, 1);
    const expected = location === "surface" ? cellCount : cellCount * levelCount;
    if (valueOffset !== nextOffset || valueCount !== expected) fail(`field ${id} has invalid shape or offset`);
    nextOffset += valueCount;
    return Object.freeze({ id, unit: text(value.unit, "field unit"), location, valueOffset, valueCount });
  });
  if (!Array.isArray(parsed.periods)) fail("periods are invalid");
  let previousIndex = -1;
  const periods = parsed.periods.map((value, index) => {
    if (!record(value)) fail(`period ${index} is invalid`);
    const descriptor: DatasetPeriod = Object.freeze({
      index: integer(value.index, "period index"),
      scheduledStartS: finite(value.scheduledStartS, "scheduledStartS"),
      scheduledEndS: finite(value.scheduledEndS, "scheduledEndS"),
      actualStartS: finite(value.actualStartS, "actualStartS"),
      actualEndS: finite(value.actualEndS, "actualEndS"),
      weightS: finite(value.weightS, "weightS"), coverage: finite(value.coverage, "coverage"),
      complete: value.complete === true, path: validateDatasetPath(value.path, "period path"),
      byteLength: integer(value.byteLength, "period byteLength", headerBytes),
    });
    if (descriptor.index <= previousIndex || descriptor.scheduledEndS <= descriptor.scheduledStartS ||
        descriptor.actualEndS <= descriptor.actualStartS || descriptor.weightS <= 0 ||
        descriptor.coverage <= 0 || descriptor.coverage > 1 || descriptor.byteLength > maximumDatasetBytes) {
      fail("period bounds, ordering, or size are invalid");
    }
    previousIndex = descriptor.index;
    return descriptor;
  });
  return Object.freeze({ schemaVersion: 1, modelKind: "dry_hydrostatic", cellsPerPanel,
    cellCount, levelCount, arrayOrder: "field-cell-level", planet, hybridAHalfPa,
    hybridBHalf, configFingerprint: text(parsed.configFingerprint, "configFingerprint"),
    terrainSourceFingerprint: text(parsed.terrainSourceFingerprint, "terrainSourceFingerprint"),
    averagingRule: "accepted_step_end_time_weighted_rectangle", terrain,
    fields: Object.freeze(fields), periods: Object.freeze(periods) });
}

function decodePayload(bytes: ArrayBuffer, expectedMagic: string, expectedFields: number,
                       expectedValues: number, expectedBytes: number): Float64Array {
  if (bytes.byteLength !== expectedBytes || bytes.byteLength !== headerBytes + 8 * expectedValues ||
      bytes.byteLength > maximumDatasetBytes) fail("binary byte length is invalid");
  const view = new DataView(bytes);
  const magic = new TextDecoder().decode(new Uint8Array(bytes, 0, 8));
  if (magic !== expectedMagic || view.getUint32(8, true) !== 1 ||
      view.getUint32(12, true) !== expectedFields ||
      Number(view.getBigUint64(16, true)) !== expectedValues) fail("binary header is invalid");
  const values = new Float64Array(expectedValues);
  for (let index = 0; index < expectedValues; index += 1) {
    values[index] = finite(view.getFloat64(headerBytes + 8 * index, true), "binary value");
  }
  return values;
}

export function decodeTerrain(manifest: VisualDatasetManifest, bytes: ArrayBuffer): TerrainDataset {
  const values = decodePayload(bytes, "MPSTERR1", manifest.terrain.fields.length,
    manifest.cellCount * manifest.terrain.fields.length, manifest.terrain.byteLength);
  const fields = new Map<string, Float64Array>();
  manifest.terrain.fields.forEach((field, index) => {
    if (fields.has(field.id)) fail("terrain field identifiers must be unique");
    fields.set(field.id, values.slice(index * manifest.cellCount, (index + 1) * manifest.cellCount));
  });
  const x = fields.get("center_x"), y = fields.get("center_y"), z = fields.get("center_z");
  if (!x || !y || !z) fail("terrain is missing cell centers");
  const centers = Object.freeze(Array.from({ length: manifest.cellCount }, (_, cell) => {
    const length = Math.hypot(x[cell], y[cell], z[cell]);
    if (Math.abs(length - 1) > 1e-10) fail("terrain cell center is not a unit vector");
    return Object.freeze([x[cell], y[cell], z[cell]] as const);
  }));
  return Object.freeze({ fields, centers });
}

export function decodePeriod(manifest: VisualDatasetManifest, descriptor: DatasetPeriod,
                             bytes: ArrayBuffer): PeriodDataset {
  const expectedValues = manifest.fields.at(-1)!.valueOffset + manifest.fields.at(-1)!.valueCount;
  return Object.freeze({ descriptor, values: decodePayload(bytes, "MPSMEAN1",
    manifest.fields.length, expectedValues, descriptor.byteLength) });
}

export function periodField(manifest: VisualDatasetManifest, period: PeriodDataset,
                            fieldId: string): Float64Array {
  const field = manifest.fields.find((candidate) => candidate.id === fieldId);
  if (!field) fail(`unknown field ${fieldId}`);
  return period.values.slice(field.valueOffset, field.valueOffset + field.valueCount);
}

export function periodLevelSlice(manifest: VisualDatasetManifest, period: PeriodDataset,
                                 fieldId: string, level: number): Float64Array {
  const field = manifest.fields.find((candidate) => candidate.id === fieldId);
  if (!field) fail(`unknown field ${fieldId}`);
  const values = periodField(manifest, period, fieldId);
  if (field.location === "surface") return values;
  if (!Number.isInteger(level) || level < 0 || level >= manifest.levelCount) fail("level is out of range");
  const slice = new Float64Array(manifest.cellCount);
  for (let cell = 0; cell < manifest.cellCount; cell += 1) slice[cell] = values[cell * manifest.levelCount + level];
  return slice;
}

export class VisualDatasetReader {
  readonly manifest: VisualDatasetManifest;
  readonly terrain: TerrainDataset;
  private readonly files: DatasetFiles;
  private readonly cache = new Map<number, PeriodDataset>();
  private readonly cacheSize: number;

  private constructor(files: DatasetFiles, manifest: VisualDatasetManifest,
                      terrain: TerrainDataset, cacheSize: number) {
    this.files = files; this.manifest = manifest; this.terrain = terrain; this.cacheSize = cacheSize;
  }

  static async open(files: DatasetFiles, cacheSize = 2): Promise<VisualDatasetReader> {
    if (!Number.isInteger(cacheSize) || cacheSize < 1 || cacheSize > 16) fail("cache size is invalid");
    const manifest = parseVisualDatasetManifest(await files.readText("manifest.json"));
    const terrain = decodeTerrain(manifest, await files.readBytes(manifest.terrain.path));
    return new VisualDatasetReader(files, manifest, terrain, cacheSize);
  }

  async period(index: number): Promise<PeriodDataset> {
    const cached = this.cache.get(index);
    if (cached) { this.cache.delete(index); this.cache.set(index, cached); return cached; }
    const descriptor = this.manifest.periods.find((period) => period.index === index);
    if (!descriptor) fail(`unknown period ${index}`);
    const decoded = decodePeriod(this.manifest, descriptor, await this.files.readBytes(descriptor.path));
    this.cache.set(index, decoded);
    while (this.cache.size > this.cacheSize) this.cache.delete(this.cache.keys().next().value!);
    return decoded;
  }
}
