import { describe, expect, it } from "vitest";
import { decodePeriod, decodeTerrain, parseVisualDatasetManifest, periodLevelSlice,
  ProtocolError, VisualDatasetReader } from "./index";

const cellsPerPanel = 1, cellCount = 6, levelCount = 2;
const fields = [
  { id: "surface_pressure", unit: "Pa", location: "surface", valueOffset: 0, valueCount: 6 },
  { id: "temperature", unit: "K", location: "atmosphere", valueOffset: 6, valueCount: 12 },
];
const period = { index: 0, scheduledStartS: 0, scheduledEndS: 60, actualStartS: 0,
  actualEndS: 60, weightS: 60, coverage: 1, complete: true,
  path: "means/period_000000.bin", byteLength: 24 + 18 * 8 };
const manifestValue = { schemaVersion: 1, modelKind: "dry_hydrostatic", cellsPerPanel,
  cellCount, levelCount, arrayOrder: "field-cell-level", planet: { radiusM: 2, gravityMS2: 10 },
  hybridAHalfPa: [1000, 500, 0], hybridBHalf: [0, .5, 1], configFingerprint: "fp",
  terrainSourceFingerprint: "terrain-fp", averagingRule: "accepted_step_end_time_weighted_rectangle",
  terrain: { path: "terrain.bin", byteLength: 24 + 5 * cellCount * 8, fields: [
    { id: "center_x", unit: "1" }, { id: "center_y", unit: "1" }, { id: "center_z", unit: "1" },
    { id: "elevation", unit: "m" }, { id: "surface_geopotential", unit: "m2 s-2" }] },
  fields, periods: [period] };

function binary(magic: string, fieldCount: number, values: readonly number[]): ArrayBuffer {
  const bytes = new ArrayBuffer(24 + values.length * 8), view = new DataView(bytes);
  new Uint8Array(bytes, 0, 8).set(new TextEncoder().encode(magic));
  view.setUint32(8, 1, true); view.setUint32(12, fieldCount, true);
  view.setBigUint64(16, BigInt(values.length), true);
  values.forEach((value, index) => view.setFloat64(24 + 8 * index, value, true));
  return bytes;
}

const terrainValues = [1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1,
  ...Array(12).fill(0)];
const meanValues = Array.from({ length: 18 }, (_, index) => index + 1);

describe("period-mean dataset", () => {
  it("decodes terrain and cell-major level fields", () => {
    const manifest = parseVisualDatasetManifest(JSON.stringify(manifestValue));
    const terrain = decodeTerrain(manifest, binary("MPSTERR1", 5, terrainValues));
    const mean = decodePeriod(manifest, manifest.periods[0], binary("MPSMEAN1", 2, meanValues));
    expect(terrain.centers[3]).toEqual([-1, 0, 0]);
    expect([...periodLevelSlice(manifest, mean, "temperature", 1)]).toEqual([8, 10, 12, 14, 16, 18]);
    expect([...periodLevelSlice(manifest, mean, "surface_pressure", 999)]).toEqual([1, 2, 3, 4, 5, 6]);
  });

  it("loads lazily and rejects missing or corrupt files", async () => {
    const files = new Map<string, string | ArrayBuffer>([["manifest.json", JSON.stringify(manifestValue)],
      ["terrain.bin", binary("MPSTERR1", 5, terrainValues)],
      [period.path, binary("MPSMEAN1", 2, meanValues)]]);
    const reader = await VisualDatasetReader.open({
      readText: async (path) => files.get(path) as string,
      readBytes: async (path) => files.get(path) as ArrayBuffer,
    });
    expect((await reader.period(0)).values[17]).toBe(18);
    await expect(reader.period(1)).rejects.toThrow(ProtocolError);
    const corrupt = binary("MPSMEAN9", 2, meanValues);
    expect(() => decodePeriod(reader.manifest, reader.manifest.periods[0], corrupt)).toThrow(/header/);
    expect(() => parseVisualDatasetManifest(JSON.stringify({ ...manifestValue,
      terrain: { ...manifestValue.terrain, path: "../terrain.bin" } }))).toThrow(/safe relative/);
  });
});
