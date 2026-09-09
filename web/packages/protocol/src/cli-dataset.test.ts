import { describe, expect, it } from "vitest";
import { readFile } from "node:fs/promises";
import { join } from "node:path";
import { periodLevelSlice, VisualDatasetReader } from "./index";

const datasetRoot = process.env.MPS_VISUAL_DATASET;

describe.skipIf(!datasetRoot)("C++ visual dataset", () => {
  it("decodes the checked CLI fixture without a browser", async () => {
    const root = datasetRoot!;
    const reader = await VisualDatasetReader.open({
      readText: async (relativePath) => readFile(join(root, relativePath), "utf8"),
      readBytes: async (relativePath) => {
        const bytes = await readFile(join(root, relativePath));
        return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
      },
    });

    expect(reader.manifest.cellsPerPanel).toBe(4);
    expect(reader.manifest.levelCount).toBe(8);
    expect(reader.manifest.periods).toHaveLength(2);
    expect(reader.manifest.periods.every((period) => period.complete)).toBe(true);
    expect([...reader.terrain.fields.get("elevation")!].every((value) => value === 0)).toBe(true);

    const period = await reader.period(0);
    const temperature = periodLevelSlice(reader.manifest, period, "temperature", 0);
    expect(temperature).toHaveLength(reader.manifest.cellCount);
    expect([...temperature].every(Number.isFinite)).toBe(true);
  });
});
