import { describe, expect, it } from "vitest";
import { createGlobeGeometry, createGridGeometry, globeSubdivisionsPerCell, globeTrianglesPerCell, pickedCellFromFace, updateGlobeColors } from "./globe";
import { VisualDatasetV1 } from "@myplanetsim/protocol";

function dataset(cellsPerPanel: number, firstValue = 2): VisualDatasetV1 {
  const values = new Float64Array(6 * cellsPerPanel ** 2).fill(2);
  values[0] = firstValue;
  return {
    schemaVersion: 1,
    sourceKind: "visual_dataset",
    grid: { topology: "cubed_sphere", mapping: "equiangular_gnomonic_v1",
      cellsPerPanel, flattenOrder: "panel_major_then_j_then_i" },
    frame: { timeSeconds: null, step: null, configFingerprint: "fixture" },
    fields: [{ id: "depth", label: "Depth", unit: "m", kind: "scalar",
      provenance: "fixture", values }],
  };
}

describe("3D globe buffers", () => {
  it("creates finite one-geometry cell buffers and batched grid lines", () => {
    const fixture = dataset(2);
    const geometry = createGlobeGeometry(fixture, "depth");
    const subdivisions = globeSubdivisionsPerCell(2);
    expect(geometry.getAttribute("position").count).toBe(24 * (subdivisions + 1) ** 2);
    expect([...geometry.getAttribute("position").array].every(Number.isFinite)).toBe(true);
    expect(createGridGeometry(2, "panel_seams").getAttribute("position")?.count).toBeGreaterThan(48);
    expect(pickedCellFromFace(globeTrianglesPerCell(2) * 3 + 7, globeTrianglesPerCell(2))).toBe(3);
  });

  it("recolours an existing geometry in place when the grid is unchanged", () => {
    const geometry = createGlobeGeometry(dataset(2), "depth");
    const positions = geometry.getAttribute("position");
    const varied = dataset(2, 9);
    expect(updateGlobeColors(geometry, varied, "depth")).toBe(true);
    expect(geometry.getAttribute("position")).toBe(positions);
    const colors = geometry.getAttribute("color").array;
    expect([...colors.slice(0, 3)]).not.toEqual([...colors.slice(colors.length - 3)]);
    expect(updateGlobeColors(geometry, dataset(1), "depth")).toBe(false);
  });

  it("tessellates sparse panel cells onto the unit sphere", () => {
    const positions = createGlobeGeometry(dataset(2), "depth").getAttribute("position").array;
    for (let offset = 0; offset < positions.length; offset += 3) {
      expect(Math.hypot(positions[offset], positions[offset + 1], positions[offset + 2])).toBeCloseTo(1, 6);
    }
  });
});
