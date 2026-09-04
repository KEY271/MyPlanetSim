import { describe, expect, it } from "vitest";
import { createGlobeGeometry, createGridGeometry, globeSubdivisionsPerCell, globeTrianglesPerCell, pickedCellFromFace, updateGlobeColors } from "./globe";
import { parseShallowWaterCsv } from "@myplanetsim/protocol";

const csv = ["panel,i,j,depth_m,momentum_x,momentum_y,momentum_z",
  ...["PX", "PY", "NX", "NY", "PZ", "NZ"].flatMap((panel) =>
    [0, 1].flatMap((j) => [0, 1].map((i) => `${panel},${i},${j},2,0,0,0`)))].join("\n");

describe("3D globe buffers", () => {
  it("creates finite one-geometry cell buffers and batched grid lines", () => {
    const dataset = parseShallowWaterCsv(csv, 2);
    const geometry = createGlobeGeometry(dataset, "depth");
    const subdivisions = globeSubdivisionsPerCell(2);
    expect(geometry.getAttribute("position").count).toBe(24 * (subdivisions + 1) ** 2);
    expect([...geometry.getAttribute("position").array].every(Number.isFinite)).toBe(true);
    expect(createGridGeometry(2, "panel_seams").getAttribute("position")?.count).toBeGreaterThan(48);
    expect(pickedCellFromFace(globeTrianglesPerCell(2) * 3 + 7, globeTrianglesPerCell(2))).toBe(3);
  });

  it("recolours an existing geometry in place when the grid is unchanged", () => {
    const geometry = createGlobeGeometry(parseShallowWaterCsv(csv, 2), "depth");
    const positions = geometry.getAttribute("position");
    const varied = parseShallowWaterCsv(csv.replace("PX,0,0,2", "PX,0,0,9"), 2);
    expect(updateGlobeColors(geometry, varied, "depth")).toBe(true);
    expect(geometry.getAttribute("position")).toBe(positions);
    const colors = geometry.getAttribute("color").array;
    expect([...colors.slice(0, 3)]).not.toEqual([...colors.slice(colors.length - 3)]);
    const coarse = ["panel,i,j,depth_m,momentum_x,momentum_y,momentum_z",
      ...["PX", "PY", "NX", "NY", "PZ", "NZ"].map((panel) => `${panel},0,0,2,0,0,0`)].join("\n");
    expect(updateGlobeColors(geometry, parseShallowWaterCsv(coarse, 1), "depth")).toBe(false);
  });

  it("tessellates sparse panel cells onto the unit sphere", () => {
    const dataset = parseShallowWaterCsv(csv, 2);
    const positions = createGlobeGeometry(dataset, "depth").getAttribute("position").array;
    for (let offset = 0; offset < positions.length; offset += 3) {
      expect(Math.hypot(positions[offset], positions[offset + 1], positions[offset + 2])).toBeCloseTo(1, 6);
    }
  });
});
