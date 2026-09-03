import { describe, expect, it } from "vitest";
import { createGlobeGeometry, createGridGeometry, pickedCellFromFace } from "./globe";
import { parseShallowWaterCsv } from "@myplanetsim/protocol";

const csv = ["panel,i,j,depth_m,momentum_x,momentum_y,momentum_z",
  ...["PX", "PY", "NX", "NY", "PZ", "NZ"].flatMap((panel) =>
    [0, 1].flatMap((j) => [0, 1].map((i) => `${panel},${i},${j},2,0,0,0`)))].join("\n");

describe("3D globe buffers", () => {
  it("creates finite one-geometry cell buffers and batched grid lines", () => {
    const dataset = parseShallowWaterCsv(csv, 2);
    const geometry = createGlobeGeometry(dataset, "depth");
    expect(geometry.getAttribute("position").count).toBe(96);
    expect([...geometry.getAttribute("position").array].every(Number.isFinite)).toBe(true);
    expect(createGridGeometry(2, "panel_seams").getAttribute("position")?.count).toBe(48);
    expect(pickedCellFromFace(7)).toBe(3);
  });
});
