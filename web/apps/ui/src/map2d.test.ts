import { describe, expect, it } from "vitest";
import { geoPath } from "d3-geo";
import { cubedSphereCells } from "@myplanetsim/protocol";
import { equirectangularProjection, geoCellPolygon, inverseProject, pickMapCell, projectUnit } from "./projection";

describe("2D map projection", () => {
  it("round trips known unit vectors", () => {
    for (const point of [[1, 0, 0], [0, 0, 1], [0, -1, 0]] as const) {
      const projected = projectUnit(point, 800, 400); const restored = inverseProject(projected[0], projected[1], 800, 400);
      expect(restored[0]).toBeCloseTo(point[0], 12); expect(restored[1]).toBeCloseTo(point[1], 12); expect(restored[2]).toBeCloseTo(point[2], 12);
    }
  });

  it("picks the containing cell instead of requiring proximity to its center", () => {
    const point = [1, 0, 0] as const;
    const projected = projectUnit(point, 800, 400);
    expect(pickMapCell(projected[0], projected[1], 800, 400, 2))
      .toEqual({ panel: "PX", i: 1, j: 1 });
  });

  it("clips cells that cross the antimeridian into separate projected pieces", () => {
    const cell = cubedSphereCells(1).find((candidate) => candidate.panel === "NX")!;
    const path = geoPath(equirectangularProjection(800, 400))(geoCellPolygon(cell.corners)) ?? "";
    expect(path.match(/M/g)?.length).toBeGreaterThan(1);
  });
});
