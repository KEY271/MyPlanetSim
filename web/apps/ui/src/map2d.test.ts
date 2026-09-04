import { describe, expect, it } from "vitest";
import { geoContains, geoPath } from "d3-geo";
import { cubedSphereCells, panelNames } from "@myplanetsim/protocol";
import { equirectangularProjection, geoCellPolygon, inverseProject, isPointInsideProjectedMap, pickMapCell, projectUnit, unitToCellIndex, unitToLongitudeLatitude } from "./projection";

describe("2D map projection", () => {
  it("round trips known unit vectors", () => {
    for (const point of [[1, 0, 0], [0, 0, 1], [0, -1, 0]] as const) {
      const projected = projectUnit(point, 800, 400); const restored = inverseProject(projected[0], projected[1], 800, 400);
      expect(restored[0]).toBeCloseTo(point[0], 12); expect(restored[1]).toBeCloseTo(point[1], 12); expect(restored[2]).toBeCloseTo(point[2], 12);
    }
  });

  it("preserves an arbitrary clicked location without snapping to a cell center", () => {
    const point = [0.8123, -0.3154, 0.4901] as const;
    const length = Math.hypot(...point);
    const normalized: readonly [number, number, number] = [point[0] / length, point[1] / length, point[2] / length];
    const projected = projectUnit(normalized, 930, 700, 1.4, [21, -13]);
    const restored = inverseProject(projected[0], projected[1], 930, 700, 1.4, [21, -13]);
    expect(restored[0]).toBeCloseTo(normalized[0], 12);
    expect(restored[1]).toBeCloseTo(normalized[1], 12);
    expect(restored[2]).toBeCloseTo(normalized[2], 12);
  });

  it("fits the whole globe inside either a wide or tall viewport", () => {
    for (const [width, height] of [[1000, 300], [300, 1000]]) {
      for (const point of [[1, 0, 0], [-1, 0, 0], [0, 0, 1], [0, 0, -1]] as const) {
        const [x, y] = projectUnit(point, width, height);
        expect(x).toBeGreaterThanOrEqual(0); expect(x).toBeLessThanOrEqual(width);
        expect(y).toBeGreaterThanOrEqual(0); expect(y).toBeLessThanOrEqual(height);
      }
    }
  });

  it("rejects clicks in the canvas margins outside the projected map", () => {
    expect(isPointInsideProjectedMap(500, 500, 1000, 1000)).toBe(true);
    expect(isPointInsideProjectedMap(500, 100, 1000, 1000)).toBe(false);
    expect(isPointInsideProjectedMap(500, 900, 1000, 1000)).toBe(false);
    expect(isPointInsideProjectedMap(500, 100, 1000, 1000, 2, [0, 0])).toBe(true);
  });

  it("picks the containing cell instead of requiring proximity to its center", () => {
    const point = [1, 0, 0] as const;
    const projected = projectUnit(point, 800, 400);
    expect(pickMapCell(projected[0], projected[1], 800, 400, 2))
      .toEqual({ panel: "PX", i: 1, j: 1 });
  });

  it("matches the shared hit test in the allocation-free raster lookup", () => {
    for (const cellsPerPanel of [1, 4, 7]) {
      for (let y = 3; y < 400; y += 7) for (let x = 3; x < 800; x += 11) {
        const point = inverseProject(x + 0.5, y + 0.5, 800, 400);
        const picked = pickMapCell(x + 0.5, y + 0.5, 800, 400, cellsPerPanel);
        const expected = (panelNames.indexOf(picked.panel) * cellsPerPanel + picked.j) * cellsPerPanel + picked.i;
        expect(unitToCellIndex(point[0], point[1], point[2], cellsPerPanel), `${x},${y},N=${cellsPerPanel}`).toBe(expected);
      }
    }
  });

  it("clips cells that cross the antimeridian into separate projected pieces", () => {
    const cell = cubedSphereCells(1).find((candidate) => candidate.panel === "NX")!;
    const path = geoPath(equirectangularProjection(800, 400))(geoCellPolygon(cell.corners)) ?? "";
    expect(path.match(/M/g)?.length).toBeGreaterThan(1);
  });

  it("uses the small spherical polygon containing each cubed-sphere cell center", () => {
    for (const cell of cubedSphereCells(4)) {
      expect(geoContains(geoCellPolygon(cell.corners), unitToLongitudeLatitude(cell.center)), `${cell.panel}/${cell.i}/${cell.j}`).toBe(true);
    }
  });

  it("covers the sphere exactly once away from cell boundaries", () => {
    const polygons = cubedSphereCells(4).map((cell) => geoCellPolygon(cell.corners));
    for (let latitude = -80; latitude <= 80; latitude += 20) for (let longitude = -170; longitude < 180; longitude += 20) {
      expect(polygons.filter((polygon) => geoContains(polygon, [longitude + 0.37, latitude + 0.23])).length).toBe(1);
    }
  });

  it("projects a cubed-sphere edge as a curved spherical path", () => {
    const cell = cubedSphereCells(2).find((candidate) => candidate.panel === "PZ" && candidate.i === 0 && candidate.j === 0)!;
    const path = geoPath(equirectangularProjection(800, 400))({
      type: "LineString",
      coordinates: [unitToLongitudeLatitude(cell.corners[0]), unitToLongitudeLatitude(cell.corners[1])],
    }) ?? "";
    expect(path.match(/L/g)?.length).toBeGreaterThan(1);
  });
});
