import { describe, expect, it } from "vitest";
import { inverseProject, projectUnit } from "./projection";

describe("2D map projection", () => {
  it("round trips known unit vectors", () => {
    for (const point of [[1, 0, 0], [0, 0, 1], [0, -1, 0]] as const) {
      const projected = projectUnit(point, 800, 400); const restored = inverseProject(projected[0], projected[1], 800, 400);
      expect(restored[0]).toBeCloseTo(point[0], 12); expect(restored[1]).toBeCloseTo(point[1], 12); expect(restored[2]).toBeCloseTo(point[2], 12);
    }
  });
});
