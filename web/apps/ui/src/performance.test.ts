import { describe, expect, it } from "vitest";
import { cubedSphereCells, gridEdges } from "@myplanetsim/protocol";

describe("dataset scale gates", () => {
  it("keeps N=48/96 geometry structural counts bounded", () => {
    for (const cellsPerPanel of [48, 96]) {
      expect(cubedSphereCells(cellsPerPanel)).toHaveLength(6 * cellsPerPanel ** 2);
      expect(gridEdges(cellsPerPanel)).toHaveLength(12 * cellsPerPanel ** 2);
    }
  });
});
