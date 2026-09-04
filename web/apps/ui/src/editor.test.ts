import { describe, expect, it } from "vitest";
import { appendEdit, createGaussianEdit, DraftRun, toRunRequest, undoEdit } from "./editor";

const draft: DraftRun = { presetId: "rest", cellsPerPanel: 8, endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 2, edits: [] };

describe("initial-condition editor", () => {
  it("uses the same unit-vector origin for a request and supports undo", () => {
    const edit = createGaussianEdit([0, 0, 1], 2, 0.2, "preserve_global");
    const edited = appendEdit(draft, edit);
    expect(toRunRequest(edited).grid.cellsPerPanel).toBe(8);
    expect(toRunRequest(edited).initialCondition.edits[0].centerUnit).toEqual([0, 0, 1]);
    expect(undoEdit(edited).edits).toHaveLength(0);
  });

  it("sends a vertical resolution only when the draft sets one", () => {
    // A shallow-water draft must not carry a levels key at all, because the gateway
    // rejects an override for a preset with no column.
    expect("levels" in toRunRequest(draft).grid).toBe(false);
    expect(toRunRequest({ ...draft, levels: 16 }).grid.levels).toBe(16);
    expect(toRunRequest({ ...draft, levels: 1 }).grid.levels).toBe(1);
    expect(toRunRequest({ ...draft, levels: 30 }).grid.levels).toBe(30);
    for (const levels of [0, 31, 2.5, -1]) {
      expect(() => toRunRequest({ ...draft, levels })).toThrow(/levels is outside the supported range/);
    }
  });
});
