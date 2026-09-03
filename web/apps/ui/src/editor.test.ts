import { describe, expect, it } from "vitest";
import { appendEdit, createGaussianEdit, DraftRun, toRunRequest, undoEdit } from "./editor";

const draft: DraftRun = { presetId: "rest", endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 2, edits: [] };

describe("initial-condition editor", () => {
  it("uses the same unit-vector origin for a request and supports undo", () => {
    const edit = createGaussianEdit([0, 0, 1], 2, 0.2, "preserve_global");
    const edited = appendEdit(draft, edit);
    expect(toRunRequest(edited).initialCondition.edits[0].centerUnit).toEqual([0, 0, 1]);
    expect(undoEdit(edited).edits).toHaveLength(0);
  });
});
