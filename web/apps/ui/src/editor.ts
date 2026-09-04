import { GaussianDepthEditV1, RunRequestV1, validateRunRequest, UnitVector } from "@myplanetsim/protocol";

export interface DraftRun {
  readonly presetId: string;
  readonly cellsPerPanel: number;
  readonly endTimeSeconds: number;
  readonly maximumTimeStepSeconds: number;
  readonly frameIntervalSteps: number;
  readonly edits: readonly GaussianDepthEditV1[];
}

export function createGaussianEdit(centerUnit: UnitVector, amplitudeMeters: number, sigmaRadians: number, massPolicy: GaussianDepthEditV1["massPolicy"], id = "edit-1"): GaussianDepthEditV1 {
  return { id, kind: "gaussian_depth", centerUnit: [...centerUnit] as [number, number, number], amplitudeMeters, sigmaRadians, massPolicy };
}

export function appendEdit(draft: DraftRun, edit: GaussianDepthEditV1): DraftRun {
  return { ...draft, edits: [...draft.edits, edit] };
}

export function undoEdit(draft: DraftRun): DraftRun { return { ...draft, edits: draft.edits.slice(0, -1) }; }
export function clearEdits(draft: DraftRun): DraftRun { return { ...draft, edits: [] }; }

export function toRunRequest(draft: DraftRun): RunRequestV1 {
  const request: RunRequestV1 = { protocolVersion: 1, presetId: draft.presetId,
    grid: { cellsPerPanel: draft.cellsPerPanel },
    run: { endTimeSeconds: draft.endTimeSeconds, maximumTimeStepSeconds: draft.maximumTimeStepSeconds, frameIntervalSteps: draft.frameIntervalSteps },
    initialCondition: { edits: draft.edits.map((edit) => ({ ...edit, centerUnit: [...edit.centerUnit] as [number, number, number] })) } };
  return validateRunRequest(request);
}
