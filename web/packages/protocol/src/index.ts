export const protocolVersion = 1 as const;

export type MassPolicy = "preserve_global" | "allow_change";

export interface GaussianDepthEditV1 {
  id: string;
  kind: "gaussian_depth";
  centerUnit: readonly [number, number, number];
  amplitudeMeters: number;
  sigmaRadians: number;
  massPolicy: MassPolicy;
}

export interface RunRequestV1 {
  protocolVersion: typeof protocolVersion;
  presetId: string;
  grid: { cellsPerPanel: number };
  run: {
    endTimeSeconds: number;
    maximumTimeStepSeconds: number;
    frameIntervalSteps: number;
  };
  initialCondition: { edits: readonly GaussianDepthEditV1[] };
}

export type RunState =
  | "draft"
  | "submitting"
  | "running"
  | "cancelling"
  | "cancelled"
  | "completed"
  | "failed";

export interface EventBaseV1 {
  protocolVersion: typeof protocolVersion;
  runId: string;
  sequence: number;
}

export type EventV1 =
  | (EventBaseV1 & { type: "run.accepted" })
  | (EventBaseV1 & { type: "run.started" })
  | (EventBaseV1 & {
      type: "frame.ready";
      frameSequence: number;
      timeSeconds: number;
      step: number;
      relativePath: string;
      byteLength: number;
      frameSchemaVersion?: 1 | 2;
    })
  | (EventBaseV1 & { type: "diagnostics.sample"; timeSeconds: number; step: number; mass?: number; energy?: number; potentialEnstrophy?: number; axialAngularMomentum?: number; maximumCfl?: number })
  | (EventBaseV1 & { type: "run.completed" | "run.cancelled" })
  | (EventBaseV1 & { type: "run.failed"; code: string; message: string });

export interface FrameMetadataV1 {
  schemaVersion: 1;
  mapping: "equiangular_gnomonic_v1";
  cellsPerPanel: number;
  flattenOrder: "panel_major_then_j_then_i";
  timeSeconds: number;
  step: number;
  configFingerprint: string;
  fields: readonly ["depth", "momentum_x", "momentum_y", "momentum_z"];
}

export class ProtocolError extends Error {
  readonly code: string;

  constructor(code: string, message: string) {
    super(message);
    this.name = "ProtocolError";
    this.code = code;
  }
}

function record(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

function finite(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value);
}

function require(condition: boolean, message: string): asserts condition {
  if (!condition) throw new ProtocolError("invalid_contract", message);
}

// Every published frame is retained by the gateway on disk and decoded into browser memory,
// so the run budget is bounded before the native process is started.
export const maxPublishedFrames = 512;

export function estimatedFrameCount(run: RunRequestV1["run"]): number {
  const steps = Math.ceil(run.endTimeSeconds / run.maximumTimeStepSeconds);
  return Math.floor(steps / run.frameIntervalSteps) + 2;
}

export function validateRunRequest(value: unknown): RunRequestV1 {
  require(record(value), "request must be an object");
  const object = value as Record<string, unknown>;
  require(object.protocolVersion === protocolVersion, "unsupported protocol version");
  require(typeof object.presetId === "string" && object.presetId.length > 0, "presetId is required");
  require(record(object.grid) && record(object.run) && record(object.initialCondition), "request sections are required");
  const grid = object.grid as Record<string, unknown>;
  require(Number.isInteger(grid.cellsPerPanel) && (grid.cellsPerPanel as number) >= 1 &&
    (grid.cellsPerPanel as number) <= 96, "cellsPerPanel is outside the supported range");
  const run = object.run as Record<string, unknown>;
  require(finite(run.endTimeSeconds) && run.endTimeSeconds > 0 && run.endTimeSeconds <= 31536000,
    "endTimeSeconds is outside the supported range");
  require(finite(run.maximumTimeStepSeconds) && run.maximumTimeStepSeconds > 0 &&
    run.maximumTimeStepSeconds <= 86400, "maximumTimeStepSeconds is outside the supported range");
  require(Number.isInteger(run.frameIntervalSteps) && (run.frameIntervalSteps as number) >= 1 &&
    (run.frameIntervalSteps as number) <= 1000000, "frameIntervalSteps is outside the supported range");
  require(estimatedFrameCount(run as RunRequestV1["run"]) <= maxPublishedFrames,
    `the request would publish more than ${maxPublishedFrames} frames; raise the frame interval or shorten the run`);
  const initialCondition = object.initialCondition as Record<string, unknown>;
  const edits = initialCondition.edits;
  require(Array.isArray(edits) && edits.length <= 64, "initialCondition.edits is invalid");
  for (const edit of edits) {
    require(record(edit) && edit.kind === "gaussian_depth", "unsupported initial edit");
    require(typeof edit.id === "string" && edit.id.length > 0, "edit id is required");
    require(Array.isArray(edit.centerUnit) && edit.centerUnit.length === 3 &&
      edit.centerUnit.every(finite), "edit center must be a finite 3-vector");
    const norm = Math.hypot(...edit.centerUnit);
    require(Math.abs(norm - 1) <= 1e-12, "edit center must be a unit vector");
    require(finite(edit.amplitudeMeters) && Math.abs(edit.amplitudeMeters) <= 1e6,
      "edit amplitude is outside the supported range");
    require(finite(edit.sigmaRadians) && edit.sigmaRadians > 1e-6 &&
      edit.sigmaRadians <= Math.PI, "edit sigma is outside the supported range");
    require(edit.massPolicy === "preserve_global" || edit.massPolicy === "allow_change",
      "unsupported edit mass policy");
  }
  return value as unknown as RunRequestV1;
}

export function validateEvent(value: unknown): EventV1 {
  require(record(value), "event must be an object");
  const object = value as Record<string, unknown>;
  require(object.protocolVersion === protocolVersion, "unsupported event version");
  require(typeof object.runId === "string" && object.runId.length > 0, "event runId is required");
  require(Number.isInteger(object.sequence) && (object.sequence as number) >= 0, "event sequence is invalid");
  require(typeof object.type === "string" && object.type.length > 0, "event type is required");
  if (object.type === "frame.ready") {
    require(Number.isInteger(object.frameSequence) && (object.frameSequence as number) >= 0 &&
      finite(object.timeSeconds) && Number.isInteger(object.step) && (object.step as number) >= 0 &&
      typeof object.relativePath === "string" && !(object.relativePath as string).includes("..") &&
      Number.isInteger(object.byteLength) && (object.byteLength as number) >= 0, "frame event is invalid");
  }
  if (object.type === "run.failed") {
    require(typeof object.code === "string" && typeof object.message === "string",
      "failure event is invalid");
  }
  return value as unknown as EventV1;
}

const transitions: Record<RunState, readonly RunState[]> = {
  draft: ["submitting"],
  submitting: ["running", "failed", "cancelled"],
  running: ["cancelling", "completed", "failed", "cancelled"],
  cancelling: ["cancelled", "failed"],
  cancelled: [],
  completed: [],
  failed: [],
};

export function transitionRunState(state: RunState, next: RunState): RunState {
  if (!transitions[state].includes(next)) {
    throw new ProtocolError("invalid_state_transition", `${state} cannot transition to ${next}`);
  }
  return next;
}

export interface RunBundleV1 {
  readonly protocolVersion: typeof protocolVersion;
  readonly runId: string;
  readonly status: RunState;
  readonly request: RunRequestV1;
  readonly controlRequest: string;
  readonly events: readonly EventV1[];
  readonly stderr: string;
  readonly buildMetadata?: Record<string, string>;
}

export function diagnosticEvents(events: readonly EventV1[]) {
  return events.filter((event): event is Extract<EventV1, { type: "diagnostics.sample" }> => event.type === "diagnostics.sample");
}

export * from "./visual";
export * from "./geometry";
