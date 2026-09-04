import {
  EventV1,
  ProtocolError,
  RunBundleV1,
  RunRequestV1,
  validateEvent,
  validateRunRequest,
} from "./index";
import { cubedSphereCells, UnitVector } from "./geometry";

export interface CapabilitiesV1 {
  protocolVersion: 1;
  presets: readonly string[];
  supportedEdits: readonly ["gaussian_depth"];
}

export interface SimulationClient {
  capabilities(): Promise<CapabilitiesV1>;
  submit(request: RunRequestV1): Promise<{ runId: string }>;
  events(runId: string, afterSequence?: number): AsyncIterable<EventV1>;
  frame(runId: string, sequence: number): Promise<ArrayBuffer>;
  bundle(runId: string): Promise<RunBundleV1>;
  cancel(runId: string): Promise<void>;
}

export class HttpSimulationClient implements SimulationClient {
  private readonly baseUrl: string;

  constructor(
    baseUrl: string,
    private readonly sessionToken: string,
    private readonly fetcher: typeof fetch = fetch,
  ) {
    if (!sessionToken) throw new ProtocolError("missing_token", "a gateway session token is required");
    const parsed = new URL(baseUrl);
    if (parsed.protocol !== "http:" ||
        (parsed.hostname !== "127.0.0.1" && parsed.hostname !== "localhost" && parsed.hostname !== "[::1]")) {
      throw new ProtocolError("invalid_gateway", "the gateway URL must use HTTP on a loopback host");
    }
    this.baseUrl = parsed.href.replace(/\/$/, "");
  }

  private headers(values: HeadersInit = {}): Headers {
    const headers = new Headers(values);
    headers.set("authorization", `Bearer ${this.sessionToken}`);
    return headers;
  }

  async capabilities(): Promise<CapabilitiesV1> {
    return this.json<CapabilitiesV1>("/api/v1/capabilities", { method: "GET" });
  }

  async submit(request: RunRequestV1): Promise<{ runId: string }> {
    validateRunRequest(request);
    return this.json("/api/v1/runs", {
      method: "POST",
      headers: this.headers({ "content-type": "application/json" }),
      body: JSON.stringify(request),
    });
  }

  async *events(runId: string, afterSequence = -1): AsyncIterable<EventV1> {
    const response = await this.fetcher(`${this.baseUrl}/api/v1/runs/${encodeURIComponent(runId)}/events`, {
      headers: this.headers({ accept: "application/x-ndjson", "x-after-sequence": String(afterSequence) }),
    });
    if (!response.ok || !response.body) throw new ProtocolError("http_error", `event stream failed: ${response.status}`);
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let pending = "";
    for (;;) {
      const chunk = await reader.read();
      pending += decoder.decode(chunk.value ?? new Uint8Array(), { stream: !chunk.done });
      const lines = pending.split("\n");
      pending = lines.pop() ?? "";
      for (const line of lines) {
        if (line.trim()) yield validateEvent(JSON.parse(line));
      }
      if (chunk.done) break;
    }
    if (pending.trim()) yield validateEvent(JSON.parse(pending));
  }

  async cancel(runId: string): Promise<void> {
    const response = await this.fetcher(`${this.baseUrl}/api/v1/runs/${encodeURIComponent(runId)}/cancel`, {
      method: "POST", headers: this.headers(),
    });
    if (!response.ok) throw new ProtocolError("http_error", `cancel failed: ${response.status}`);
  }

  async frame(runId: string, sequence: number): Promise<ArrayBuffer> {
    const response = await this.fetcher(`${this.baseUrl}/api/v1/runs/${encodeURIComponent(runId)}/frames/${sequence}`, { headers: this.headers() });
    if (!response.ok) throw new ProtocolError("http_error", `frame failed: ${response.status}`);
    return response.arrayBuffer();
  }

  async bundle(runId: string): Promise<RunBundleV1> {
    return this.json<RunBundleV1>(`/api/v1/runs/${encodeURIComponent(runId)}/bundle`, { method: "GET" });
  }

  private async json<T>(path: string, init: RequestInit): Promise<T> {
    const response = await this.fetcher(`${this.baseUrl}${path}`, { ...init, headers: this.headers(init.headers) });
    if (!response.ok) throw new ProtocolError("http_error", `request failed: ${response.status}`);
    return (await response.json()) as T;
  }
}

export class DeterministicMockSimulationClient implements SimulationClient {
  private readonly requests = new Map<string, RunRequestV1>();
  private readonly cancelled = new Set<string>();
  private nextRun = 1;

  async capabilities(): Promise<CapabilitiesV1> {
    return { protocolVersion: 1, presets: ["rest"], supportedEdits: ["gaussian_depth"] };
  }

  async submit(request: RunRequestV1): Promise<{ runId: string }> {
    validateRunRequest(request);
    const runId = `mock-${this.nextRun++}`;
    this.requests.set(runId, request);
    return { runId };
  }

  async *events(runId: string, afterSequence = -1): AsyncIterable<EventV1> {
    const request = this.requests.get(runId);
    if (!request) throw new ProtocolError("unknown_run", runId);
    const event = (sequence: number, type: EventV1["type"]): EventV1 =>
      ({ protocolVersion: 1, runId, sequence, type } as EventV1);
    const steps = mockFrameSteps(request);
    const fingerprint = `mock-${runId}`;
    const byteLength = 52 + fingerprint.length + 32 * 6 * request.grid.cellsPerPanel ** 2;
    const events: EventV1[] = [
      event(0, "run.accepted"),
      event(1, "run.started"),
      ...steps.map((step, frameSequence): EventV1 => ({
        protocolVersion: 1, runId, sequence: frameSequence + 2, type: "frame.ready", frameSequence,
        timeSeconds: mockFrameTime(request, step), step,
        relativePath: `frame_${frameSequence}.bin`, byteLength,
      })),
      event(steps.length + 2, "run.completed"),
    ];
    for (const value of events) {
      if (value.sequence <= afterSequence) continue;
      yield this.cancelled.has(runId) && value.type === "run.completed"
        ? { protocolVersion: 1, runId, sequence: value.sequence, type: "run.cancelled" }
        : value;
    }
  }

  async cancel(runId: string): Promise<void> {
    if (!this.requests.has(runId)) throw new ProtocolError("unknown_run", runId);
    this.cancelled.add(runId);
  }

  async frame(runId: string, sequence: number): Promise<ArrayBuffer> {
    const request = this.requests.get(runId);
    if (!request) throw new ProtocolError("unknown_run", runId);
    const steps = mockFrameSteps(request);
    const step = steps[sequence];
    if (step === undefined) throw new ProtocolError("unknown_frame", String(sequence));
    const cellsPerPanel = request.grid.cellsPerPanel; const cellCount = 6 * cellsPerPanel ** 2; const fingerprint = `mock-${runId}`;
    const bytes = new ArrayBuffer(52 + fingerprint.length + 32 * cellCount);
    const view = new DataView(bytes); const magic = "MPSFRAM1";
    for (let index = 0; index < magic.length; index += 1) view.setUint8(index, magic.charCodeAt(index));
    view.setUint32(8, 1, true); view.setUint32(12, 0, true); view.setBigUint64(16, BigInt(cellsPerPanel), true);
    view.setFloat64(24, mockFrameTime(request, step), true); view.setBigUint64(32, BigInt(step), true); view.setBigUint64(40, BigInt(cellCount), true); view.setUint32(48, fingerprint.length, true);
    new TextEncoder().encode(fingerprint).forEach((value, index) => view.setUint8(52 + index, value));
    const offset = 52 + fingerprint.length;
    const fields = mockFields(cellsPerPanel, request, step);
    for (let field = 0; field < fields.length; field += 1) for (let cell = 0; cell < cellCount; cell += 1) {
      view.setFloat64(offset + (field * cellCount + cell) * 8, fields[field][cell], true);
    }
    return bytes;
  }

  async bundle(runId: string): Promise<RunBundleV1> {
    const request = this.requests.get(runId);
    if (!request) throw new ProtocolError("unknown_run", runId);
    return { protocolVersion: 1, runId, status: this.cancelled.has(runId) ? "cancelled" : "completed", request, controlRequest: "mock", events: [], stderr: "" };
  }
}

function mockFrameSteps(request: RunRequestV1): number[] {
  const totalSteps = Math.max(1, Math.ceil(request.run.endTimeSeconds / request.run.maximumTimeStepSeconds));
  const interval = request.run.frameIntervalSteps;
  const requestedFrameCount = Math.floor(totalSteps / interval) + 1;
  const stride = requestedFrameCount <= 121 ? interval : interval * Math.ceil(requestedFrameCount / 121);
  const steps: number[] = [];
  for (let step = 0; step <= totalSteps; step += stride) steps.push(step);
  if (steps.at(-1) !== totalSteps) steps.push(totalSteps);
  return steps;
}

function mockFrameTime(request: RunRequestV1, step: number): number {
  return Math.min(request.run.endTimeSeconds, step * request.run.maximumTimeStepSeconds);
}

function dot(first: UnitVector, second: UnitVector): number {
  return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

function mockFields(cellsPerPanel: number, request: RunRequestV1, step: number): readonly Float64Array[] {
  const cells = cubedSphereCells(cellsPerPanel);
  const depth = new Float64Array(cells.length).fill(3000);
  const momentumX = new Float64Array(cells.length);
  const momentumY = new Float64Array(cells.length);
  const momentumZ = new Float64Array(cells.length);
  const elapsedTimeSeconds = mockFrameTime(request, step);
  const waveTravel = elapsedTimeSeconds * Math.PI / (3 * 3600);
  const waveCycle = waveTravel % (2 * Math.PI);
  const waveRadius = waveCycle <= Math.PI ? waveCycle : 2 * Math.PI - waveCycle;
  const velocityRamp = 1 - Math.exp(-elapsedTimeSeconds / 600);
  for (const edit of request.initialCondition.edits) {
    const anomaly = new Float64Array(cells.length);
    for (let cell = 0; cell < cells.length; cell += 1) {
      const center = cells[cell].center;
      const cosine = Math.max(-1, Math.min(1, dot(edit.centerUnit, center)));
      const distance = Math.acos(cosine);
      anomaly[cell] = edit.amplitudeMeters * Math.exp(-((distance - waveRadius) ** 2) / (2 * edit.sigmaRadians ** 2));
      if (elapsedTimeSeconds > 0) {
        const tangent = [edit.centerUnit[0] - cosine * center[0], edit.centerUnit[1] - cosine * center[1], edit.centerUnit[2] - cosine * center[2]] as const;
        const tangentLength = Math.hypot(...tangent);
        if (tangentLength > 1e-12) {
          const speed = Math.min(25, Math.abs(edit.amplitudeMeters) * 0.04) * velocityRamp * Math.exp(-((distance - waveRadius) ** 2) / (2 * (edit.sigmaRadians * 1.5) ** 2));
          momentumX[cell] -= 3000 * speed * tangent[0] / tangentLength;
          momentumY[cell] -= 3000 * speed * tangent[1] / tangentLength;
          momentumZ[cell] -= 3000 * speed * tangent[2] / tangentLength;
        }
      }
    }
    const mean = edit.massPolicy === "preserve_global" ? anomaly.reduce((sum, value) => sum + value, 0) / anomaly.length : 0;
    for (let cell = 0; cell < cells.length; cell += 1) depth[cell] += anomaly[cell] - mean;
  }
  return [depth, momentumX, momentumY, momentumZ];
}
