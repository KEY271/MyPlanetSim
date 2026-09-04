import {
  EventV1,
  ProtocolError,
  RunBundleV1,
  RunRequestV1,
  validateEvent,
  validateRunRequest,
} from "./index";

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
  constructor(private readonly baseUrl: string, private readonly fetcher: typeof fetch = fetch) {}

  async capabilities(): Promise<CapabilitiesV1> {
    return this.json<CapabilitiesV1>("/api/v1/capabilities", { method: "GET" });
  }

  async submit(request: RunRequestV1): Promise<{ runId: string }> {
    validateRunRequest(request);
    return this.json("/api/v1/runs", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify(request),
    });
  }

  async *events(runId: string, afterSequence = -1): AsyncIterable<EventV1> {
    const response = await this.fetcher(`${this.baseUrl}/api/v1/runs/${encodeURIComponent(runId)}/events`, {
      headers: { accept: "application/x-ndjson", "x-after-sequence": String(afterSequence) },
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
      method: "POST",
    });
    if (!response.ok) throw new ProtocolError("http_error", `cancel failed: ${response.status}`);
  }

  async frame(runId: string, sequence: number): Promise<ArrayBuffer> {
    const response = await this.fetcher(`${this.baseUrl}/api/v1/runs/${encodeURIComponent(runId)}/frames/${sequence}`);
    if (!response.ok) throw new ProtocolError("http_error", `frame failed: ${response.status}`);
    return response.arrayBuffer();
  }

  async bundle(runId: string): Promise<RunBundleV1> {
    return this.json<RunBundleV1>(`/api/v1/runs/${encodeURIComponent(runId)}/bundle`, { method: "GET" });
  }

  private async json<T>(path: string, init: RequestInit): Promise<T> {
    const response = await this.fetcher(`${this.baseUrl}${path}`, init);
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
    if (!this.requests.has(runId)) throw new ProtocolError("unknown_run", runId);
    const event = (sequence: number, type: EventV1["type"]): EventV1 =>
      ({ protocolVersion: 1, runId, sequence, type } as EventV1);
    const events: EventV1[] = [
      event(0, "run.accepted"),
      event(1, "run.started"),
      { protocolVersion: 1, runId, sequence: 2, type: "frame.ready", frameSequence: 0,
        timeSeconds: 0, step: 0, relativePath: "frame_0.bin", byteLength: 64 },
      { protocolVersion: 1, runId, sequence: 3, type: "run.completed" },
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

  async frame(_runId: string, _sequence: number): Promise<ArrayBuffer> {
    const cellsPerPanel = 1; const cellCount = 6; const fingerprint = "mock";
    const bytes = new ArrayBuffer(52 + fingerprint.length + 32 * cellCount);
    const view = new DataView(bytes); const magic = "MPSFRAM1";
    for (let index = 0; index < magic.length; index += 1) view.setUint8(index, magic.charCodeAt(index));
    view.setUint32(8, 1, true); view.setUint32(12, 0, true); view.setBigUint64(16, BigInt(cellsPerPanel), true);
    view.setFloat64(24, 0, true); view.setBigUint64(32, 0n, true); view.setBigUint64(40, BigInt(cellCount), true); view.setUint32(48, fingerprint.length, true);
    new TextEncoder().encode(fingerprint).forEach((value, index) => view.setUint8(52 + index, value));
    const offset = 52 + fingerprint.length;
    for (let index = 0; index < cellCount; index += 1) view.setFloat64(offset + index * 8, 1, true);
    return bytes;
  }

  async bundle(runId: string): Promise<RunBundleV1> {
    return { protocolVersion: 1, runId, status: "completed", request: { protocolVersion: 1, presetId: "rest", run: { endTimeSeconds: 1, maximumTimeStepSeconds: 1, frameIntervalSteps: 1 }, initialCondition: { edits: [] } }, controlRequest: "", events: [], stderr: "" };
  }
}
