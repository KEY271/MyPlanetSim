import { addShallowWaterDerivedFields, decodeFrameV1, EventV1, ProtocolError, RunBundleV1, RunRequestV1, RunState, transitionRunState } from "@myplanetsim/protocol";
import { SimulationClient } from "@myplanetsim/protocol/client";

export interface ControllerSnapshot {
  readonly state: RunState;
  readonly runId: string | null;
  readonly events: readonly EventV1[];
  readonly frames: readonly ReturnType<typeof decodeFrameV1>[];
  readonly currentFrame: number;
  readonly error: string | null;
  readonly request: RunRequestV1 | null;
}

export class SimulationController {
  private snapshot: ControllerSnapshot = { state: "draft", runId: null, events: [], frames: [], currentFrame: 0, error: null, request: null };
  private listeners = new Set<(snapshot: ControllerSnapshot) => void>();
  private client: SimulationClient;

  constructor(client: SimulationClient) { this.client = client; }
  get value() { return this.snapshot; }
  subscribe(listener: (snapshot: ControllerSnapshot) => void) { this.listeners.add(listener); return () => this.listeners.delete(listener); }
  private update(next: Partial<ControllerSnapshot>) { this.snapshot = { ...this.snapshot, ...next }; for (const listener of this.listeners) listener(this.snapshot); }

  async start(request: RunRequestV1) {
    if (this.snapshot.state !== "draft" && this.snapshot.state !== "completed" && this.snapshot.state !== "failed" && this.snapshot.state !== "cancelled") throw new ProtocolError("run_in_progress", "a run is already active");
    this.update({ state: "submitting", events: [], frames: [], currentFrame: 0, error: null, request });
    try {
      const submitted = await this.client.submit(request); this.update({ runId: submitted.runId });
      await this.consume(submitted.runId, -1);
    } catch (error) { this.update({ state: "failed", error: error instanceof Error ? error.message : "run failed" }); }
  }

  private async consume(runId: string, afterSequence: number) {
    try {
      for await (const event of this.client.events(runId, afterSequence)) {
        const nextState = event.type === "run.started" ? "running" : event.type === "run.completed" ? "completed" : event.type === "run.cancelled" ? "cancelled" : event.type === "run.failed" ? "failed" : this.snapshot.state;
        if (nextState !== this.snapshot.state) this.update({ state: transitionRunState(this.snapshot.state, nextState) });
        this.update({ events: [...this.snapshot.events, event] });
        if (event.type === "frame.ready") {
          const frame = addShallowWaterDerivedFields(decodeFrameV1(await this.client.frame(runId, event.frameSequence)));
          this.update({ frames: [...this.snapshot.frames, frame], currentFrame: this.snapshot.frames.length });
        }
        if (event.type === "run.failed") throw new ProtocolError(event.code, event.message);
      }
    } catch (error) { this.update({ state: "failed", error: error instanceof Error ? error.message : "event stream failed" }); }
  }

  async reconnect() { if (this.snapshot.runId) await this.consume(this.snapshot.runId, this.snapshot.events.at(-1)?.sequence ?? -1); }
  async cancel() { if (this.snapshot.runId && this.snapshot.state === "running") { this.update({ state: transitionRunState("running", "cancelling") }); await this.client.cancel(this.snapshot.runId); } }
  setCurrentFrame(index: number) { if (index >= 0 && index < this.snapshot.frames.length) this.update({ currentFrame: index }); }
  async bundle(): Promise<RunBundleV1 | null> { return this.snapshot.runId ? this.client.bundle(this.snapshot.runId) : null; }
}
