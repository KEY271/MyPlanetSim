import { describe, expect, it } from "vitest";
import { EventV1, RunBundleV1, RunRequestV1 } from "@myplanetsim/protocol";
import { CapabilitiesV1, DeterministicMockSimulationClient, SimulationClient } from "@myplanetsim/protocol/client";
import { SimulationController } from "./controller";

const request: RunRequestV1 = { protocolVersion: 1, presetId: "dry", grid: { cellsPerPanel: 2 }, run: { endTimeSeconds: 20, maximumTimeStepSeconds: 10, frameIntervalSteps: 1 }, initialCondition: { edits: [] } };

// A dry hydrostatic gateway stands in for the native binary: it publishes FrameV2 bytes
// and reports the frame schema on the lifecycle event.
class DryHydrostaticStubClient implements SimulationClient {
  readonly levels = 4;
  private readonly cells = 24;

  async capabilities(): Promise<CapabilitiesV1> {
    return { protocolVersion: 1, presets: ["dry"], supportedEdits: ["gaussian_depth"],
      presetDetails: [{ id: "dry", modelKind: "dry_hydrostatic", frameSchemaVersion: 2, levels: this.levels, supportedEdits: [], maximumCellsPerPanel: 24 }] };
  }

  async submit(): Promise<{ runId: string }> { return { runId: "dry-1" }; }
  async cancel(): Promise<void> {}
  async bundle(): Promise<RunBundleV1> { return { protocolVersion: 1, runId: "dry-1", status: "completed", request, controlRequest: "", events: [], stderr: "" }; }

  async *events(): AsyncIterable<EventV1> {
    yield { protocolVersion: 1, runId: "dry-1", sequence: 0, type: "run.accepted" };
    yield { protocolVersion: 1, runId: "dry-1", sequence: 1, type: "run.started" };
    for (const frameSequence of [0, 1]) {
      yield { protocolVersion: 1, runId: "dry-1", sequence: 2 + frameSequence, type: "frame.ready", frameSequence, timeSeconds: 10 * frameSequence, step: frameSequence, relativePath: `frame_${frameSequence}.bin`, byteLength: 0, frameSchemaVersion: 2 };
    }
    yield { protocolVersion: 1, runId: "dry-1", sequence: 4, type: "run.completed" };
  }

  async frame(_runId: string, sequence: number): Promise<ArrayBuffer> {
    const volume = this.cells * this.levels;
    const fingerprint = new TextEncoder().encode("stub");
    const bytes = new ArrayBuffer(60 + fingerprint.length + 8 * (this.cells + 7 * volume));
    const view = new DataView(bytes);
    new Uint8Array(bytes, 0, 8).set(new TextEncoder().encode("MPSFRAM2"));
    view.setUint32(8, 2, true);
    view.setBigUint64(16, 2n, true);
    view.setBigUint64(24, BigInt(this.levels), true);
    view.setFloat64(32, 10 * sequence, true);
    view.setBigUint64(40, BigInt(sequence), true);
    view.setBigUint64(48, BigInt(this.cells), true);
    view.setUint32(56, fingerprint.length, true);
    new Uint8Array(bytes, 60, fingerprint.length).set(fingerprint);
    let offset = 60 + fingerprint.length;
    for (let index = 0; index < this.cells + 7 * volume; index += 1) {
      view.setFloat64(offset, 1 + sequence + index / 1000, true);
      offset += 8;
    }
    return bytes;
  }
}

describe("simulation controller", () => {
  it("runs a mock request and stores the authoritative frame", async () => {
    const controller = new SimulationController(new DeterministicMockSimulationClient());
    await controller.start({ protocolVersion: 1, presetId: "rest", grid: { cellsPerPanel: 2 }, run: { endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 2 }, initialCondition: { edits: [{ id: "wave", kind: "gaussian_depth", centerUnit: [1, 0, 0], amplitudeMeters: 100, sigmaRadians: 0.2, massPolicy: "preserve_global" }] } });
    expect(controller.value.state).toBe("completed");
    expect(controller.value.frames).toHaveLength(6);
    const frames = controller.value.frames.filter((frame) => frame.schemaVersion === 1);
    expect(frames).toHaveLength(6);
    expect(frames.every((frame) => frame.grid.cellsPerPanel === 2)).toBe(true);
    expect(frames[0].sourceKind).toBe("live_run");
    expect(controller.value.currentFrame).toBe(controller.value.frames.length - 1);
    expect([...frames[0].fields.find((field) => field.id === "depth_anomaly_m")!.values])
      .toEqual(new Array(24).fill(0));
    expect(Math.max(...frames.at(-1)!.fields.find((field) => field.id === "speed_m_s")!.values)).toBeGreaterThan(0);
    expect([...frames.at(-1)!.fields.find((field) => field.id === "depth_anomaly_m")!.values].some((value) => Math.abs(value) > 0)).toBe(true);
  });

  it("stores dry hydrostatic frames without applying shallow-water derivation", async () => {
    const controller = new SimulationController(new DryHydrostaticStubClient());
    await controller.start(request);
    expect(controller.value.state).toBe("completed");
    expect(controller.value.frames).toHaveLength(2);
    // A FrameV2 has no depth field, so deriving speed or a depth anomaly from it would
    // throw; the controller must keep the decoded frame as published.
    const [first, second] = controller.value.frames;
    expect(first.schemaVersion).toBe(2);
    expect(second.schemaVersion).toBe(2);
    if (first.schemaVersion !== 2 || second.schemaVersion !== 2) throw new Error("expected FrameV2");
    expect(first.levels).toBe(4);
    expect(first.timeSeconds).toBe(0);
    expect(second.timeSeconds).toBe(10);
    expect(second.step).toBe(1);
    expect(first.configFingerprint).toBe("stub");
    expect(controller.value.currentFrame).toBe(1);
    const capabilities = await controller.capabilities();
    expect(capabilities.presetDetails?.[0].modelKind).toBe("dry_hydrostatic");
  });
});
