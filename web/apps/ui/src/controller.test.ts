import { describe, expect, it } from "vitest";
import { DeterministicMockSimulationClient } from "@myplanetsim/protocol/client";
import { SimulationController } from "./controller";

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
});
