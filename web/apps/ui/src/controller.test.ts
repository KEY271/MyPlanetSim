import { describe, expect, it } from "vitest";
import { DeterministicMockSimulationClient } from "@myplanetsim/protocol/client";
import { SimulationController } from "./controller";

describe("simulation controller", () => {
  it("runs a mock request and stores the authoritative frame", async () => {
    const controller = new SimulationController(new DeterministicMockSimulationClient());
    await controller.start({ protocolVersion: 1, presetId: "rest", run: { endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 1 }, initialCondition: { edits: [] } });
    expect(controller.value.state).toBe("completed");
    expect(controller.value.frames).toHaveLength(1);
    expect(controller.value.frames[0].sourceKind).toBe("live_run");
  });
});
