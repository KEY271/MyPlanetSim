import { describe, expect, it } from "vitest";
import { ProtocolError, transitionRunState, validateRunRequest } from "./index";
import { DeterministicMockSimulationClient as Mock } from "./client";

const request = {
  protocolVersion: 1 as const,
  presetId: "rest",
  run: { endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 1 },
  initialCondition: { edits: [] },
};

describe("protocol contracts", () => {
  it("validates a request and rejects a version mismatch", () => {
    expect(validateRunRequest(request)).toEqual(request);
    expect(() => validateRunRequest({ ...request, protocolVersion: 2 })).toThrow(ProtocolError);
  });

  it("enforces lifecycle transitions", () => {
    expect(transitionRunState("draft", "submitting")).toBe("submitting");
    expect(() => transitionRunState("completed", "running")).toThrow(ProtocolError);
  });

  it("emits deterministic mock events and supports cancellation", async () => {
    const client = new Mock();
    const { runId } = await client.submit(request);
    await client.cancel(runId);
    const events = [];
    for await (const event of client.events(runId)) events.push(event.type);
    expect(events).toEqual(["run.accepted", "run.started", "frame.ready", "run.cancelled"]);
  });
});
