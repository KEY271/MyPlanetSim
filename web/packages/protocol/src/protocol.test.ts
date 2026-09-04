import { describe, expect, it } from "vitest";
import { ProtocolError, transitionRunState, validateRunRequest } from "./index";
import { DeterministicMockSimulationClient as Mock, HttpSimulationClient } from "./client";
import { addShallowWaterDerivedFields, fieldRange, parseShallowWaterCsv } from "./visual";
import { cubedSphereCells, gridEdges, hitTest, mapToUnitSphere } from "./geometry";

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

  it("authenticates HTTP gateway requests with the session token", async () => {
    const requests: RequestInit[] = [];
    const fetcher = (async (_input: RequestInfo | URL, init?: RequestInit) => {
      requests.push(init ?? {});
      return new Response(JSON.stringify({ protocolVersion: 1, presets: ["rest"], supportedEdits: ["gaussian_depth"] }),
        { status: 200, headers: { "content-type": "application/json" } });
    }) as typeof fetch;
    const client = new HttpSimulationClient("http://127.0.0.1:3000", "session-token", fetcher);
    await client.capabilities();
    expect(new Headers(requests[0].headers).get("authorization")).toBe("Bearer session-token");
    expect(() => new HttpSimulationClient("https://example.com", "session-token", fetcher))
      .toThrow(ProtocolError);
  });

  it("normalizes shuffled shallow-water CSV rows and derives speed", () => {
    const row = (panel: string, i: number, j: number) => `${panel},${i},${j},2,2,0,0`;
    const rows = [row("PY", 1, 1), row("PX", 0, 0)];
    for (const panel of ["PX", "PY", "NX", "NY", "PZ", "NZ"]) {
      for (let j = 0; j < 2; j += 1) for (let i = 0; i < 2; i += 1) {
        const candidate = row(panel, i, j);
        if (!rows.some((value) => value.startsWith(`${panel},${i},${j},`))) rows.push(candidate);
      }
    }
    const value = addShallowWaterDerivedFields(parseShallowWaterCsv(
      `panel,i,j,depth_m,momentum_x,momentum_y,momentum_z\n${rows.join("\n")}`, 2));
    expect(fieldRange(value.fields[0])).toEqual([2, 2]);
    expect(value.fields.find((field) => field.id === "speed_m_s")?.values[0]).toBe(1);
  });

  it("reconstructs panel mapping and cubed-sphere edge counts", () => {
    const point = mapToUnitSphere("PX", 0, 0);
    expect(point).toEqual([1, 0, 0]);
    for (const cellsPerPanel of [1, 2, 4, 16]) {
      expect(cubedSphereCells(cellsPerPanel)).toHaveLength(6 * cellsPerPanel ** 2);
      const edges = gridEdges(cellsPerPanel);
      expect(edges).toHaveLength(12 * cellsPerPanel ** 2);
      expect(edges.filter((edge) => edge.panelSeam)).toHaveLength(12 * cellsPerPanel);
      expect(hitTest(cubedSphereCells(cellsPerPanel)[0].center, cellsPerPanel)).toEqual({ panel: "PX", i: 0, j: 0 });
    }
  });
});
