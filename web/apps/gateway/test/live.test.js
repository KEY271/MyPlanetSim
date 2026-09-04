import { existsSync } from "node:fs";
import { mkdtemp } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { test } from "node:test";
import assert from "node:assert/strict";
import { startGateway } from "../src/main.js";

test("native run applies the requested N and produces authoritative frames", { skip: process.env.MPS_ENABLE_LIVE !== "1" }, async () => {
  const root = resolve(process.cwd(), "../../..");
  const binary = resolve(root, "build/dev/my_planet_sim");
  const preset = resolve(root, "configs/phase3_rest_n4.cfg");
  if (!existsSync(binary) || !existsSync(preset)) return;
  const runRoot = await mkdtemp(join(tmpdir(), "myplanetsim-live-"));
  const gateway = await startGateway({ binary, presets: { rest: preset }, runRoot, sessionToken: "live-token", port: 0 });
  try {
    const port = gateway.server.address().port;
    const headers = { authorization: "Bearer live-token", origin: `http://127.0.0.1:${port}`, "content-type": "application/json" };
    const submit = (run) => fetch(`http://127.0.0.1:${port}/api/v1/runs`, { method: "POST", headers, body: JSON.stringify({ protocolVersion: 1, presetId: "rest", grid: { cellsPerPanel: 6 }, run, initialCondition: { edits: [{ id: "wave", kind: "gaussian_depth", centerUnit: [0, 0, 1], amplitudeMeters: 10, sigmaRadians: 0.2, massPolicy: "preserve_global" }] } }) });
    const waitForTerminal = async (runId) => {
      let status = "running";
      for (let attempt = 0; attempt < 100 && (status === "running" || status === "cancelling"); attempt += 1) {
        await new Promise((resolveSleep) => setTimeout(resolveSleep, 100));
        const statusResponse = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}`, { headers });
        status = (await statusResponse.json()).status;
      }
      return status;
    };

    const response = await submit({ endTimeSeconds: 120, maximumTimeStepSeconds: 60, frameIntervalSteps: 1 });
    assert.equal(response.status, 202); const { runId } = await response.json();
    assert.equal(await waitForTerminal(runId), "completed");
    const bundleResponse = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/bundle`, { headers });
    const events = (await bundleResponse.json()).events;
    assert.deepEqual(events.map((event) => event.type).filter((type) => type === "run.accepted" || type === "run.started" || type === "run.completed"), ["run.accepted", "run.started", "run.completed"]);
    const frame = events.find((event) => event.type === "frame.ready" && event.frameSequence === 0); assert.ok(frame);
    const frameResponse = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/frames/0`, { headers });
    assert.equal(frameResponse.status, 200); const bytes = new Uint8Array(await frameResponse.arrayBuffer());
    assert.deepEqual([...bytes.slice(0, 8)], [...Buffer.from("MPSFRAM1")]);
    assert.equal(new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength).getBigUint64(16, true), 6n);
    const laterFrame = events.find((event) => event.type === "frame.ready" && event.frameSequence === 1); assert.ok(laterFrame);
    const laterFrameResponse = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/frames/1`, { headers });
    assert.equal(laterFrameResponse.status, 200); const laterBytes = new Uint8Array(await laterFrameResponse.arrayBuffer());
    assert.notDeepEqual([...laterBytes], [...bytes]);

    const cancelResponse = await submit({ endTimeSeconds: 31536000, maximumTimeStepSeconds: 60, frameIntervalSteps: 1000000 });
    assert.equal(cancelResponse.status, 202); const { runId: cancelledRunId } = await cancelResponse.json();
    const cancel = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${cancelledRunId}/cancel`, { method: "POST", headers });
    assert.equal(cancel.status, 202);
    assert.equal(await waitForTerminal(cancelledRunId), "cancelled");
  } finally {
    await gateway.shutdown();
  }
});
