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

test("native dry hydrostatic runs publish decodable FrameV2 levels and columns", { skip: process.env.MPS_ENABLE_LIVE !== "1" }, async () => {
  const root = resolve(process.cwd(), "../../..");
  const binary = resolve(root, "build/dev/my_planet_sim");
  const preset = resolve(root, "configs/phase5_visualizer_rest_n4.cfg");
  if (!existsSync(binary) || !existsSync(preset)) return;
  const levels = 8;
  const runRoot = await mkdtemp(join(tmpdir(), "myplanetsim-live-v2-"));
  const gateway = await startGateway({
    binary, runRoot, sessionToken: "live-token", port: 0,
    presets: { rest: resolve(root, "configs/phase3_rest_n4.cfg"),
      dry: { path: preset, modelKind: "dry_hydrostatic", frameSchemaVersion: 2, levels, supportedEdits: [], maximumCellsPerPanel: 8 } },
  });
  try {
    const port = gateway.server.address().port;
    const headers = { authorization: "Bearer live-token", origin: `http://127.0.0.1:${port}`, "content-type": "application/json" };
    const submit = (body) => fetch(`http://127.0.0.1:${port}/api/v1/runs`, { method: "POST", headers, body: JSON.stringify(body) });
    const waitForTerminal = async (runId) => {
      let status = "running";
      for (let attempt = 0; attempt < 100 && (status === "running" || status === "cancelling"); attempt += 1) {
        await new Promise((resolveSleep) => setTimeout(resolveSleep, 100));
        status = (await (await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}`, { headers })).json()).status;
      }
      return status;
    };

    const capabilities = await (await fetch(`http://127.0.0.1:${port}/api/v1/capabilities`, { headers })).json();
    const descriptor = capabilities.presetDetails.find((value) => value.id === "dry");
    assert.deepEqual(descriptor, { id: "dry", modelKind: "dry_hydrostatic", frameSchemaVersion: 2, levels, maximumLevels: 30, supportedEdits: [], maximumCellsPerPanel: 8 });

    const dryRequest = { protocolVersion: 1, presetId: "dry", grid: { cellsPerPanel: 4 }, run: { endTimeSeconds: 30, maximumTimeStepSeconds: 10, frameIntervalSteps: 1 }, initialCondition: { edits: [] } };
    // An edit that the shallow-water preset accepts must be rejected for this preset.
    const rejected = await submit({ ...dryRequest, initialCondition: { edits: [{ id: "wave", kind: "gaussian_depth", centerUnit: [0, 0, 1], amplitudeMeters: 10, sigmaRadians: 0.2, massPolicy: "preserve_global" }] } });
    assert.equal(rejected.status, 400);

    const response = await submit(dryRequest);
    assert.equal(response.status, 202); const { runId } = await response.json();
    assert.equal(await waitForTerminal(runId), "completed");
    const events = (await (await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/bundle`, { headers })).json()).events;
    const frames = events.filter((event) => event.type === "frame.ready");
    assert.ok(frames.length >= 2);
    assert.ok(frames.every((event) => event.frameSchemaVersion === 2));
    assert.equal(events.at(-1).type, "run.completed");

    const decode = async (frameSequence) => {
      const frameResponse = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/frames/${frameSequence}`, { headers });
      assert.equal(frameResponse.status, 200);
      const bytes = new Uint8Array(await frameResponse.arrayBuffer());
      assert.deepEqual([...bytes.slice(0, 8)], [...Buffer.from("MPSFRAM2")]);
      const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
      const cellsPerPanel = Number(view.getBigUint64(16, true));
      const frameLevels = Number(view.getBigUint64(24, true));
      const cells = Number(view.getBigUint64(48, true));
      const fingerprintLength = view.getUint32(56, true);
      assert.equal(cellsPerPanel, 4);
      assert.equal(frameLevels, levels);
      assert.equal(cells, 6 * cellsPerPanel ** 2);
      // The exact fixed payload the browser decoder expects: 8 * C * (1 + 7K) after the header.
      assert.equal(bytes.byteLength, 60 + fingerprintLength + 8 * cells * (1 + 7 * frameLevels));
      const payload = 60 + fingerprintLength;
      const value = (field, cell, level) => view.getFloat64(payload + 8 * (cells + field * cells * frameLevels + cell * frameLevels + level), true);
      return { cells, levels: frameLevels, timeSeconds: view.getFloat64(32, true), step: Number(view.getBigUint64(40, true)),
        surfacePressurePa: (cell) => view.getFloat64(payload + 8 * cell, true), value };
    };

    const first = await decode(0);
    assert.equal(first.timeSeconds, 0);
    assert.equal(first.step, 0);
    for (let cell = 0; cell < first.cells; cell += 1) {
      assert.equal(first.surfacePressurePa(cell), 100000);
      let previous = 0;
      for (let level = 0; level < first.levels; level += 1) {
        const pressurePa = first.value(0, cell, level);
        assert.ok(pressurePa > previous, "pressure must increase downward through the column");
        previous = pressurePa;
        assert.ok(Math.abs(first.value(2, cell, level) - 288) < 1, "the rest preset is isothermal at 288 K");
        assert.equal(first.value(3, cell, level), 0);
      }
    }
    const last = await decode(frames.length - 1);
    assert.ok(last.timeSeconds > 0);
    assert.equal(last.step, frames.length - 1);

    // ADR 0007: the same preset re-resolved vertically, checked end to end because the
    // coordinate is generated by the solver rather than written in the configuration.
    const refined = await submit({ ...dryRequest, grid: { cellsPerPanel: 4, levels: 20 } });
    assert.equal(refined.status, 202); const { runId: refinedRunId } = await refined.json();
    assert.equal(await waitForTerminal(refinedRunId), "completed");
    const refinedFrame = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${refinedRunId}/frames/0`, { headers });
    const refinedBytes = new Uint8Array(await refinedFrame.arrayBuffer());
    const refinedView = new DataView(refinedBytes.buffer, refinedBytes.byteOffset, refinedBytes.byteLength);
    assert.equal(Number(refinedView.getBigUint64(24, true)), 20);
    const refinedFingerprintLength = refinedView.getUint32(56, true);
    assert.equal(refinedBytes.byteLength, 60 + refinedFingerprintLength + 8 * 96 * (1 + 7 * 20));
    const refinedFingerprint = new TextDecoder().decode(refinedBytes.subarray(60, 60 + refinedFingerprintLength));
    const baseFingerprintLength = new DataView((await (await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/frames/0`, { headers })).arrayBuffer())).getUint32(56, true);
    const baseBytes = new Uint8Array(await (await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/frames/0`, { headers })).arrayBuffer());
    const baseFingerprint = new TextDecoder().decode(baseBytes.subarray(60, 60 + baseFingerprintLength));
    // A different vertical resolution is a different configuration, not the preset's.
    assert.notEqual(refinedFingerprint, baseFingerprint);
    // The generated coordinate keeps the preset's model top and still increases downward.
    let previous = 0;
    for (let level = 0; level < 20; level += 1) {
      const pressurePa = refinedView.getFloat64(60 + refinedFingerprintLength + 8 * (96 + level), true);
      assert.ok(pressurePa > previous);
      previous = pressurePa;
    }
    // Beyond the interactive bound the gateway refuses before the solver is started.
    assert.equal((await submit({ ...dryRequest, grid: { cellsPerPanel: 4, levels: 31 } })).status, 400);
    assert.equal((await submit({ ...dryRequest, grid: { cellsPerPanel: 4, levels: 0 } })).status, 400);

    const cancelResponse = await submit({ ...dryRequest, run: { endTimeSeconds: 31536000, maximumTimeStepSeconds: 10, frameIntervalSteps: 1000000 } });
    assert.equal(cancelResponse.status, 202); const { runId: cancelledRunId } = await cancelResponse.json();
    assert.equal((await fetch(`http://127.0.0.1:${port}/api/v1/runs/${cancelledRunId}/cancel`, { method: "POST", headers })).status, 202);
    assert.equal(await waitForTerminal(cancelledRunId), "cancelled");
    // A cancelled run must leave no partial frame behind: every announced frame is servable.
    const cancelledEvents = (await (await fetch(`http://127.0.0.1:${port}/api/v1/runs/${cancelledRunId}/bundle`, { headers })).json()).events;
    for (const event of cancelledEvents.filter((value) => value.type === "frame.ready")) {
      const frameResponse = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${cancelledRunId}/frames/${event.frameSequence}`, { headers });
      assert.equal(frameResponse.status, 200);
      assert.equal((await frameResponse.arrayBuffer()).byteLength, event.byteLength);
    }
  } finally {
    await gateway.shutdown();
  }
});
