import { mkdtemp, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { test } from "node:test";
import assert from "node:assert/strict";
import { controlRequestText, createGatewayServer, validateRunRequest } from "../src/main.js";

const request = { protocolVersion: 1, presetId: "rest", run: { endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 2 }, initialCondition: { edits: [] } };

test("control translation is strict and path-free", () => {
  const text = controlRequestText(request, "run-1");
  assert.match(text, /control\.frame_directory = frames/);
  assert.doesNotMatch(text, /config|output|shell/);
  assert.throws(() => validateRunRequest({ ...request, presetId: "../../escape" }, new Map([["rest", "preset.cfg"]])));
});

test("gateway uses loopback auth and shell-free argv", async () => {
  const root = await mkdtemp(join(tmpdir(), "myplanetsim-gateway-"));
  const preset = join(root, "preset.cfg"); await writeFile(preset, "fixture");
  let spawnOptions; let spawnArgs;
  const fakeSpawn = (_binary, args, options) => {
    spawnArgs = args; spawnOptions = options;
    const listeners = new Map();
    const child = { stdout: null, stderr: null, on: (event, callback) => { listeners.set(event, callback); return child; }, kill: () => true };
    return child;
  };
  const gateway = createGatewayServer({ binary: "/configured/simulator", presets: { rest: preset }, runRoot: root, sessionToken: "token", spawn: fakeSpawn });
  await new Promise((resolve) => gateway.server.listen(0, "127.0.0.1", resolve));
  const port = gateway.server.address().port;
  const response = await fetch(`http://127.0.0.1:${port}/api/v1/runs`, { method: "POST", headers: { authorization: "Bearer token", host: `127.0.0.1:${port}`, origin: `http://127.0.0.1:${port}`, "content-type": "application/json" }, body: JSON.stringify(request) });
  assert.equal(response.status, 202); assert.equal(spawnOptions.shell, false); assert.deepEqual(spawnArgs.slice(-2), ["--event-stream", "ndjson"]);
  gateway.server.close();
});

test("gateway serves only published frames and supports reconnect sequence", async () => {
  const root = await mkdtemp(join(tmpdir(), "myplanetsim-gateway-"));
  const preset = join(root, "preset.cfg"); await writeFile(preset, "fixture");
  const listeners = new Map();
  const child = { stdout: { on: () => child.stdout }, stderr: null, on: (event, callback) => { listeners.set(event, callback); return child; }, kill: () => true };
  const gateway = createGatewayServer({ binary: "/configured/simulator", presets: { rest: preset }, runRoot: root, sessionToken: "token", spawn: () => child });
  await new Promise((resolve) => gateway.server.listen(0, "127.0.0.1", resolve));
  const port = gateway.server.address().port;
  const response = await fetch(`http://127.0.0.1:${port}/api/v1/runs`, { method: "POST", headers: { authorization: "Bearer token", origin: `http://127.0.0.1:${port}`, "content-type": "application/json" }, body: JSON.stringify(request) });
  const { runId } = await response.json();
  const run = gateway.runs.get(runId);
  const framePath = join(run.directory, "frames", "frame_0.bin"); await writeFile(framePath, Buffer.from([1, 2, 3]));
  child.stdout.on = (_event, callback) => { child.stdoutCallback = callback; return child.stdout; };
  // The fake child emits a complete event line through the registered stdout callback.
  run.pending = ""; run.lastSequence = -1;
  run.events.push({ protocolVersion: 1, runId, sequence: 0, type: "run.accepted" }); run.lastSequence = 0;
  run.events.push({ protocolVersion: 1, runId, sequence: 1, type: "frame.ready", frameSequence: 0, relativePath: "frame_0.bin", byteLength: 3, timeSeconds: 0, step: 0 }); run.lastSequence = 1;
  run.terminal = true;
  const events = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/events`, { headers: { authorization: "Bearer token", "x-after-sequence": "0" } });
  assert.equal(events.status, 200);
  const eventText = await events.text();
  assert.match(eventText, /frame.ready/);
  const frame = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/frames/0`, { headers: { authorization: "Bearer token" } });
  assert.deepEqual([...new Uint8Array(await frame.arrayBuffer())], [1, 2, 3]);
  gateway.server.closeAllConnections?.(); await gateway.shutdown();
});
