import { mkdtemp, writeFile } from "node:fs/promises";
import { EventEmitter } from "node:events";
import { request as httpRequest } from "node:http";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { PassThrough } from "node:stream";
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
  const unauthorized = await fetch(`http://127.0.0.1:${port}/api/v1/capabilities`);
  assert.equal(unauthorized.status, 401);
  await unauthorized.arrayBuffer();
  const preflight = await fetch(`http://127.0.0.1:${port}/api/v1/runs`, { method: "OPTIONS", headers: { origin: "http://localhost:5173", "access-control-request-method": "POST", "access-control-request-headers": "authorization,content-type" } });
  assert.equal(preflight.status, 204);
  assert.equal(preflight.headers.get("access-control-allow-origin"), "http://localhost:5173");
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
  const bundle = await fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/bundle`, { headers: { authorization: "Bearer token" } });
  assert.equal((await bundle.json()).request.presetId, "rest");
  gateway.server.closeAllConnections?.(); await gateway.shutdown();
});

test("gateway keeps a live event stream subscribed through the terminal event", async () => {
  const root = await mkdtemp(join(tmpdir(), "myplanetsim-gateway-"));
  const preset = join(root, "preset.cfg"); await writeFile(preset, "fixture");
  const child = new EventEmitter(); child.stdout = new PassThrough(); child.stderr = new PassThrough(); child.kill = () => true;
  const gateway = createGatewayServer({ binary: "/configured/simulator", presets: { rest: preset }, runRoot: root, sessionToken: "token", spawn: () => child });
  await new Promise((resolve) => gateway.server.listen(0, "127.0.0.1", resolve));
  const port = gateway.server.address().port;
  const headers = { authorization: "Bearer token", "content-type": "application/json" };
  const response = await fetch(`http://127.0.0.1:${port}/api/v1/runs`, { method: "POST", headers, body: JSON.stringify(request) });
  const { runId } = await response.json();
  const streamResponsePromise = fetch(`http://127.0.0.1:${port}/api/v1/runs/${runId}/events`, { headers });
  await new Promise((resolve) => setImmediate(resolve));
  const emit = (sequence, type) => child.stdout.write(`${JSON.stringify({ protocolVersion: 1, runId, sequence, type })}\n`);
  emit(0, "run.accepted");
  const streamResponse = await streamResponsePromise; const reader = streamResponse.body.getReader(); const decoder = new TextDecoder();
  assert.match(decoder.decode((await reader.read()).value), /run\.accepted/);
  emit(1, "run.started");
  assert.match(decoder.decode((await reader.read()).value), /run\.started/);
  emit(2, "run.completed");
  assert.match(decoder.decode((await reader.read()).value), /run\.completed/);
  assert.equal((await reader.read()).done, true);
  child.emit("close", 0, null);
  await gateway.shutdown();
});

test("gateway reserves its single run slot before asynchronous startup", async () => {
  const root = await mkdtemp(join(tmpdir(), "myplanetsim-gateway-"));
  const preset = join(root, "preset.cfg"); await writeFile(preset, "fixture");
  const children = [];
  const gateway = createGatewayServer({ binary: "/configured/simulator", presets: { rest: preset }, runRoot: root, sessionToken: "token", spawn: () => {
    const child = new EventEmitter(); child.stdout = null; child.stderr = null; child.kill = () => true; children.push(child); return child;
  } });
  await new Promise((resolve) => gateway.server.listen(0, "127.0.0.1", resolve));
  const port = gateway.server.address().port; const payload = JSON.stringify(request); const split = Math.floor(payload.length / 2);
  const startRequest = () => {
    let outgoing;
    const response = new Promise((resolveResponse, reject) => {
      outgoing = httpRequest({ hostname: "127.0.0.1", port, path: "/api/v1/runs", method: "POST", headers: { authorization: "Bearer token", "content-type": "application/json", "content-length": Buffer.byteLength(payload) } }, (incoming) => {
        incoming.resume(); incoming.on("end", () => resolveResponse(incoming.statusCode));
      });
      outgoing.on("error", reject); outgoing.write(payload.slice(0, split));
    });
    return { response, outgoing };
  };
  const first = startRequest(); const second = startRequest();
  await new Promise((resolve) => setTimeout(resolve, 10));
  first.outgoing.end(payload.slice(split)); second.outgoing.end(payload.slice(split));
  assert.deepEqual((await Promise.all([first.response, second.response])).sort(), [202, 409]);
  assert.equal(children.length, 1);
  children[0].emit("close", 0, null);
  await gateway.shutdown();
});
