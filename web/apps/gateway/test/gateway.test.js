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
