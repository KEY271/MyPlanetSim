import { randomBytes } from "node:crypto";
import { createServer } from "node:http";
import { readFileSync } from "node:fs";
import { lstat, mkdir, mkdtemp, readFile, realpath, rm, writeFile } from "node:fs/promises";
import { spawn as defaultSpawn } from "node:child_process";
import { basename, join, resolve } from "node:path";

const protocolVersion = 1;
const maxBodyBytes = 64 * 1024;
// Mirrors @myplanetsim/protocol maxPublishedFrames. The native solver may shrink its time step
// below the requested maximum, so the budget is enforced again while the run streams events.
export const maxPublishedFrames = 512;
export const maxPublishedBytes = 256 * 1024 * 1024;
export const maxRetainedRuns = 8;
// Mirrors kControlMaxLevels in the C++ control contract and the FrameV2 allowlist.
export const maxInteractiveLevels = 30;

export function estimatedFrameCount(run) {
  return Math.floor(Math.ceil(run.endTimeSeconds / run.maximumTimeStepSeconds) / run.frameIntervalSteps) + 2;
}

function protocolError(message) { const error = new Error(message); error.code = "invalid_request"; return error; }

function finite(value) { return typeof value === "number" && Number.isFinite(value); }

// A preset is configured either as a bare path, which keeps the Phase 3 shallow-water
// contract, or as a descriptor that states its model kind, frame schema, levels, edits,
// and interactive N limit.
export function presetDescriptor(preset, id) {
  const described = typeof preset === "string" ? { path: preset } : preset;
  const descriptor = { id, modelKind: "shallow_water", frameSchemaVersion: 1, levels: null,
    supportedEdits: ["gaussian_depth"], maximumCellsPerPanel: 96, ...described };
  // ADR 0007: a dry preset may be re-resolved vertically. `levels` stays the preset default
  // and `maximumLevels` bounds what a request may ask for.
  if (descriptor.modelKind === "dry_hydrostatic" && descriptor.maximumLevels === undefined) {
    descriptor.maximumLevels = maxInteractiveLevels;
  }
  return descriptor;
}

// FrameV2 stores surface pressure plus seven volume fields as little-endian float64, so the
// published size of one frame is exact and can be checked before the run starts. A run that
// overrides the vertical resolution is sized by the levels it asked for.
export function frameBytes(details, cellsPerPanel, levels = details.levels) {
  const cells = 6 * cellsPerPanel ** 2;
  return details.frameSchemaVersion === 2 ? 8 * cells * (1 + 7 * levels) : 32 * cells;
}

export function validateRunRequest(value, presets) {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw protocolError("request must be an object");
  if (value.protocolVersion !== protocolVersion || typeof value.presetId !== "string" || !presets.has(value.presetId)) throw protocolError("unsupported protocol or preset");
  const grid = value.grid; const run = value.run; const initialCondition = value.initialCondition;
  if (!grid || !Number.isInteger(grid.cellsPerPanel) || grid.cellsPerPanel < 1 || grid.cellsPerPanel > 96 ||
      !run || !initialCondition || !finite(run.endTimeSeconds) || run.endTimeSeconds <= 0 || run.endTimeSeconds > 31536000 ||
      !finite(run.maximumTimeStepSeconds) || run.maximumTimeStepSeconds <= 0 || run.maximumTimeStepSeconds > 86400 ||
      !Number.isInteger(run.frameIntervalSteps) || run.frameIntervalSteps < 1 || run.frameIntervalSteps > 1000000 ||
      !Array.isArray(initialCondition.edits) || initialCondition.edits.length > 64) throw protocolError("run parameters are outside the supported range");
  if (estimatedFrameCount(run) > maxPublishedFrames) throw protocolError(`the run would publish more than ${maxPublishedFrames} frames`);
  const details = presetDescriptor(presets.get(value.presetId), value.presetId);
  if (grid.cellsPerPanel > details.maximumCellsPerPanel) throw protocolError("grid exceeds preset limit");
  if (initialCondition.edits.length !== 0 && !details.supportedEdits.includes("gaussian_depth")) {
    throw protocolError(`${details.modelKind} presets do not support initial edits`);
  }
  // A vertical resolution override only exists for a preset that declares levels; a
  // shallow-water preset has no column to re-resolve.
  if (grid.levels !== undefined) {
    if (details.levels === null) throw protocolError(`${details.modelKind} presets have no vertical levels to override`);
    if (!Number.isInteger(grid.levels) || grid.levels < 1 || grid.levels > details.maximumLevels) {
      throw protocolError(`levels must be an integer in [1, ${details.maximumLevels}]`);
    }
  }
  const levels = grid.levels ?? details.levels;
  if (estimatedFrameCount(run) * frameBytes(details, grid.cellsPerPanel, levels) > maxPublishedBytes) {
    throw protocolError("the run would exceed the published byte budget");
  }
  for (const edit of initialCondition.edits) {
    if (!edit || edit.kind !== "gaussian_depth" || typeof edit.id !== "string" || !Array.isArray(edit.centerUnit) || edit.centerUnit.length !== 3 ||
        !edit.centerUnit.every(finite) || Math.abs(Math.hypot(...edit.centerUnit) - 1) > 1e-12 || !finite(edit.amplitudeMeters) ||
        Math.abs(edit.amplitudeMeters) > 1e6 || !finite(edit.sigmaRadians) || edit.sigmaRadians <= 1e-6 || edit.sigmaRadians > Math.PI ||
        (edit.massPolicy !== "preserve_global" && edit.massPolicy !== "allow_change")) throw protocolError("invalid initial edit");
  }
  return value;
}

export function controlRequestText(request, runId) {
  const lines = [
    "control.format_version = 1",
    `control.run_id = ${runId}`,
    `control.cells_per_panel = ${request.grid.cellsPerPanel}`,
    ...(request.grid.levels === undefined ? [] : [`control.levels = ${request.grid.levels}`]),
    `control.end_time_s = ${request.run.endTimeSeconds}`,
    `control.maximum_time_step_s = ${request.run.maximumTimeStepSeconds}`,
    `control.frame_interval_steps = ${request.run.frameIntervalSteps}`,
    "control.frame_directory = frames",
    `initial_edits.count = ${request.initialCondition.edits.length}`,
  ];
  request.initialCondition.edits.forEach((edit, index) => {
    lines.push(`initial_edits.${index}.kind = ${edit.kind}`,
      `initial_edits.${index}.center_x = ${edit.centerUnit[0]}`,
      `initial_edits.${index}.center_y = ${edit.centerUnit[1]}`,
      `initial_edits.${index}.center_z = ${edit.centerUnit[2]}`,
      `initial_edits.${index}.amplitude_m = ${edit.amplitudeMeters}`,
      `initial_edits.${index}.sigma_rad = ${edit.sigmaRadians}`,
      `initial_edits.${index}.mass_policy = ${edit.massPolicy}`);
  });
  return `${lines.join("\n")}\n`;
}

function allowedHost(host) {
  if (!host) return false;
  try { const hostname = new URL(`http://${host}`).hostname; return hostname === "127.0.0.1" || hostname === "localhost" || hostname === "[::1]" || hostname === "::1"; }
  catch { return false; }
}

function allowedOrigin(origin) {
  if (!origin) return true;
  try { const hostname = new URL(origin).hostname; return hostname === "127.0.0.1" || hostname === "localhost" || hostname === "[::1]" || hostname === "::1"; }
  catch { return false; }
}

function json(response, status, value) {
  response.statusCode = status; response.setHeader("content-type", "application/json; charset=utf-8");
  response.end(JSON.stringify(value));
}

async function body(request) {
  const chunks = []; let size = 0;
  for await (const chunk of request) { size += chunk.length; if (size > maxBodyBytes) throw protocolError("request body is too large"); chunks.push(chunk); }
  try { return JSON.parse(Buffer.concat(chunks).toString("utf8")); } catch { throw protocolError("request body is not valid JSON"); }
}

export function createGatewayServer(options) {
  const binary = resolve(options.binary);
  const presets = new Map(Object.entries(options.presets ?? {}));
  const presetDetails = [...presets.entries()].map(([id, value]) => {
    const { path, ...descriptor } = presetDescriptor(value, id);
    return descriptor;
  });
  const runRoot = resolve(options.runRoot);
  const token = options.sessionToken ?? randomBytes(32).toString("hex");
  const spawn = options.spawn ?? defaultSpawn;
  const runs = new Map();
  const retainedRunLimit = options.maxRetainedRuns ?? maxRetainedRuns;
  if (!Number.isInteger(retainedRunLimit) || retainedRunLimit < 0) {
    throw new Error("maxRetainedRuns must be a non-negative integer");
  }
  let activeRunId = null;
  let accessOrder = 0;
  let cleanupQueue = Promise.resolve();
  const touch = (run) => { run.lastAccessOrder = accessOrder++; };
  const scheduleRetention = () => {
    cleanupQueue = cleanupQueue.then(async () => {
      const terminal = [...runs.values()].filter((run) => run.terminal)
        .sort((left, right) => left.lastAccessOrder - right.lastAccessOrder);
      while (terminal.length > retainedRunLimit) {
        const victim = terminal.shift();
        runs.delete(victim.runId);
        victim.events.length = 0; victim.stderr = ""; victim.pending = "";
        victim.waiters.clear(); victim.child = null;
        await rm(victim.directory, { recursive: true, force: true });
      }
    });
  };
  const server = createServer(async (request, response) => {
    try {
      if (!allowedHost(request.headers.host) || !allowedOrigin(request.headers.origin)) return json(response, 403, { error: "forbidden" });
      if (request.headers.origin) {
        response.setHeader("access-control-allow-origin", request.headers.origin);
        response.setHeader("vary", "Origin");
        response.setHeader("access-control-allow-headers", "authorization, content-type, x-after-sequence");
        response.setHeader("access-control-allow-methods", "GET, POST, OPTIONS");
      }
      if (request.method === "OPTIONS") { response.statusCode = 204; response.end(); return; }
      if (request.headers.authorization !== `Bearer ${token}`) return json(response, 401, { error: "unauthorized" });
      const url = new URL(request.url ?? "/", "http://127.0.0.1");
      if (request.method === "GET" && url.pathname === "/api/v1/capabilities") return json(response, 200, { protocolVersion, presets: [...presets.keys()], supportedEdits: ["gaussian_depth"], presetDetails });
      if (request.method === "POST" && url.pathname === "/api/v1/runs") {
        if (activeRunId !== null) return json(response, 409, { error: "run_in_progress" });
        const requestValue = validateRunRequest(await body(request), presets);
        if (activeRunId !== null) return json(response, 409, { error: "run_in_progress" });
        let run;
        const runId = randomBytes(16).toString("hex");
        activeRunId = runId;
        try {
          await mkdir(runRoot, { recursive: true });
          const directory = await mkdtemp(join(runRoot, "run-"));
          await mkdir(join(directory, "frames"));
          const configPath = resolve(presetDescriptor(presets.get(requestValue.presetId), requestValue.presetId).path);
          const controlPath = join(directory, "control.request");
          await writeFile(controlPath, controlRequestText(requestValue, runId), "utf8");
          const child = spawn(binary, ["--config", configPath, "--control-request", controlPath, "--event-stream", "ndjson"], { cwd: directory, shell: false, stdio: ["ignore", "pipe", "pipe"] });
          run = { runId, directory, status: "running", child, request: requestValue, controlRequest: controlRequestText(requestValue, runId), events: [], stderr: "", pending: "", waiters: new Set(), lastSequence: -1, cancellationRequested: false };
          touch(run);
          runs.set(runId, run);
        } catch (error) {
          if (activeRunId === runId) activeRunId = null;
          throw error;
        }
        const { child } = run;
        const publish = (event) => {
          if (!event || event.protocolVersion !== protocolVersion || event.runId !== runId || !Number.isInteger(event.sequence) || event.sequence <= run.lastSequence) throw protocolError("event sequence is invalid");
          run.lastSequence = event.sequence; run.events.push(event);
          if (event.type === "frame.ready") {
            run.frameCount = (run.frameCount ?? 0) + 1;
            run.publishedBytes = (run.publishedBytes ?? 0) + (event.byteLength ?? 0);
            if (run.frameCount > (options.maxPublishedFrames ?? maxPublishedFrames) && !run.terminal && !run.cancellationRequested) {
              run.cancellationRequested = true; run.status = "cancelling"; run.frameBudgetExceeded = true;
              run.stderr = `${run.stderr}gateway stopped the run after the published frame budget was exceeded\n`.slice(-64 * 1024);
              run.child.kill?.("SIGTERM");
              run.forceTimer = setTimeout(() => { if (!run.terminal) run.child.kill?.("SIGKILL"); }, options.graceMilliseconds ?? 2000);
            }
            if (run.publishedBytes > (options.maxPublishedBytes ?? maxPublishedBytes) && !run.terminal && !run.cancellationRequested) { run.cancellationRequested = true; run.status = "cancelling"; run.frameBudgetExceeded = true; run.child.kill?.("SIGTERM"); }
          }
          if (event.type === "run.completed" || event.type === "run.cancelled" || event.type === "run.failed") { run.status = event.type === "run.completed" ? "completed" : event.type === "run.cancelled" ? "cancelled" : "failed"; run.terminal = true; }
          for (const waiter of run.waiters) waiter(event);
        };
        child.stdout?.on("data", (chunk) => {
          run.pending += `${chunk}`;
          const lines = run.pending.split("\n"); run.pending = lines.pop() ?? "";
          for (const line of lines) if (line.trim()) { try { publish(JSON.parse(line)); } catch { run.status = "failed"; run.protocolError = true; child.kill?.("SIGTERM"); } }
        });
        child.on("close", (code, signal) => { if (!run.terminal) run.status = run.cancellationRequested ? "cancelled" : code === 0 ? "completed" : "failed"; run.exitCode = code; run.signal = signal; run.terminal = true; if (activeRunId === runId) activeRunId = null; if (run.forceTimer) clearTimeout(run.forceTimer); for (const waiter of run.waiters) waiter(); run.waiters.clear(); touch(run); scheduleRetention(); });
        child.on("error", (error) => { run.stderr = `${run.stderr}${error.message}`.slice(-64 * 1024); run.status = "failed"; run.terminal = true; if (activeRunId === runId) activeRunId = null; for (const waiter of run.waiters) waiter(); run.waiters.clear(); touch(run); scheduleRetention(); });
        child.stderr?.on("data", (chunk) => { run.stderr = `${run.stderr}${chunk}`.slice(-64 * 1024); });
        return json(response, 202, { protocolVersion, runId, status: "running" });
      }
      const match = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)$/);
      if (request.method === "GET" && match) {
        const run = runs.get(match[1]); if (!run) return json(response, 404, { error: "not_found" }); touch(run);
        return json(response, 200, { protocolVersion, runId: run.runId, status: run.status, exitCode: run.exitCode ?? null });
      }
      const eventsMatch = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)\/events$/);
      if (request.method === "GET" && eventsMatch) {
        const run = runs.get(eventsMatch[1]); if (!run) return json(response, 404, { error: "not_found" }); touch(run);
        const after = Number(request.headers["x-after-sequence"] ?? url.searchParams.get("after") ?? -1);
        if (!Number.isInteger(after)) return json(response, 400, { error: "invalid_sequence" });
        response.statusCode = 200; response.setHeader("content-type", "application/x-ndjson; charset=utf-8");
        if (run.terminal) {
          response.end(run.events.filter((event) => event.sequence > after).map((event) => JSON.stringify(event)).join("\n") + "\n");
          return;
        }
        const send = (event) => { if (event && event.sequence > after) response.write(`${JSON.stringify(event)}\n`); };
        for (const event of run.events) send(event);
        const waiter = (event) => { send(event); if (run.terminal) response.end(); };
        run.waiters.add(waiter); request.on("close", () => run.waiters.delete(waiter));
        return;
      }
      const frameMatch = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)\/frames\/(\d+)$/);
      if (request.method === "GET" && frameMatch) {
        const run = runs.get(frameMatch[1]); if (!run) return json(response, 404, { error: "not_found" }); touch(run);
        const sequence = Number(frameMatch[2]); const event = run.events.find((candidate) => candidate.type === "frame.ready" && candidate.frameSequence === sequence);
        if (!event) return json(response, 404, { error: "frame_not_published" });
        const path = resolve(run.directory, "frames", event.relativePath); const root = resolve(run.directory, "frames");
        if (!path.startsWith(`${root}/`) || event.relativePath.includes("..")) return json(response, 404, { error: "frame_not_published" });
        const stat = await lstat(path); if (stat.isSymbolicLink()) return json(response, 404, { error: "frame_not_published" });
        const realRoot = await realpath(root); if (!(await realpath(path)).startsWith(`${realRoot}/`)) return json(response, 404, { error: "frame_not_published" });
        const bytes = await readFile(path); response.statusCode = 200; response.setHeader("content-type", "application/octet-stream"); response.setHeader("content-length", bytes.byteLength); response.end(bytes); return;
      }
      const cancelMatch = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)\/cancel$/);
      if (request.method === "POST" && cancelMatch) {
        const run = runs.get(cancelMatch[1]); if (!run) return json(response, 404, { error: "not_found" }); touch(run);
        if (!run.terminal) { run.cancellationRequested = true; run.status = "cancelling"; run.child.kill?.("SIGTERM"); run.forceTimer = setTimeout(() => { if (!run.terminal) run.child.kill?.("SIGKILL"); }, options.graceMilliseconds ?? 2000); }
        return json(response, 202, { protocolVersion, runId: run.runId, status: run.status });
      }
      const bundleMatch = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)\/bundle$/);
      if (request.method === "GET" && bundleMatch) {
        const run = runs.get(bundleMatch[1]); if (!run) return json(response, 404, { error: "not_found" }); touch(run);
        return json(response, 200, { protocolVersion, runId: run.runId, status: run.status, request: run.request, controlRequest: run.controlRequest, events: run.events, stderr: run.stderr });
      }
      return json(response, 404, { error: "not_found" });
    } catch (error) { return json(response, error.code === "invalid_request" ? 400 : 500, { error: error.code ?? "gateway_error", message: error.code === "invalid_request" ? error.message : "request failed" }); }
  });
  const shutdown = async () => {
    for (const run of runs.values()) if (!run.terminal) { run.cancellationRequested = true; run.child?.kill?.("SIGTERM"); }
    server.closeAllConnections?.();
    await new Promise((resolveClose) => server.close(resolveClose));
    await cleanupQueue;
  };
  return { server, token, runs, shutdown, flushRetention: async () => cleanupQueue };
}

export async function describeSimulator(options) {
  const spawn = options.spawn ?? defaultSpawn;
  return new Promise((resolveDescription, reject) => {
    const child = spawn(resolve(options.binary), ["--describe-control"], { shell: false, stdio: ["ignore", "pipe", "pipe"] });
    let output = ""; child.stdout?.on("data", (chunk) => { output += chunk; });
    child.on("error", reject); child.on("close", (code) => { if (code !== 0) return reject(new Error("simulator capability probe failed")); try { const value = JSON.parse(output); if (value.protocolVersion !== protocolVersion || !value.supportedEdits?.includes("gaussian_depth")) return reject(new Error("simulator protocol is incompatible")); resolveDescription(value); } catch { reject(new Error("simulator capability probe was not JSON")); } });
  });
}

// A preset descriptor promises the browser a model kind, a frame schema, and an N limit.
// If the simulator on disk cannot honour that promise the gateway refuses to start, rather
// than failing once a reader has already submitted a run.
export function assertPresetsAreSupported(description, presets) {
  for (const [id, value] of Object.entries(presets ?? {})) {
    const details = presetDescriptor(value, id);
    const kinds = description.modelKinds ?? ["shallow_water"];
    const schemas = description.frameSchemaVersions ?? [1];
    if (!kinds.includes(details.modelKind)) throw new Error(`the simulator does not support ${details.modelKind} preset ${id}`);
    if (!schemas.includes(details.frameSchemaVersion)) throw new Error(`the simulator does not publish frame schema ${details.frameSchemaVersion} for preset ${id}`);
    if (details.modelKind !== "dry_hydrostatic") continue;
    const limits = description.dryHydrostatic ?? {};
    if (!Number.isInteger(details.levels) || details.levels < 1 || details.levels > (limits.maxLevels ?? 0)) {
      throw new Error(`preset ${id} declares an unsupported level count`);
    }
    if (details.maximumCellsPerPanel > (limits.maxCellsPerPanel ?? 0)) {
      throw new Error(`preset ${id} exceeds the interactive dry hydrostatic grid limit`);
    }
    if (details.supportedEdits.length !== 0) throw new Error(`preset ${id} cannot support initial edits`);
  }
}

export async function startGateway(options) {
  assertPresetsAreSupported(await describeSimulator(options), options.presets);
  const gateway = createGatewayServer(options);
  const host = options.host ?? "127.0.0.1";
  if (host !== "127.0.0.1" && host !== "::1") throw new Error("gateway host must be loopback");
  await new Promise((resolveListen, reject) => { gateway.server.once("error", reject); gateway.server.listen(options.port ?? 0, host, resolveListen); });
  return gateway;
}

// The descriptor's level count must equal the preset's own vertical.levels or the byte
// budget is computed against a column the solver will not produce, so it is read from the
// configuration rather than repeated in an environment variable.
export function readPresetLevels(path) {
  const match = /^\s*vertical\.levels\s*=\s*(\d+)\s*$/m.exec(readFileSync(path, "utf8"));
  if (!match) throw new Error(`${path} does not declare vertical.levels`);
  return Number(match[1]);
}

export function dryPresetsFromEnvironment(environment) {
  const paths = (environment.MPS_DRY_PRESETS ?? environment.MPS_DRY_PRESET ?? "")
    .split(",").map((value) => value.trim()).filter(Boolean);
  const presets = {};
  for (const path of paths) {
    const id = basename(path).replace(/\.cfg$/, "");
    presets[id] = { path, modelKind: "dry_hydrostatic", frameSchemaVersion: 2,
      levels: readPresetLevels(path), supportedEdits: [],
      maximumCellsPerPanel: Number(environment.MPS_DRY_PRESET_MAX_N ?? 24) };
  }
  return presets;
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const binary = process.env.MPS_SIMULATOR_BINARY; const preset = process.env.MPS_REST_PRESET;
  if (!binary || !preset) throw new Error("MPS_SIMULATOR_BINARY and MPS_REST_PRESET are required");
  // Dry presets are optional so an existing shallow-water session keeps working.
  const presets = { rest: preset, ...dryPresetsFromEnvironment(process.env) };
  const gateway = await startGateway({ binary, presets, runRoot: process.env.MPS_RUN_ROOT ?? ".runs", sessionToken: process.env.MPS_SESSION_TOKEN, port: Number(process.env.MPS_GATEWAY_PORT ?? 0) });
  const address = gateway.server.address();
  console.log(JSON.stringify({ host: "127.0.0.1", port: typeof address === "object" ? address.port : address, token: gateway.token }));
}
