import { randomBytes } from "node:crypto";
import { createServer } from "node:http";
import { lstat, mkdir, mkdtemp, readFile, realpath, writeFile } from "node:fs/promises";
import { spawn as defaultSpawn } from "node:child_process";
import { join, resolve } from "node:path";

const protocolVersion = 1;
const maxBodyBytes = 64 * 1024;
// Mirrors @myplanetsim/protocol maxPublishedFrames. The native solver may shrink its time step
// below the requested maximum, so the budget is enforced again while the run streams events.
export const maxPublishedFrames = 512;

export function estimatedFrameCount(run) {
  return Math.floor(Math.ceil(run.endTimeSeconds / run.maximumTimeStepSeconds) / run.frameIntervalSteps) + 2;
}

function protocolError(message) { const error = new Error(message); error.code = "invalid_request"; return error; }

function finite(value) { return typeof value === "number" && Number.isFinite(value); }

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
  const runRoot = resolve(options.runRoot);
  const token = options.sessionToken ?? randomBytes(32).toString("hex");
  const spawn = options.spawn ?? defaultSpawn;
  const runs = new Map();
  let activeRun = false;
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
      if (request.method === "GET" && url.pathname === "/api/v1/capabilities") return json(response, 200, { protocolVersion, presets: [...presets.keys()], supportedEdits: ["gaussian_depth"] });
      if (request.method === "POST" && url.pathname === "/api/v1/runs") {
        if (activeRun) return json(response, 409, { error: "run_in_progress" });
        const requestValue = validateRunRequest(await body(request), presets);
        if (activeRun) return json(response, 409, { error: "run_in_progress" });
        activeRun = true;
        let run;
        try {
          const runId = randomBytes(16).toString("hex");
          await mkdir(runRoot, { recursive: true });
          const directory = await mkdtemp(join(runRoot, "run-"));
          await mkdir(join(directory, "frames"));
          const configPath = resolve(presets.get(requestValue.presetId));
          const controlPath = join(directory, "control.request");
          await writeFile(controlPath, controlRequestText(requestValue, runId), "utf8");
          const child = spawn(binary, ["--config", configPath, "--control-request", controlPath, "--event-stream", "ndjson"], { cwd: directory, shell: false, stdio: ["ignore", "pipe", "pipe"] });
          run = { runId, directory, status: "running", child, request: requestValue, controlRequest: controlRequestText(requestValue, runId), events: [], stderr: "", pending: "", waiters: new Set(), lastSequence: -1, cancellationRequested: false };
          runs.set(runId, run);
        } catch (error) {
          activeRun = false;
          throw error;
        }
        const { child, runId } = run;
        const publish = (event) => {
          if (!event || event.protocolVersion !== protocolVersion || event.runId !== runId || !Number.isInteger(event.sequence) || event.sequence <= run.lastSequence) throw protocolError("event sequence is invalid");
          run.lastSequence = event.sequence; run.events.push(event);
          if (event.type === "frame.ready") {
            run.frameCount = (run.frameCount ?? 0) + 1;
            if (run.frameCount > (options.maxPublishedFrames ?? maxPublishedFrames) && !run.terminal && !run.cancellationRequested) {
              run.cancellationRequested = true; run.status = "cancelling"; run.frameBudgetExceeded = true;
              run.stderr = `${run.stderr}gateway stopped the run after the published frame budget was exceeded\n`.slice(-64 * 1024);
              run.child.kill?.("SIGTERM");
              run.forceTimer = setTimeout(() => { if (!run.terminal) run.child.kill?.("SIGKILL"); }, options.graceMilliseconds ?? 2000);
            }
          }
          if (event.type === "run.completed" || event.type === "run.cancelled" || event.type === "run.failed") { run.status = event.type === "run.completed" ? "completed" : event.type === "run.cancelled" ? "cancelled" : "failed"; run.terminal = true; }
          for (const waiter of run.waiters) waiter(event);
        };
        child.stdout?.on("data", (chunk) => {
          run.pending += `${chunk}`;
          const lines = run.pending.split("\n"); run.pending = lines.pop() ?? "";
          for (const line of lines) if (line.trim()) { try { publish(JSON.parse(line)); } catch { run.status = "failed"; run.protocolError = true; child.kill?.("SIGTERM"); } }
        });
        child.on("close", (code, signal) => { if (!run.terminal) run.status = run.cancellationRequested ? "cancelled" : code === 0 ? "completed" : "failed"; run.exitCode = code; run.signal = signal; run.terminal = true; activeRun = false; if (run.forceTimer) clearTimeout(run.forceTimer); for (const waiter of run.waiters) waiter(); run.waiters.clear(); });
        child.on("error", (error) => { run.stderr = `${run.stderr}${error.message}`.slice(-64 * 1024); run.status = "failed"; run.terminal = true; activeRun = false; for (const waiter of run.waiters) waiter(); run.waiters.clear(); });
        child.stderr?.on("data", (chunk) => { run.stderr = `${run.stderr}${chunk}`.slice(-64 * 1024); });
        return json(response, 202, { protocolVersion, runId, status: "running" });
      }
      const match = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)$/);
      if (request.method === "GET" && match) {
        const run = runs.get(match[1]); if (!run) return json(response, 404, { error: "not_found" });
        return json(response, 200, { protocolVersion, runId: run.runId, status: run.status, exitCode: run.exitCode ?? null });
      }
      const eventsMatch = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)\/events$/);
      if (request.method === "GET" && eventsMatch) {
        const run = runs.get(eventsMatch[1]); if (!run) return json(response, 404, { error: "not_found" });
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
        const run = runs.get(frameMatch[1]); if (!run) return json(response, 404, { error: "not_found" });
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
        const run = runs.get(cancelMatch[1]); if (!run) return json(response, 404, { error: "not_found" });
        if (!run.terminal) { run.cancellationRequested = true; run.status = "cancelling"; run.child.kill?.("SIGTERM"); run.forceTimer = setTimeout(() => { if (!run.terminal) run.child.kill?.("SIGKILL"); }, options.graceMilliseconds ?? 2000); }
        return json(response, 202, { protocolVersion, runId: run.runId, status: run.status });
      }
      const bundleMatch = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)\/bundle$/);
      if (request.method === "GET" && bundleMatch) {
        const run = runs.get(bundleMatch[1]); if (!run) return json(response, 404, { error: "not_found" });
        return json(response, 200, { protocolVersion, runId: run.runId, status: run.status, request: run.request, controlRequest: run.controlRequest, events: run.events, stderr: run.stderr });
      }
      return json(response, 404, { error: "not_found" });
    } catch (error) { return json(response, error.code === "invalid_request" ? 400 : 500, { error: error.code ?? "gateway_error", message: error.code === "invalid_request" ? error.message : "request failed" }); }
  });
  const shutdown = async () => {
    for (const run of runs.values()) if (!run.terminal) { run.cancellationRequested = true; run.child.kill?.("SIGTERM"); }
    server.closeAllConnections?.();
    await new Promise((resolveClose) => server.close(resolveClose));
  };
  return { server, token, runs, shutdown };
}

export async function describeSimulator(options) {
  const spawn = options.spawn ?? defaultSpawn;
  return new Promise((resolveDescription, reject) => {
    const child = spawn(resolve(options.binary), ["--describe-control"], { shell: false, stdio: ["ignore", "pipe", "pipe"] });
    let output = ""; child.stdout?.on("data", (chunk) => { output += chunk; });
    child.on("error", reject); child.on("close", (code) => { if (code !== 0) return reject(new Error("simulator capability probe failed")); try { const value = JSON.parse(output); if (value.protocolVersion !== protocolVersion || !value.supportedEdits?.includes("gaussian_depth")) return reject(new Error("simulator protocol is incompatible")); resolveDescription(value); } catch { reject(new Error("simulator capability probe was not JSON")); } });
  });
}

export async function startGateway(options) {
  await describeSimulator(options);
  const gateway = createGatewayServer(options);
  const host = options.host ?? "127.0.0.1";
  if (host !== "127.0.0.1" && host !== "::1") throw new Error("gateway host must be loopback");
  await new Promise((resolveListen, reject) => { gateway.server.once("error", reject); gateway.server.listen(options.port ?? 0, host, resolveListen); });
  return gateway;
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const binary = process.env.MPS_SIMULATOR_BINARY; const preset = process.env.MPS_REST_PRESET;
  if (!binary || !preset) throw new Error("MPS_SIMULATOR_BINARY and MPS_REST_PRESET are required");
  const gateway = await startGateway({ binary, presets: { rest: preset }, runRoot: process.env.MPS_RUN_ROOT ?? ".runs", sessionToken: process.env.MPS_SESSION_TOKEN, port: Number(process.env.MPS_GATEWAY_PORT ?? 0) });
  const address = gateway.server.address();
  console.log(JSON.stringify({ host: "127.0.0.1", port: typeof address === "object" ? address.port : address, token: gateway.token }));
}
