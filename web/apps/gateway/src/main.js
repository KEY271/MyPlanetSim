import { randomBytes } from "node:crypto";
import { createServer } from "node:http";
import { mkdir, mkdtemp, writeFile } from "node:fs/promises";
import { spawn as defaultSpawn } from "node:child_process";
import { join, resolve } from "node:path";

const protocolVersion = 1;
const maxBodyBytes = 64 * 1024;

function protocolError(message) { const error = new Error(message); error.code = "invalid_request"; return error; }

function finite(value) { return typeof value === "number" && Number.isFinite(value); }

export function validateRunRequest(value, presets) {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw protocolError("request must be an object");
  if (value.protocolVersion !== protocolVersion || typeof value.presetId !== "string" || !presets.has(value.presetId)) throw protocolError("unsupported protocol or preset");
  const run = value.run; const initialCondition = value.initialCondition;
  if (!run || !initialCondition || !finite(run.endTimeSeconds) || run.endTimeSeconds <= 0 || run.endTimeSeconds > 31536000 ||
      !finite(run.maximumTimeStepSeconds) || run.maximumTimeStepSeconds <= 0 || run.maximumTimeStepSeconds > 86400 ||
      !Number.isInteger(run.frameIntervalSteps) || run.frameIntervalSteps < 1 || run.frameIntervalSteps > 1000000 ||
      !Array.isArray(initialCondition.edits) || initialCondition.edits.length > 64) throw protocolError("run parameters are outside the supported range");
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
      if (request.headers.authorization !== `Bearer ${token}`) return json(response, 401, { error: "unauthorized" });
      const url = new URL(request.url ?? "/", "http://127.0.0.1");
      if (request.method === "GET" && url.pathname === "/api/v1/capabilities") return json(response, 200, { protocolVersion, presets: [...presets.keys()], supportedEdits: ["gaussian_depth"] });
      if (request.method === "POST" && url.pathname === "/api/v1/runs") {
        if (activeRun) return json(response, 409, { error: "run_in_progress" });
        const requestValue = validateRunRequest(await body(request), presets);
        const runId = randomBytes(16).toString("hex");
        await mkdir(runRoot, { recursive: true });
        const directory = await mkdtemp(join(runRoot, "run-"));
        await mkdir(join(directory, "frames"));
        const configPath = resolve(presets.get(requestValue.presetId));
        const controlPath = join(directory, "control.request");
        await writeFile(controlPath, controlRequestText(requestValue, runId), "utf8");
        const child = spawn(binary, ["--config", configPath, "--control-request", controlPath, "--event-stream", "ndjson"], { cwd: directory, shell: false, stdio: ["ignore", "pipe", "pipe"] });
        const run = { runId, directory, status: "running", child, events: [], stderr: "" };
        runs.set(runId, run); activeRun = true;
        child.on("close", (code, signal) => { run.status = code === 0 ? "completed" : "failed"; run.exitCode = code; run.signal = signal; activeRun = false; });
        child.stderr?.on("data", (chunk) => { run.stderr = `${run.stderr}${chunk}`.slice(-64 * 1024); });
        return json(response, 202, { protocolVersion, runId, status: "running" });
      }
      const match = url.pathname.match(/^\/api\/v1\/runs\/([^/]+)$/);
      if (request.method === "GET" && match) {
        const run = runs.get(match[1]); if (!run) return json(response, 404, { error: "not_found" });
        return json(response, 200, { protocolVersion, runId: run.runId, status: run.status, exitCode: run.exitCode ?? null });
      }
      return json(response, 404, { error: "not_found" });
    } catch (error) { return json(response, error.code === "invalid_request" ? 400 : 500, { error: error.code ?? "gateway_error", message: error.code === "invalid_request" ? error.message : "request failed" }); }
  });
  return { server, token, runs };
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
  const gateway = await startGateway({ binary, presets: { rest: preset }, runRoot: process.env.MPS_RUN_ROOT ?? ".runs", port: Number(process.env.MPS_GATEWAY_PORT ?? 0) });
  const address = gateway.server.address();
  console.log(JSON.stringify({ host: "127.0.0.1", port: typeof address === "object" ? address.port : address, token: gateway.token }));
}
