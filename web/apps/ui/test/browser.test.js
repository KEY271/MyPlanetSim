import { existsSync } from "node:fs";
import { mkdtemp } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawn } from "node:child_process";
import { test } from "node:test";
import assert from "node:assert/strict";
import { dryPresetsFromEnvironment, startGateway } from "../../gateway/src/main.js";

// Everything else in this repository exercises the viewer through node. Node accepts API
// shapes a browser rejects and never renders, so several defects reached the viewer
// unnoticed: fetch called with the wrong receiver, the tokenless dev-server URL silently
// selecting the offline demo, and a "visualizer" preset that never moves. This gate opens
// the real page in a real browser against the real gateway.
const enabled = process.env.MPS_ENABLE_BROWSER === "1";
const root = resolve(process.cwd(), "../../..");

async function loadPlaywright() {
  try {
    return await import("playwright");
  } catch {
    return null;
  }
}

async function waitForHttp(url) {
  for (let attempt = 0; attempt < 160; attempt += 1) {
    try {
      if ((await fetch(url)).ok) return true;
    } catch {
      // not listening yet
    }
    await new Promise((resolveSleep) => setTimeout(resolveSleep, 250));
  }
  return false;
}

function startDevServer(port, environment) {
  return spawn("npm", ["run", "dev", "--", "--host", "127.0.0.1", "--port", String(port), "--strictPort"],
    { cwd: process.cwd(), shell: false, stdio: ["ignore", "pipe", "pipe"], env: { ...process.env, ...environment } });
}

test("the viewer reaches the gateway, offers moving presets, and renders their evolution", { skip: !enabled }, async () => {
  const playwright = await loadPlaywright();
  if (!playwright) return;
  const binary = resolve(root, "build/dev/my_planet_sim");
  const moving = resolve(root, "configs/phase5_umjs14_baroclinic.cfg");
  if (!existsSync(binary) || !existsSync(moving)) return;

  const runRoot = await mkdtemp(join(tmpdir(), "myplanetsim-browser-"));
  const token = "browser-gate-token";
  const gateway = await startGateway({
    binary, runRoot, sessionToken: token, port: 0,
    presets: {
      rest: resolve(root, "configs/phase3_rest_n4.cfg"),
      ...dryPresetsFromEnvironment({ MPS_DRY_PRESETS: `${moving},${resolve(root, "configs/phase5_visualizer_rest_n4.cfg")}` }),
    },
  });
  const gatewayPort = gateway.server.address().port;
  // One dev server carries the session `just dev` exports; the other has none, which is the
  // offline case a reader gets from `just ui`.
  const livePort = 5290;
  const demoPort = 5291;
  const live = startDevServer(livePort, { MPS_SESSION_TOKEN: token, MPS_GATEWAY_PORT: String(gatewayPort) });
  const demo = startDevServer(demoPort, { MPS_SESSION_TOKEN: "", MPS_GATEWAY_PORT: "" });
  let browser;
  try {
    assert.ok(await waitForHttp(`http://127.0.0.1:${livePort}/`), "the live dev server did not start");
    assert.ok(await waitForHttp(`http://127.0.0.1:${demoPort}/`), "the demo dev server did not start");
    try {
      browser = await playwright.chromium.launch();
    } catch (error) {
      // The browser binary is a separate download; report the command rather than failing
      // a checkout that has not run it.
      console.log(`skipping: chromium is not installed (${error.message.split("\n")[0]}); run "npx playwright install chromium"`);
      return;
    }
    const failures = [];
    const open = async (url) => {
      const page = await browser.newPage();
      page.on("pageerror", (error) => failures.push(String(error.message)));
      await page.goto(url, { waitUntil: "networkidle" });
      return page;
    };

    // Without an injected session the plain URL is the offline demo, and says so.
    const offline = await open(`http://127.0.0.1:${demoPort}/`);
    assert.match(await offline.textContent(".engine-status"), /deterministic demo/);
    assert.match(await offline.textContent(".demo-banner"), /Offline demo/);
    assert.deepEqual(await offline.$$eval("select", (nodes) => [...nodes[0].options].map((o) => o.value)), ["rest"]);

    // The URL `just dev` prints carries no token, so this is the exact address a reader
    // opens. It must reach the gateway anyway.
    const page = await open(`http://127.0.0.1:${livePort}/`);
    await page.waitForFunction(() => document.querySelectorAll("select")[0]?.options.length > 1, null, { timeout: 15000 });
    assert.equal(await page.$(".demo-banner"), null);
    assert.equal(await page.$(".preset-error"), null);
    assert.match(await page.textContent(".engine-status"), /native C\+\+/);
    const presetOptions = await page.$$eval("select", (nodes) => [...nodes[0].options].map((o) => o.value));
    assert.ok(presetOptions.includes("phase5_umjs14_baroclinic"), `expected a moving dry preset, saw ${JSON.stringify(presetOptions)}`);

    // Selecting a dry preset must reveal the vertical controls before any run starts.
    await page.selectOption("select", "phase5_umjs14_baroclinic");
    assert.match(await page.textContent(".engine-status"), /dry_hydrostatic/);
    assert.match(await page.textContent(".preset-hint"), /8 model levels/);
    const kControl = await page.$("input[aria-label='K vertical resolution']");
    assert.ok(kControl, "expected a K control");

    // Run it, then require the displayed field to actually change over the run. A preset
    // that renders a still image is not a visualizer.
    await kControl.fill("12");
    await page.fill("input[type=number][max='31536000']", "1800");
    await page.fill("input[type=number][max='86400']", "60");
    await page.fill("input[type=number][max='1000000']", "1");
    await page.click("button:text('Run')");
    await page.waitForFunction(() => document.querySelector(".run-status")?.textContent?.includes("completed"), null, { timeout: 120000 });
    assert.match(await page.textContent(".engine-status"), /frame schema: 2/);

    const levelSlider = await page.$("input[aria-label='Model level']");
    assert.equal(await levelSlider.getAttribute("max"), "11", "the level slider must span the K that ran");
    assert.equal(await levelSlider.isDisabled(), false);
    assert.ok(await page.$(".column-profile"), "expected the column profile");

    const frameSlider = await page.$("input[aria-label='Frame']");
    const frameCount = Number(await frameSlider.getAttribute("max"));
    assert.ok(frameCount >= 2, `expected several frames, saw ${frameCount + 1}`);
    const inspectorAt = async (frame) => {
      await frameSlider.fill(String(frame));
      await page.waitForTimeout(50);
      return page.textContent(".inspector");
    };
    const atStart = await inspectorAt(0);
    const atEnd = await inspectorAt(frameCount);
    assert.notEqual(atStart, atEnd, "the inspected value must change over the run");

    // The level selector must select a genuinely different column entry.
    await levelSlider.fill("11");
    await page.waitForFunction(() => document.querySelector(".inspector")?.textContent?.includes("level 11/11"), null, { timeout: 10000 });
    assert.notEqual(await page.textContent(".inspector"), atEnd);
    assert.deepEqual(failures, []);
  } finally {
    await browser?.close();
    live.kill("SIGTERM");
    demo.kill("SIGTERM");
    await gateway.shutdown();
  }
});
