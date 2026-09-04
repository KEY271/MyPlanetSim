import { existsSync } from "node:fs";
import { mkdtemp } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawn } from "node:child_process";
import { test } from "node:test";
import assert from "node:assert/strict";
import { startGateway } from "../../gateway/src/main.js";

// Everything else in this repository exercises the viewer through node. Node accepts API
// shapes a browser rejects, and it never renders, so two defects reached the viewer
// unnoticed: fetch called with the wrong receiver, and the tokenless dev-server URL
// silently selecting the offline demo. This gate opens the real page in a real browser.
const enabled = process.env.MPS_ENABLE_BROWSER === "1";

async function loadPlaywright() {
  try {
    return await import("playwright");
  } catch {
    return null;
  }
}

async function waitForHttp(url, headers = {}) {
  for (let attempt = 0; attempt < 120; attempt += 1) {
    try {
      const response = await fetch(url, { headers });
      if (response.ok) return true;
    } catch {
      // not listening yet
    }
    await new Promise((resolveSleep) => setTimeout(resolveSleep, 250));
  }
  return false;
}

test("the real browser reaches the gateway and offers the dry hydrostatic preset", { skip: !enabled }, async () => {
  const playwright = await loadPlaywright();
  if (!playwright) return;
  const root = resolve(process.cwd(), "../../..");
  const binary = resolve(root, "build/dev/my_planet_sim");
  const dryPreset = resolve(root, "configs/phase5_visualizer_rest_n4.cfg");
  if (!existsSync(binary) || !existsSync(dryPreset)) return;

  const runRoot = await mkdtemp(join(tmpdir(), "myplanetsim-browser-"));
  const token = "browser-gate-token";
  const gateway = await startGateway({
    binary, runRoot, sessionToken: token, port: 0,
    presets: {
      rest: resolve(root, "configs/phase3_rest_n4.cfg"),
      dry_hydrostatic_rest: { path: dryPreset, modelKind: "dry_hydrostatic", frameSchemaVersion: 2,
        levels: 8, supportedEdits: [], maximumCellsPerPanel: 24 },
    },
  });
  const gatewayPort = gateway.server.address().port;
  const uiPort = 5290;
  const vite = spawn("npm", ["run", "dev", "--", "--host", "127.0.0.1", "--port", String(uiPort), "--strictPort"],
    { cwd: resolve(process.cwd()), shell: false, stdio: ["ignore", "pipe", "pipe"] });
  let browser;
  try {
    assert.ok(await waitForHttp(`http://127.0.0.1:${uiPort}/`), "the dev server did not start");
    try {
      browser = await playwright.chromium.launch();
    } catch (error) {
      // The browser binary is a separate download; report the command rather than failing
      // a checkout that has not run it.
      console.log(`skipping: chromium is not installed (${error.message.split("\n")[0]}); run "npx playwright install chromium"`);
      return;
    }
    const failures = [];
    // The session token is read once at module load, so each URL needs its own document.
    const open = async (path) => {
      const page = await browser.newPage();
      page.on("pageerror", (error) => failures.push(String(error.message)));
      await page.goto(`http://127.0.0.1:${uiPort}${path}`, { waitUntil: "networkidle" });
      return page;
    };

    // The tokenless URL must announce the offline fallback rather than looking like the
    // dry hydrostatic model is missing.
    const demo = await open("/");
    assert.match(await demo.textContent(".engine-status"), /deterministic demo/);
    assert.match(await demo.textContent(".demo-banner"), /Offline demo/);
    assert.deepEqual(await demo.$$eval("select", (nodes) => [...nodes[0].options].map((o) => o.value)), ["rest"]);

    // The tokenised URL must reach the gateway and expose the dry preset and its levels.
    const live = await open(`/#token=${token}&gateway=${encodeURIComponent(`http://127.0.0.1:${gatewayPort}`)}`);
    await live.waitForFunction(() => document.querySelectorAll("select")[0]?.options.length > 1, null, { timeout: 15000 });
    assert.equal(await live.$(".demo-banner"), null);
    assert.match(await live.textContent(".engine-status"), /native C\+\+/);
    assert.equal(await live.$(".preset-error"), null);
    const presetOptions = await live.$$eval("select", (nodes) => [...nodes[0].options].map((o) => o.value));
    assert.deepEqual(presetOptions, ["rest", "dry_hydrostatic_rest"]);

    // Selecting the dry preset must reveal the vertical controls before any run starts.
    await live.selectOption("select", "dry_hydrostatic_rest");
    assert.match(await live.textContent(".engine-status"), /dry_hydrostatic/);
    assert.match(await live.textContent(".preset-hint"), /8 model levels/);
    const ranges = await live.$$eval("input[type=range]", (nodes) => nodes.map((n) => n.getAttribute("aria-label")));
    assert.ok(ranges.includes("Model level"), `expected a model level slider, saw ${JSON.stringify(ranges)}`);
    const kControl = await live.$("input[aria-label='K vertical resolution']");
    assert.ok(kControl, "expected a K control");

    // Run the preset at a resolution the reader chose and confirm the vertical view is
    // actually driven by the published frames, not just announced.
    await kControl.fill("12");
    await live.fill("input[type=number][max='31536000']", "60");
    await live.fill("input[type=number][max='86400']", "10");
    await live.fill("input[type=number][max='1000000']", "1");
    await live.click("button:text('Run')");
    await live.waitForFunction(() => document.querySelector(".run-status")?.textContent?.includes("completed"), null, { timeout: 60000 });
    assert.match(await live.textContent(".engine-status"), /frame schema: 2/);
    assert.equal(await live.$(".preset-hint"), null);
    const levelSlider = await live.$("input[aria-label='Model level']");
    assert.equal(await levelSlider.getAttribute("max"), "11", "the level slider must span the K that ran");
    assert.equal(await levelSlider.isDisabled(), false);
    assert.ok(await live.$("section.column-profile, article.column-profile"), "expected the column profile");
    assert.ok((await live.textContent(".inspector")).includes("level 0/11"));

    // Moving the slider must change the inspected level and the reported pressure.
    const pressureAtTop = await live.textContent(".inspector");
    await levelSlider.fill("11");
    await live.waitForFunction(() => document.querySelector(".inspector")?.textContent?.includes("level 11/11"), null, { timeout: 10000 });
    assert.notEqual(await live.textContent(".inspector"), pressureAtTop);
    assert.deepEqual(failures, []);
  } finally {
    await browser?.close();
    vite.kill("SIGTERM");
    await gateway.shutdown();
  }
});
