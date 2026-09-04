import { StrictMode, useCallback, useEffect, useMemo, useRef, useState } from "react";
import { createRoot } from "react-dom/client";
import { addShallowWaterDerivedFields, diagnosticEvents, parseShallowWaterCsv, UnitVector } from "@myplanetsim/protocol";
import { DeterministicMockSimulationClient, HttpSimulationClient } from "@myplanetsim/protocol/client";
import { Globe3D } from "./Globe3D";
import { Map2D } from "./Map2D";
import { appendEdit, clearEdits, createGaussianEdit, DraftRun, toRunRequest, undoEdit } from "./editor";
import { SimulationController } from "./controller";
import "./style.css";

function demoDataset(cellsPerPanel: number) {
  const rows = ["panel,i,j,depth_m,momentum_x,momentum_y,momentum_z"];
  for (const panel of ["PX", "PY", "NX", "NY", "PZ", "NZ"]) for (let j = 0; j < cellsPerPanel; j += 1) for (let i = 0; i < cellsPerPanel; i += 1) rows.push(`${panel},${i},${j},3000,0,0,0`);
  return addShallowWaterDerivedFields(parseShallowWaterCsv(rows.join("\n"), cellsPerPanel));
}

const startupParameters = new URLSearchParams(window.location.hash.slice(1));
const startupToken = startupParameters.get("token");
const startupGateway = startupParameters.get("gateway");
if (startupToken) window.history.replaceState(null, "", `${window.location.pathname}${window.location.search}`);

function App() {
  const [draft, setDraft] = useState<DraftRun>({ presetId: "rest", cellsPerPanel: 4, endTimeSeconds: 3600, maximumTimeStepSeconds: 60, frameIntervalSteps: 10, edits: [] });
  const dataset = useMemo(() => demoDataset(draft.cellsPerPanel), [draft.cellsPerPanel]);
  const controller = useMemo(() => {
    const client = startupToken
      ? new HttpSimulationClient(startupGateway ?? window.location.origin, startupToken)
      : new DeterministicMockSimulationClient();
    return new SimulationController(client);
  }, []);
  const [run, setRun] = useState(controller.value);
  const [playing, setPlaying] = useState(false);
  useEffect(() => { const unsubscribe = controller.subscribe(setRun); return () => { unsubscribe(); }; }, [controller]);
  useEffect(() => {
    if (!playing || run.frames.length < 2) return undefined;
    const timer = window.setInterval(() => {
      const currentFrame = controller.value.currentFrame;
      if (currentFrame >= controller.value.frames.length - 1) { setPlaying(false); return; }
      controller.setCurrentFrame(currentFrame + 1);
    }, 250);
    return () => window.clearInterval(timer);
  }, [controller, playing, run.frames.length]);
  const [fieldId, setFieldId] = useState("depth");
  const [gridMode, setGridMode] = useState<"off" | "panel_seams" | "all_cells">("panel_seams");
  const [amplitude, setAmplitude] = useState(100);
  const [sigma, setSigma] = useState(0.2);
  const [massPolicy, setMassPolicy] = useState<"preserve_global" | "allow_change">("preserve_global");
  const [message, setMessage] = useState("Choose a cell in either view to create a pending edit.");
  const editCounter = useRef(1);
  const currentDataset = run.frames[run.currentFrame] ?? dataset;
  const addFromCell = useCallback((cell: number | null, origin?: UnitVector) => {
    if (cell === null || !origin) return;
    // The updater must stay pure: React may replay it, so the identifier and the status
    // message are derived once here instead of inside setDraft.
    const edit = createGaussianEdit(origin, amplitude, sigma, massPolicy, `edit-${editCounter.current++}`);
    setDraft((current) => appendEdit(current, edit));
    setMessage(`${edit.id} added at the clicked location (cell ${cell}); the ring shows sigma.`);
  }, [amplitude, massPolicy, sigma]);
  const requestSummary = () => { try { return JSON.stringify(toRunRequest(draft), null, 2); } catch (error) { return error instanceof Error ? error.message : "invalid request"; } };
  const start = () => { try { void controller.start(toRunRequest(draft)); } catch (error) { setMessage(error instanceof Error ? error.message : "invalid request"); } };
  const togglePlayback = () => { if (playing) { setPlaying(false); return; } if (run.currentFrame >= run.frames.length - 1) controller.setCurrentFrame(0); setPlaying(true); };
  const displayed = run.frames.length > 0 ? currentDataset : dataset;
  const diagnostics = diagnosticEvents(run.events);
  const downloadBundle = async () => { const bundle = await controller.bundle(); if (!bundle) return; const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([JSON.stringify(bundle, null, 2)], { type: "application/json" })); link.download = `${bundle.runId}.json`; link.click(); URL.revokeObjectURL(link.href); };
  return <main className="app-shell">
    <header className="top-panel">
      <div className="title-row"><p className="eyebrow">MyPlanetSim · Phase 3</p><h1>Interactive shallow-water visualizer</h1><p className="message">{message}</p></div>
      <section className="control-panel" aria-label="Initial condition editor">
        <label>Field <select value={fieldId} onChange={(event) => setFieldId(event.target.value)}>{dataset.fields.map((field) => <option key={field.id} value={field.id}>{field.label}</option>)}</select></label>
        <label>Grid <select value={gridMode} onChange={(event) => setGridMode(event.target.value as typeof gridMode)}><option value="off">Off</option><option value="panel_seams">Panel seams</option><option value="all_cells">All cells</option></select></label>
        <label>N (cells/face) <select value={draft.cellsPerPanel} onChange={(event) => setDraft((current) => ({ ...current, cellsPerPanel: Number(event.target.value) }))}>{[2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96].map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
        <label>End time (s) <input type="number" min="0.000001" max="31536000" value={draft.endTimeSeconds} onChange={(event) => setDraft((current) => ({ ...current, endTimeSeconds: Number(event.target.value) }))} /></label>
        <label>Simulation Δt max (s) <input type="number" min="0.000001" max="86400" value={draft.maximumTimeStepSeconds} onChange={(event) => setDraft((current) => ({ ...current, maximumTimeStepSeconds: Number(event.target.value) }))} /></label>
        <label>Display every (steps) <input type="number" min="1" max="1000000" step="1" value={draft.frameIntervalSteps} onChange={(event) => setDraft((current) => ({ ...current, frameIntervalSteps: Number(event.target.value) }))} /></label>
        <label>Amplitude (m) <input type="number" value={amplitude} onChange={(event) => setAmplitude(Number(event.target.value))} /></label>
        <label>Sigma (rad) <input type="number" min="0.000001" step="0.01" value={sigma} onChange={(event) => setSigma(Number(event.target.value))} /></label>
        <label>Mass policy <select value={massPolicy} onChange={(event) => setMassPolicy(event.target.value as typeof massPolicy)}><option value="preserve_global">Preserve global</option><option value="allow_change">Allow change</option></select></label>
        <button type="button" onClick={() => setDraft((current) => undoEdit(current))}>Undo edit</button><button type="button" onClick={() => setDraft((current) => clearEdits(current))}>Clear edits</button><button type="button" onClick={start} disabled={run.state === "submitting" || run.state === "running" || run.state === "cancelling"}>Run</button><button type="button" onClick={() => void controller.cancel()} disabled={run.state !== "running"}>Cancel</button>
      </section>
      <div className="run-tools">
        <p className="engine-status">Engine: {startupToken ? "native C++" : "deterministic demo"}</p>
        <p className="run-status" role="status">Run state: {run.state}{run.error ? ` · ${run.error}` : ""}</p>
        <p className="edit-status" aria-live="polite">Configured edits: {draft.edits.length}{run.request ? ` · current run applied: ${run.request.initialCondition.edits.length}` : ""} · rings show center and sigma</p>
        <section className="timeline" aria-label="Frame timeline"><button type="button" onClick={() => controller.setCurrentFrame(Math.max(0, run.currentFrame - 1))}>Previous</button><button type="button" onClick={togglePlayback} disabled={run.frames.length < 2}>{playing ? "Pause" : run.currentFrame >= run.frames.length - 1 ? "Replay" : "Play"}</button><input aria-label="Frame" type="range" min="0" max={Math.max(0, run.frames.length - 1)} value={run.currentFrame} onChange={(event) => controller.setCurrentFrame(Number(event.target.value))} /><span>{run.frames.length ? `${run.currentFrame + 1}/${run.frames.length} · t=${currentDataset.frame.timeSeconds ?? "?"}s · step ${currentDataset.frame.step ?? "?"}` : "No frames"}</span></section>
        <section className="diagnostics" aria-label="Diagnostics"><button type="button" onClick={() => void downloadBundle()} disabled={!run.runId}>Download bundle</button>{diagnostics.slice(-1).map((sample) => <p key={sample.sequence}>step {sample.step} · mass {sample.mass ?? "n/a"} · energy {sample.energy ?? "n/a"}</p>)}</section>
        <details className="request-panel"><summary>Request</summary><pre>{requestSummary()}</pre></details>
      </div>
    </header>
    <section className="view-grid">
      <article className="view-panel"><h2>2D global map</h2><Map2D dataset={displayed} fieldId={fieldId} gridMode={gridMode} edits={draft.edits} onPick={addFromCell} /></article>
      <article className="view-panel"><h2>3D globe</h2><Globe3D dataset={displayed} fieldId={fieldId} gridMode={gridMode} edits={draft.edits} onPick={addFromCell} /></article>
    </section>
  </main>;
}

createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
