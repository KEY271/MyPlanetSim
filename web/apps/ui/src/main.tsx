import { StrictMode, useEffect, useMemo, useState } from "react";
import { createRoot } from "react-dom/client";
import { addShallowWaterDerivedFields, diagnosticEvents, parseShallowWaterCsv, UnitVector, cubedSphereCells } from "@myplanetsim/protocol";
import { DeterministicMockSimulationClient } from "@myplanetsim/protocol/client";
import { Globe3D } from "./Globe3D";
import { Map2D } from "./Map2D";
import { appendEdit, clearEdits, createGaussianEdit, DraftRun, toRunRequest, undoEdit } from "./editor";
import { SimulationController } from "./controller";
import "./style.css";

function demoDataset() {
  const rows = ["panel,i,j,depth_m,momentum_x,momentum_y,momentum_z"];
  for (const panel of ["PX", "PY", "NX", "NY", "PZ", "NZ"]) for (let j = 0; j < 2; j += 1) for (let i = 0; i < 2; i += 1) rows.push(`${panel},${i},${j},3000,0,0,0`);
  return addShallowWaterDerivedFields(parseShallowWaterCsv(rows.join("\n"), 2));
}

function App() {
  const dataset = useMemo(demoDataset, []);
  const controller = useMemo(() => new SimulationController(new DeterministicMockSimulationClient()), []);
  const [run, setRun] = useState(controller.value);
  const [playing, setPlaying] = useState(false);
  useEffect(() => { const unsubscribe = controller.subscribe(setRun); return () => { unsubscribe(); }; }, [controller]);
  useEffect(() => {
    if (!playing || run.frames.length < 2) return undefined;
    const timer = window.setInterval(() => controller.setCurrentFrame((controller.value.currentFrame + 1) % controller.value.frames.length), 250);
    return () => window.clearInterval(timer);
  }, [controller, playing, run.frames.length, run.currentFrame]);
  const [fieldId, setFieldId] = useState("depth");
  const [gridMode, setGridMode] = useState<"off" | "panel_seams" | "all_cells">("panel_seams");
  const [selectedCell, setSelectedCell] = useState<number | null>(null);
  const [amplitude, setAmplitude] = useState(100);
  const [sigma, setSigma] = useState(0.2);
  const [massPolicy, setMassPolicy] = useState<"preserve_global" | "allow_change">("preserve_global");
  const [draft, setDraft] = useState<DraftRun>({ presetId: "rest", endTimeSeconds: 3600, maximumTimeStepSeconds: 60, frameIntervalSteps: 10, edits: [] });
  const [message, setMessage] = useState("Choose a cell in either view to create a pending edit.");
  const currentDataset = run.frames[run.currentFrame] ?? dataset;
  const selectedOrigin: UnitVector | undefined = selectedCell === null ? undefined : cubedSphereCells(currentDataset.grid.cellsPerPanel)[selectedCell]?.center;
  const addFromCell = (cell: number | null) => { if (cell === null) return; setSelectedCell(cell); const origin = cubedSphereCells(dataset.grid.cellsPerPanel)[cell].center; const edit = createGaussianEdit(origin, amplitude, sigma, massPolicy, `edit-${draft.edits.length + 1}`); setDraft((current) => appendEdit(current, edit)); setMessage(`Pending edit ${edit.id} is ready; native C++ remains authoritative after submit.`); };
  const requestSummary = () => { try { return JSON.stringify(toRunRequest(draft), null, 2); } catch (error) { return error instanceof Error ? error.message : "invalid request"; } };
  const start = () => { try { void controller.start(toRunRequest(draft)); } catch (error) { setMessage(error instanceof Error ? error.message : "invalid request"); } };
  const displayed = run.frames.length > 0 ? currentDataset : dataset;
  const diagnostics = diagnosticEvents(run.events);
  const downloadBundle = async () => { const bundle = await controller.bundle(); if (!bundle) return; const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([JSON.stringify(bundle, null, 2)], { type: "application/json" })); link.download = `${bundle.runId}.json`; link.click(); URL.revokeObjectURL(link.href); };
  return <main className="app-shell">
    <p className="eyebrow">MyPlanetSim · Phase 3</p><h1>Interactive shallow-water visualizer</h1><p>{message}</p>
    <section className="control-panel" aria-label="Initial condition editor">
      <label>Field <select value={fieldId} onChange={(event) => setFieldId(event.target.value)}>{dataset.fields.map((field) => <option key={field.id} value={field.id}>{field.label}</option>)}</select></label>
      <label>Grid <select value={gridMode} onChange={(event) => setGridMode(event.target.value as typeof gridMode)}><option value="off">Off</option><option value="panel_seams">Panel seams</option><option value="all_cells">All cells</option></select></label>
      <label>Amplitude (m) <input type="number" value={amplitude} onChange={(event) => setAmplitude(Number(event.target.value))} /></label>
      <label>Sigma (rad) <input type="number" min="0.000001" step="0.01" value={sigma} onChange={(event) => setSigma(Number(event.target.value))} /></label>
      <label>Mass policy <select value={massPolicy} onChange={(event) => setMassPolicy(event.target.value as typeof massPolicy)}><option value="preserve_global">Preserve global</option><option value="allow_change">Allow change</option></select></label>
      <button type="button" onClick={() => setDraft((current) => undoEdit(current))}>Undo edit</button><button type="button" onClick={() => setDraft((current) => clearEdits(current))}>Clear edits</button><button type="button" onClick={start} disabled={run.state === "submitting" || run.state === "running" || run.state === "cancelling"}>Run</button><button type="button" onClick={() => void controller.cancel()} disabled={run.state !== "running"}>Cancel</button>
    </section>
    <p role="status">Run state: {run.state}{run.error ? ` · ${run.error}` : ""}</p>
    <section className="view-grid"><div><h2>3D globe</h2><Globe3D dataset={displayed} fieldId={fieldId} gridMode={gridMode} pendingOrigin={run.frames.length === 0 ? selectedOrigin : undefined} onPick={addFromCell} /></div><div><h2>2D global map</h2><Map2D dataset={displayed} fieldId={fieldId} gridMode={gridMode} pendingOrigin={run.frames.length === 0 ? selectedOrigin : undefined} onPick={addFromCell} /></div></section>
    <section className="timeline"><button type="button" onClick={() => controller.setCurrentFrame(Math.max(0, run.currentFrame - 1))}>Previous frame</button><button type="button" onClick={() => setPlaying((value) => !value)} disabled={run.frames.length < 2}>{playing ? "Pause" : "Play"}</button><input aria-label="Frame" type="range" min="0" max={Math.max(0, run.frames.length - 1)} value={run.currentFrame} onChange={(event) => controller.setCurrentFrame(Number(event.target.value))} /><span>{run.frames.length ? `${run.currentFrame + 1}/${run.frames.length}` : "No frames"}</span></section>
    <section className="diagnostics"><h2>Diagnostics</h2><button type="button" onClick={() => void downloadBundle()} disabled={!run.runId}>Download run bundle</button>{diagnostics.map((sample) => <p key={sample.sequence}>step {sample.step} · mass {sample.mass ?? "n/a"} · energy {sample.energy ?? "n/a"}</p>)}</section>
    <section className="request-panel"><h2>Request summary</h2><pre>{requestSummary()}</pre></section>
  </main>;
}

createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
