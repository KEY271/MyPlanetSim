import { StrictMode, useCallback, useEffect, useMemo, useRef, useState } from "react";
import { createRoot } from "react-dom/client";
import { addShallowWaterDerivedFields, diagnosticEvents, frameV2Field, frameV2Fields, FrameV2FieldId, frameV2LevelDataset, frameV2Sample, panelNames, parseShallowWaterCsv, UnitVector, VisualFrameV2 } from "@myplanetsim/protocol";
import { DeterministicMockSimulationClient, HttpSimulationClient, PresetDescriptorV1, presetDescriptors, shallowWaterPresetDescriptor } from "@myplanetsim/protocol/client";
import { Globe3D } from "./Globe3D";
import { Map2D } from "./Map2D";
import { ColumnProfile } from "./ColumnProfile";
import { appendEdit, clearEdits, createGaussianEdit, DraftRun, toRunRequest, undoEdit } from "./editor";
import { SimulationController } from "./controller";
import { formatValue } from "./format";
import "./style.css";

function demoDataset(cellsPerPanel: number) {
  const rows = ["panel,i,j,depth_m,momentum_x,momentum_y,momentum_z"];
  for (const panel of ["PX", "PY", "NX", "NY", "PZ", "NZ"]) for (let j = 0; j < cellsPerPanel; j += 1) for (let i = 0; i < cellsPerPanel; i += 1) rows.push(`${panel},${i},${j},3000,0,0,0`);
  return addShallowWaterDerivedFields(parseShallowWaterCsv(rows.join("\n"), cellsPerPanel));
}

const gridSizes = [2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96];

// Cell order is the ADR 0001 panel-major flat order, so the inspector can name the cell a
// reader clicked without a second lookup table.
function cellLabel(cell: number, cellsPerPanel: number) {
  const perPanel = cellsPerPanel * cellsPerPanel;
  const panel = panelNames[Math.floor(cell / perPanel)] ?? "??";
  const remainder = cell % perPanel;
  return `${panel} i=${remainder % cellsPerPanel} j=${Math.floor(remainder / cellsPerPanel)}`;
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
  const [presets, setPresets] = useState<readonly PresetDescriptorV1[]>([shallowWaterPresetDescriptor("rest")]);
  // Swallowing this failure once left the viewer silently stuck on the built-in
  // shallow-water preset with no way to tell that the gateway was never reached, so the
  // error is kept and shown next to the preset selector.
  const [presetError, setPresetError] = useState<string | null>(null);
  useEffect(() => { const unsubscribe = controller.subscribe(setRun); return () => { unsubscribe(); }; }, [controller]);
  useEffect(() => {
    let cancelled = false;
    void controller.capabilities()
      .then((capabilities) => { if (!cancelled) { setPresets(presetDescriptors(capabilities)); setPresetError(null); } })
      .catch((error: unknown) => {
        if (cancelled) return;
        setPresetError(error instanceof Error ? error.message : "the gateway did not return its capabilities");
      });
    return () => { cancelled = true; };
  }, [controller]);
  useEffect(() => {
    if (!playing || run.frames.length < 2) return undefined;
    const timer = window.setInterval(() => {
      const currentFrame = controller.value.currentFrame;
      if (currentFrame >= controller.value.frames.length - 1) { setPlaying(false); return; }
      controller.setCurrentFrame(currentFrame + 1);
    }, 250);
    return () => window.clearInterval(timer);
  }, [controller, playing, run.frames.length]);
  // Shallow-water and dry-hydrostatic field identifiers are disjoint, so each model kind
  // keeps its own selection instead of being coerced into the other model's field list.
  const [shallowFieldId, setShallowFieldId] = useState("depth");
  const [dryFieldId, setDryFieldId] = useState<FrameV2FieldId>("temperature");
  const [level, setLevel] = useState(0);
  const [selectedCell, setSelectedCell] = useState(0);
  const [gridMode, setGridMode] = useState<"off" | "panel_seams" | "all_cells">("panel_seams");
  const [amplitude, setAmplitude] = useState(100);
  const [sigma, setSigma] = useState(0.2);
  const [massPolicy, setMassPolicy] = useState<"preserve_global" | "allow_change">("preserve_global");
  const [message, setMessage] = useState("Choose a cell in either view to create a pending edit.");
  const editCounter = useRef(1);
  const preset = presets.find((candidate) => candidate.id === draft.presetId) ?? shallowWaterPresetDescriptor(draft.presetId);
  const editable = preset.supportedEdits.includes("gaussian_depth");
  // Which model is on screen follows the decoded frame, not the draft preset: switching the
  // preset selector must not reinterpret the frames of the run that is already displayed.
  const currentFrame = run.frames[run.currentFrame];
  const dryFrame: VisualFrameV2 | null = currentFrame?.schemaVersion === 2 ? currentFrame : null;
  // The level selector and the profile need a decoded column, so they cannot appear before
  // the first frame. Announce them from the preset descriptor instead of leaving the panel
  // unchanged, which read as the feature being absent.
  const dryPending = preset.levels !== null && dryFrame === null;
  const dryField = dryFrame ? frameV2Field(dryFieldId) : null;
  const fieldId = dryFrame ? dryFieldId : shallowFieldId;
  const activeLevel = dryFrame ? Math.min(level, dryFrame.levels - 1) : 0;
  const activeCell = dryFrame ? Math.min(selectedCell, 6 * dryFrame.cellsPerPanel ** 2 - 1) : selectedCell;
  const currentDataset = useMemo(() => {
    if (dryFrame && dryField) return frameV2LevelDataset(dryFrame, dryField.id, dryField.volume ? activeLevel : 0);
    return currentFrame?.schemaVersion === 1 ? currentFrame : dataset;
  }, [activeLevel, currentFrame, dataset, dryField, dryFrame]);
  const addFromCell = useCallback((cell: number | null, origin?: UnitVector) => {
    if (cell === null || !origin) return;
    if (dryFrame) { setSelectedCell(cell); setMessage(`Selected column ${cellLabel(cell, dryFrame.cellsPerPanel)}.`); return; }
    if (!editable) { setMessage(`${preset.id} does not accept initial-condition edits.`); return; }
    // The updater must stay pure: React may replay it, so the identifier and the status
    // message are derived once here instead of inside setDraft.
    const edit = createGaussianEdit(origin, amplitude, sigma, massPolicy, `edit-${editCounter.current++}`);
    setDraft((current) => appendEdit(current, edit));
    setMessage(`${edit.id} added at the clicked location (cell ${cell}); the ring shows sigma.`);
  }, [amplitude, dryFrame, editable, massPolicy, preset.id, sigma]);
  const selectPreset = (presetId: string) => {
    const next = presets.find((candidate) => candidate.id === presetId) ?? shallowWaterPresetDescriptor(presetId);
    // A preset without a column must carry no levels at all, and one with a column starts
    // from its own default rather than from whatever the previous preset used.
    setDraft((current) => ({ ...current, presetId, edits: [],
      cellsPerPanel: Math.min(current.cellsPerPanel, next.maximumCellsPerPanel),
      levels: next.levels ?? undefined }));
    setMessage(next.modelKind === "dry_hydrostatic"
      ? `${presetId} runs the dry hydrostatic core; clicking a cell selects a column.`
      : `${presetId} runs the shallow-water core; clicking a cell adds an edit.`);
  };
  const requestSummary = () => { try { return JSON.stringify(toRunRequest(draft), null, 2); } catch (error) { return error instanceof Error ? error.message : "invalid request"; } };
  const start = () => { try { void controller.start(toRunRequest(draft)); } catch (error) { setMessage(error instanceof Error ? error.message : "invalid request"); } };
  const togglePlayback = () => { if (playing) { setPlaying(false); return; } if (run.currentFrame >= run.frames.length - 1) controller.setCurrentFrame(0); setPlaying(true); };
  const displayed = run.frames.length > 0 ? currentDataset : dataset;
  const diagnostics = diagnosticEvents(run.events);
  const sample = dryFrame && dryField ? frameV2Sample(dryFrame, dryField.id, activeCell, activeLevel) : null;
  const timeSeconds = dryFrame ? dryFrame.timeSeconds : currentDataset.frame.timeSeconds;
  const step = dryFrame ? dryFrame.step : currentDataset.frame.step;
  const downloadBundle = async () => { const bundle = await controller.bundle(); if (!bundle) return; const link = document.createElement("a"); link.href = URL.createObjectURL(new Blob([JSON.stringify(bundle, null, 2)], { type: "application/json" })); link.download = `${bundle.runId}.json`; link.click(); URL.revokeObjectURL(link.href); };
  return <main className="app-shell">
    <header className="top-panel">
      <div className="title-row"><p className="eyebrow">MyPlanetSim · Phase 5</p><h1>Interactive cubed-sphere visualizer</h1><p className="message">{message}</p></div>
      <section className="control-panel" aria-label="Run and view controls">
        <label>Preset <select value={draft.presetId} onChange={(event) => selectPreset(event.target.value)}>{presets.map((candidate) => <option key={candidate.id} value={candidate.id}>{candidate.id}</option>)}</select></label>
        <label>Field <select value={fieldId} onChange={(event) => dryFrame ? setDryFieldId(event.target.value as FrameV2FieldId) : setShallowFieldId(event.target.value)}>
          {dryFrame
            ? frameV2Fields.map((field) => <option key={field.id} value={field.id}>{field.label}</option>)
            : dataset.fields.map((field) => <option key={field.id} value={field.id}>{field.label}</option>)}
        </select></label>
        {dryFrame ? <label>Model level (0 = top) <input aria-label="Model level" type="range" min="0" max={dryFrame.levels - 1} step="1" value={activeLevel} disabled={!dryField?.volume} onChange={(event) => setLevel(Number(event.target.value))} /></label> : null}
        {dryPending ? <label>Model level (0 = top) <input aria-label="Model level" type="range" min="0" max={(draft.levels ?? preset.levels ?? 1) - 1} step="1" value={0} disabled readOnly /></label> : null}
        <label>Grid <select value={gridMode} onChange={(event) => setGridMode(event.target.value as typeof gridMode)}><option value="off">Off</option><option value="panel_seams">Panel seams</option><option value="all_cells">All cells</option></select></label>
        <label>N (cells/face) <select value={draft.cellsPerPanel} onChange={(event) => setDraft((current) => ({ ...current, cellsPerPanel: Number(event.target.value) }))}>{gridSizes.filter((value) => value <= preset.maximumCellsPerPanel).map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
        {preset.levels === null ? null : <label>K (vertical resolution) <input aria-label="K vertical resolution" type="number" min="1" step="1" max={preset.maximumLevels ?? preset.levels} value={draft.levels ?? preset.levels} onChange={(event) => setDraft((current) => ({ ...current, levels: Number(event.target.value) }))} /></label>}
        <label>End time (s) <input type="number" min="0.000001" max="31536000" value={draft.endTimeSeconds} onChange={(event) => setDraft((current) => ({ ...current, endTimeSeconds: Number(event.target.value) }))} /></label>
        <label>Simulation Δt max (s) <input type="number" min="0.000001" max="86400" value={draft.maximumTimeStepSeconds} onChange={(event) => setDraft((current) => ({ ...current, maximumTimeStepSeconds: Number(event.target.value) }))} /></label>
        <label>Display every (steps) <input type="number" min="1" max="1000000" step="1" value={draft.frameIntervalSteps} onChange={(event) => setDraft((current) => ({ ...current, frameIntervalSteps: Number(event.target.value) }))} /></label>
        {editable ? <>
          <label>Amplitude (m) <input type="number" value={amplitude} onChange={(event) => setAmplitude(Number(event.target.value))} /></label>
          <label>Sigma (rad) <input type="number" min="0.000001" step="0.01" value={sigma} onChange={(event) => setSigma(Number(event.target.value))} /></label>
          <label>Mass policy <select value={massPolicy} onChange={(event) => setMassPolicy(event.target.value as typeof massPolicy)}><option value="preserve_global">Preserve global</option><option value="allow_change">Allow change</option></select></label>
          <button type="button" onClick={() => setDraft((current) => undoEdit(current))}>Undo edit</button>
          <button type="button" onClick={() => setDraft((current) => clearEdits(current))}>Clear edits</button>
        </> : null}
        <button type="button" onClick={start} disabled={run.state === "submitting" || run.state === "running" || run.state === "cancelling"}>Run</button><button type="button" onClick={() => void controller.cancel()} disabled={run.state !== "running"}>Cancel</button>
      </section>
      <div className="run-tools">
        <p className="engine-status">Engine: {startupToken ? "native C++" : "deterministic demo"} · model: {preset.modelKind} · frame schema: {dryFrame ? 2 : 1}</p>
        {presetError ? <p className="preset-error" role="alert">Preset list unavailable ({presetError}); showing the built-in shallow-water preset only.</p> : null}
        {dryPending ? <p className="preset-hint">{preset.id} publishes {draft.levels ?? preset.levels} model levels. Run it to enable the level selector, the column profile, and the cell inspector.</p> : null}
        <p className="run-status" role="status">Run state: {run.state}{run.error ? ` · ${run.error}` : ""}</p>
        <p className="edit-status" aria-live="polite">{editable
          ? `Configured edits: ${draft.edits.length}${run.request ? ` · current run applied: ${run.request.initialCondition.edits.length}` : ""} · rings show center and sigma`
          : "This preset has no initial-condition edits; clicking a cell selects a column."}</p>
        <section className="timeline" aria-label="Frame timeline"><button type="button" onClick={() => controller.setCurrentFrame(Math.max(0, run.currentFrame - 1))}>Previous</button><button type="button" onClick={togglePlayback} disabled={run.frames.length < 2}>{playing ? "Pause" : run.currentFrame >= run.frames.length - 1 ? "Replay" : "Play"}</button><input aria-label="Frame" type="range" min="0" max={Math.max(0, run.frames.length - 1)} value={run.currentFrame} onChange={(event) => controller.setCurrentFrame(Number(event.target.value))} /><span>{run.frames.length ? `${run.currentFrame + 1}/${run.frames.length} · t=${timeSeconds ?? "?"}s · step ${step ?? "?"}` : "No frames"}</span></section>
        {dryFrame && dryField && sample ? <section className="inspector" aria-label="Column inspector">
          <p>Cell {activeCell} ({cellLabel(activeCell, dryFrame.cellsPerPanel)}){dryField.volume ? ` · level ${activeLevel}/${dryFrame.levels - 1}` : " · surface"}</p>
          <p>{dryField.label}: {formatValue(sample.value)} {dryField.unit}{sample.pressurePa === null ? "" : ` · pressure ${formatValue(sample.pressurePa)} Pa`}</p>
          <p>t={dryFrame.timeSeconds}s · step {dryFrame.step} · fingerprint {dryFrame.configFingerprint}</p>
        </section> : null}
        <section className="diagnostics" aria-label="Diagnostics"><button type="button" onClick={() => void downloadBundle()} disabled={!run.runId}>Download bundle</button>{diagnostics.slice(-1).map((sampleEvent) => <p key={sampleEvent.sequence}>step {sampleEvent.step} · mass {sampleEvent.mass ?? "n/a"} · energy {sampleEvent.energy ?? "n/a"}</p>)}</section>
        <details className="request-panel"><summary>Request</summary><pre>{requestSummary()}</pre></details>
      </div>
    </header>
    <section className={dryFrame && dryField?.volume ? "view-grid with-profile" : "view-grid"}>
      <article className="view-panel"><h2>2D global map</h2><Map2D dataset={displayed} fieldId={fieldId} gridMode={gridMode} edits={draft.edits} onPick={addFromCell} /></article>
      <article className="view-panel"><h2>3D globe</h2><Globe3D dataset={displayed} fieldId={fieldId} gridMode={gridMode} edits={draft.edits} onPick={addFromCell} /></article>
      {dryFrame && dryField?.volume
        ? <ColumnProfile frame={dryFrame} field={dryField} cell={activeCell} level={activeLevel} onSelectLevel={setLevel} />
        : null}
    </section>
  </main>;
}

createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
