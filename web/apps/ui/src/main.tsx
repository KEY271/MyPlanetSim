import { StrictMode, useEffect, useMemo, useRef, useState } from "react";
import { createRoot } from "react-dom/client";
import {
  DatasetField,
  DatasetPeriod,
  fieldRange,
  panelNames,
  periodField,
  periodLevelSlice,
  PeriodDataset,
  VisualDatasetReader,
  VisualDatasetV1,
} from "@myplanetsim/protocol";
import { Globe3D } from "./Globe3D";
import { Map2D } from "./Map2D";
import { ColumnProfile } from "./ColumnProfile";
import { browserDatasetFiles } from "./dataset-files";
import { formatValue } from "./format";
import "./style.css";

type GridMode = "off" | "panel_seams" | "all_cells";

function cellLabel(cell: number, cellsPerPanel: number) {
  const perPanel = cellsPerPanel * cellsPerPanel;
  const panel = panelNames[Math.floor(cell / perPanel)] ?? "??";
  const remainder = cell % perPanel;
  return `${panel} i=${remainder % cellsPerPanel} j=${Math.floor(remainder / cellsPerPanel)}`;
}

function label(id: string) {
  return id.split(/[._]/).map((word) => word.charAt(0).toUpperCase() + word.slice(1)).join(" ");
}

function renderDataset(reader: VisualDatasetReader, field: DatasetField | null,
                       period: PeriodDataset | null, terrainField: string | null,
                       level: number): VisualDatasetV1 {
  const id = field?.id ?? terrainField!;
  const unit = field?.unit ??
    reader.manifest.terrain.fields.find((candidate) => candidate.id === terrainField)!.unit;
  const values = field && period
    ? periodLevelSlice(reader.manifest, period, field.id, level)
    : reader.terrain.fields.get(terrainField!)!;
  return Object.freeze({
    schemaVersion: 1,
    sourceKind: "visual_dataset",
    grid: Object.freeze({ topology: "cubed_sphere", mapping: "equiangular_gnomonic_v1",
      cellsPerPanel: reader.manifest.cellsPerPanel,
      flattenOrder: "panel_major_then_j_then_i" }),
    frame: Object.freeze({ timeSeconds: period?.descriptor.actualEndS ?? null,
      step: null, configFingerprint: reader.manifest.configFingerprint }),
    fields: Object.freeze([Object.freeze({ id, label: label(id), unit, kind: "scalar" as const,
      provenance: "VisualDatasetV1", values })]),
  });
}

function App() {
  const loadGeneration = useRef(0);
  const [reader, setReader] = useState<VisualDatasetReader | null>(null);
  const [period, setPeriod] = useState<PeriodDataset | null>(null);
  const [periodIndex, setPeriodIndex] = useState<number | null>(null);
  const [fieldId, setFieldId] = useState("elevation");
  const [level, setLevel] = useState(0);
  const [selectedCell, setSelectedCell] = useState(0);
  const [gridMode, setGridMode] = useState<GridMode>("panel_seams");
  const [message, setMessage] = useState("Open a VisualDatasetV1 folder to inspect terrain and period means.");

  const folderInput = (node: HTMLInputElement | null) => {
    node?.setAttribute("webkitdirectory", "");
  };

  const atmosphericField = reader?.manifest.fields.find((field) => field.id === fieldId) ?? null;
  const terrainField = reader?.manifest.terrain.fields.find((field) => field.id === fieldId) ?? null;
  const activePeriodIndex = periodIndex ?? reader?.manifest.periods[0]?.index ?? null;

  useEffect(() => {
    if (!reader || !atmosphericField || activePeriodIndex === null) { setPeriod(null); return; }
    const generation = ++loadGeneration.current;
    setMessage("Loading selected period…");
    void reader.period(activePeriodIndex).then((loaded) => {
      if (generation !== loadGeneration.current) return;
      setPeriod(loaded);
      setMessage(loaded.descriptor.complete ? "Complete period loaded." :
        "Partial period loaded; coverage is shown below.");
    }).catch((error: unknown) => {
      if (generation !== loadGeneration.current) return;
      setPeriod(null);
      setMessage(error instanceof Error ? error.message : "Unable to load period.");
    });
  }, [activePeriodIndex, atmosphericField, reader]);

  const displayed = useMemo(() => {
    if (!reader || (!terrainField && (!atmosphericField || !period))) return null;
    return renderDataset(reader, atmosphericField, period, terrainField?.id ?? null,
      atmosphericField?.location === "atmosphere" ? level : 0);
  }, [atmosphericField, level, period, reader, terrainField]);

  const openDataset = async (files: FileList | null) => {
    if (!files?.length) return;
    const generation = ++loadGeneration.current;
    try {
      const opened = await VisualDatasetReader.open(browserDatasetFiles(files));
      if (generation !== loadGeneration.current) return;
      setReader(opened); setPeriod(null); setPeriodIndex(opened.manifest.periods[0]?.index ?? null);
      setFieldId(opened.manifest.periods.length && opened.manifest.fields.some((field) => field.id === "temperature")
        ? "temperature" : "elevation");
      setLevel(0); setSelectedCell(0);
      setMessage(`Opened N=${opened.manifest.cellsPerPanel}, K=${opened.manifest.levelCount} dataset.`);
    } catch (error) {
      if (generation !== loadGeneration.current) return;
      setReader(null); setPeriod(null);
      setMessage(error instanceof Error ? error.message : "Unable to open dataset.");
    }
  };

  if (!reader || !displayed) return <main className="empty-state">
    <p className="eyebrow">MyPlanetSim · Phase 15</p>
    <h1>Terrain & period-mean viewer</h1>
    <p>{message}</p>
    <label className="file-button">Open dataset folder
      <input ref={folderInput} type="file" multiple onChange={(event) => void openDataset(event.target.files)} />
    </label>
    <p className="file-hint">Select the folder containing manifest.json, terrain.bin, and means/.</p>
  </main>;

  const selectedField = atmosphericField ??
    ({ id: terrainField!.id, unit: terrainField!.unit, location: "surface",
      valueOffset: 0, valueCount: reader.manifest.cellCount } satisfies DatasetField);
  const activeLevel = selectedField.location === "atmosphere"
    ? Math.min(level, reader.manifest.levelCount - 1) : 0;
  const activeCell = Math.min(selectedCell, reader.manifest.cellCount - 1);
  const shownValues = displayed.fields[0].values;
  const [minimum, maximum] = fieldRange(displayed.fields[0]);
  const sample = shownValues[activeCell];
  const pressureDescriptor = reader.manifest.fields.find((field) => field.id === "pressure");
  const volumeValues = atmosphericField && period && atmosphericField.location === "atmosphere"
    ? periodField(reader.manifest, period, atmosphericField.id) : null;
  const pressureValues = pressureDescriptor && period
    ? periodField(reader.manifest, period, pressureDescriptor.id) : null;
  const profileValues = volumeValues
    ? volumeValues.slice(activeCell * reader.manifest.levelCount,
      (activeCell + 1) * reader.manifest.levelCount) : null;
  const profilePressure = pressureValues
    ? pressureValues.slice(activeCell * reader.manifest.levelCount,
      (activeCell + 1) * reader.manifest.levelCount) : null;
  const periodDescriptor: DatasetPeriod | null = period?.descriptor ?? null;
  const pick = (cell: number | null) => {
    if (cell === null) return;
    setSelectedCell(cell);
    setMessage(`Selected ${cellLabel(cell, reader.manifest.cellsPerPanel)}.`);
  };

  return <main className="app-shell">
    <header className="top-panel">
      <div className="title-row"><p className="eyebrow">MyPlanetSim · Phase 15</p>
        <h1>Terrain & period-mean viewer</h1><p className="message">{message}</p></div>
      <section className="control-panel" aria-label="Dataset view controls">
        <label>Dataset <span className="readonly-value">N={reader.manifest.cellsPerPanel} · K={reader.manifest.levelCount}</span></label>
        <label>Field <select value={fieldId} onChange={(event) => setFieldId(event.target.value)}>
          <optgroup label="Terrain">{reader.manifest.terrain.fields
            .filter((field) => !field.id.startsWith("center_"))
            .map((field) => <option key={field.id} value={field.id}>{label(field.id)} ({field.unit})</option>)}</optgroup>
          {reader.manifest.periods.length ? <optgroup label="Period means">{reader.manifest.fields
            .map((field) => <option key={field.id} value={field.id}>{label(field.id)} ({field.unit})</option>)}</optgroup> : null}
        </select></label>
        <label>Period <select value={activePeriodIndex ?? ""} disabled={!atmosphericField}
          onChange={(event) => setPeriodIndex(Number(event.target.value))}>
          {reader.manifest.periods.map((candidate) => <option key={candidate.index} value={candidate.index}>
            {candidate.index}: {candidate.scheduledStartS}–{candidate.scheduledEndS}s
            {candidate.complete ? "" : " (partial)"}</option>)}
        </select></label>
        <label>Model level (0 = top) <input aria-label="Model level" type="range" min="0"
          max={reader.manifest.levelCount - 1} step="1" value={activeLevel}
          disabled={selectedField.location === "surface"} onChange={(event) => setLevel(Number(event.target.value))} /></label>
        <label>Grid <select value={gridMode} onChange={(event) => setGridMode(event.target.value as GridMode)}>
          <option value="off">Off</option><option value="panel_seams">Panel seams</option>
          <option value="all_cells">All cells</option></select></label>
        <label className="file-button compact">Open another
          <input type="file" multiple ref={folderInput} onChange={(event) => void openDataset(event.target.files)} /></label>
      </section>
      <div className="run-tools">
        <p className="engine-status">Read-only local dataset · fingerprint {reader.manifest.configFingerprint}</p>
        {periodDescriptor ? <p className={periodDescriptor.complete ? "period-status" : "period-status partial"}>
          {periodDescriptor.actualStartS}–{periodDescriptor.actualEndS}s ·
          {" "}{(periodDescriptor.coverage * 100).toFixed(1)}% coverage
        </p> : <p className="period-status">Static terrain</p>}
        <p className="color-range">{label(selectedField.id)} · {formatValue(minimum)}–{formatValue(maximum)} {selectedField.unit}</p>
        <section className="inspector" aria-label="Cell inspector">
          <p>{cellLabel(activeCell, reader.manifest.cellsPerPanel)} ·
            {selectedField.location === "atmosphere" ? ` level ${activeLevel}` : " surface"} ·
            {" "}{formatValue(sample)} {selectedField.unit}</p>
        </section>
      </div>
    </header>
    <section className={profileValues && profilePressure ? "view-grid with-profile" : "view-grid"}>
      <article className="view-panel"><h2>2D global map</h2>
        <Map2D dataset={displayed} fieldId={selectedField.id} gridMode={gridMode} onPick={pick} /></article>
      <article className="view-panel"><h2>3D globe</h2>
        <Globe3D dataset={displayed} fieldId={selectedField.id} gridMode={gridMode} onPick={pick} /></article>
      {profileValues && profilePressure ? <ColumnProfile values={profileValues}
        pressurePa={profilePressure} label={label(selectedField.id)} unit={selectedField.unit}
        cell={activeCell} level={activeLevel} onSelectLevel={setLevel} /> : null}
    </section>
  </main>;
}

createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
