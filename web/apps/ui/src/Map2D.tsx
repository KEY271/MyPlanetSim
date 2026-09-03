import { useEffect, useRef, useState } from "react";
import { UnitVector, VisualDatasetV1, cubedSphereCells, gridEdges } from "@myplanetsim/protocol";
import { drawWrappedSegment, inverseProject, projectUnit } from "./projection";

interface Map2DProps {
  readonly dataset: VisualDatasetV1;
  readonly fieldId: string;
  readonly gridMode: "off" | "panel_seams" | "all_cells";
  readonly onPick?: (cell: number | null) => void;
  readonly pendingOrigin?: UnitVector;
}

export function Map2D({ dataset, fieldId, gridMode, onPick, pendingOrigin }: Map2DProps) {
  const host = useRef<HTMLDivElement>(null);
  const base = useRef<HTMLCanvasElement>(null);
  const grid = useRef<HTMLCanvasElement>(null);
  const overlay = useRef<HTMLCanvasElement>(null);
  const [view, setView] = useState({ zoom: 1, pan: [0, 0] as [number, number] });
  useEffect(() => {
    if (!host.current || !base.current || !grid.current || !overlay.current) return undefined;
    const canvases = [base.current, grid.current, overlay.current];
    const draw = (width: number, height: number) => {
      const field = dataset.fields.find((candidate) => candidate.id === fieldId);
      if (!field) return;
      const baseContext = base.current?.getContext("2d"); const gridContext = grid.current?.getContext("2d"); const overlayContext = overlay.current?.getContext("2d");
      if (!baseContext || !gridContext || !overlayContext) return;
      baseContext.clearRect(0, 0, width, height); gridContext.clearRect(0, 0, width, height);
      overlayContext.clearRect(0, 0, width, height);
      const cells = cubedSphereCells(dataset.grid.cellsPerPanel);
      const minimum = Math.min(...field.values); const range = Math.max(...field.values) - minimum || 1;
      for (let index = 0; index < cells.length; index += 1) {
        const polygon = cells[index].corners.map((corner) => projectUnit(corner, width, height, view.zoom, view.pan));
        baseContext.beginPath(); baseContext.moveTo(...polygon[0]); for (const point of polygon.slice(1)) baseContext.lineTo(...point); baseContext.closePath();
        baseContext.fillStyle = `hsl(${237 - 237 * (field.values[index] - minimum) / range} 80% 50%)`; baseContext.fill();
      }
      if (gridMode !== "off") {
        gridContext.beginPath(); gridContext.strokeStyle = "rgba(255,255,255,.65)"; gridContext.lineWidth = gridMode === "panel_seams" ? 2 : 0.75;
        for (const edge of gridEdges(dataset.grid.cellsPerPanel)) { if (gridMode === "panel_seams" && !edge.panelSeam) continue; drawWrappedSegment(gridContext, edge.first, edge.second, width, height, view.zoom, view.pan); }
        gridContext.stroke();
      }
      if (pendingOrigin) { const marker = projectUnit(pendingOrigin, width, height, view.zoom, view.pan); overlayContext.fillStyle = "#ffcf56"; overlayContext.beginPath(); overlayContext.arc(marker[0], marker[1], 6, 0, 2 * Math.PI); overlayContext.fill(); }
    };
    const resize = () => {
      const bounds = host.current?.getBoundingClientRect(); if (!bounds) return;
      for (const canvas of canvases) { const ratio = Math.min(window.devicePixelRatio, 2); canvas.width = Math.max(1, Math.round(bounds.width * ratio)); canvas.height = Math.max(1, Math.round(bounds.height * ratio)); canvas.style.width = `${bounds.width}px`; canvas.style.height = `${bounds.height}px`; }
      draw(bounds.width, bounds.height);
    };
    const observer = new ResizeObserver(resize); observer.observe(host.current); resize();
    return () => observer.disconnect();
  }, [dataset, fieldId, gridMode, view, pendingOrigin]);
  useEffect(() => {
    const element = host.current; if (!element) return undefined;
    let dragging = false; let moved = false; let last: [number, number] = [0, 0];
    const down = (event: PointerEvent) => { dragging = true; moved = false; last = [event.clientX, event.clientY]; element.setPointerCapture(event.pointerId); };
    const move = (event: PointerEvent) => { if (!dragging) return; const dx = event.clientX - last[0]; const dy = event.clientY - last[1]; if (Math.abs(dx) + Math.abs(dy) > 2) moved = true; setView((current) => ({ ...current, pan: [current.pan[0] + dx, current.pan[1] + dy] })); last = [event.clientX, event.clientY]; };
    const up = (event: PointerEvent) => { if (!moved) { const bounds = element.getBoundingClientRect(); const point = inverseProject(event.clientX - bounds.left, event.clientY - bounds.top, bounds.width, bounds.height, view.zoom, view.pan); const cells = cubedSphereCells(dataset.grid.cellsPerPanel); const picked = cells.findIndex((cell) => cell.center[0] * point[0] + cell.center[1] * point[1] + cell.center[2] * point[2] > 0.999); onPick?.(picked >= 0 ? picked : null); } dragging = false; };
    const wheel = (event: WheelEvent) => setView((current) => ({ ...current, zoom: Math.max(1, Math.min(8, current.zoom * (event.deltaY < 0 ? 1.1 : 0.9))) }));
    element.addEventListener("pointerdown", down); element.addEventListener("pointermove", move); element.addEventListener("pointerup", up); element.addEventListener("wheel", wheel, { passive: true });
    return () => { element.removeEventListener("pointerdown", down); element.removeEventListener("pointermove", move); element.removeEventListener("pointerup", up); element.removeEventListener("wheel", wheel); };
  }, [dataset, onPick, view]);
  return <div ref={host} className="map-2d" aria-label="2D global map"><canvas ref={base} /><canvas ref={grid} /><canvas ref={overlay} /><button type="button" onClick={() => setView({ zoom: 1, pan: [0, 0] })}>Reset map</button></div>;
}
