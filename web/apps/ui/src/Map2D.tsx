import { useEffect, useRef, useState } from "react";
import { geoPath } from "d3-geo";
import { fieldRange, hitTest, UnitVector, VisualDatasetV1, gridEdges, panelNames } from "@myplanetsim/protocol";
import { equirectangularProjection, inverseProject, isPointInsideProjectedMap, unitToCellIndex, unitToLongitudeLatitude } from "./projection";

interface MapView { readonly zoom: number; readonly pan: readonly [number, number] }

interface Map2DProps {
  readonly dataset: VisualDatasetV1;
  readonly fieldId: string;
  readonly gridMode: "off" | "panel_seams" | "all_cells";
  readonly onPick?: (cell: number | null, origin?: UnitVector) => void;
}

export function Map2D({ dataset, fieldId, gridMode, onPick }: Map2DProps) {
  const host = useRef<HTMLDivElement>(null);
  const base = useRef<HTMLCanvasElement>(null);
  const grid = useRef<HTMLCanvasElement>(null);
  const [view, setView] = useState<MapView>({ zoom: 1, pan: [0, 0] });
  // Draw inputs live in refs so a pan or a new frame schedules one coalesced repaint
  // instead of tearing down the effect and repainting synchronously per pointer event.
  const scene = useRef({ dataset, fieldId, gridMode, view, quality: 1 });
  scene.current = { ...scene.current, dataset, fieldId, gridMode, view };
  const onPickRef = useRef(onPick); onPickRef.current = onPick;
  const frameRequest = useRef(0);
  const size = useRef<readonly [number, number]>([0, 0]);
  const requestDraw = useRef(() => {});
  useEffect(() => {
    if (!host.current || !base.current || !grid.current) return undefined;
    const canvases = [base.current, grid.current];
    const draw = () => {
      const [width, height] = size.current;
      const { dataset, fieldId, gridMode, view, quality } = scene.current;
      if (width <= 0 || height <= 0) return;
      const field = dataset.fields.find((candidate) => candidate.id === fieldId);
      if (!field) return;
      const baseContext = base.current?.getContext("2d"); const gridContext = grid.current?.getContext("2d");
      if (!baseContext || !gridContext) return;
      baseContext.clearRect(0, 0, width, height); gridContext.clearRect(0, 0, width, height);
      const [minimum, maximum] = fieldRange(field); const range = maximum - minimum || 1;
      drawRasterField(baseContext, width, height, dataset, field.values, minimum, range, view.zoom, view.pan, quality);
      if (gridMode !== "off") {
        gridContext.beginPath(); gridContext.strokeStyle = "rgba(255,255,255,.65)"; gridContext.lineWidth = gridMode === "panel_seams" ? 2 : 0.75;
        const gridPath = geoPath(equirectangularProjection(width, height, view.zoom, view.pan), gridContext);
        for (const edge of gridEdges(dataset.grid.cellsPerPanel)) {
          if (gridMode === "panel_seams" && !edge.panelSeam) continue;
          gridPath({ type: "LineString", coordinates: [unitToLongitudeLatitude(edge.first), unitToLongitudeLatitude(edge.second)] });
        }
        gridContext.stroke();
      }
    };
    const schedule = () => {
      if (frameRequest.current !== 0) return;
      frameRequest.current = requestAnimationFrame(() => { frameRequest.current = 0; draw(); });
    };
    requestDraw.current = schedule;
    const resize = () => {
      const bounds = host.current?.getBoundingClientRect(); if (!bounds) return;
      size.current = [bounds.width, bounds.height];
      for (const canvas of canvases) { const ratio = Math.min(window.devicePixelRatio, 2); canvas.width = Math.max(1, Math.round(bounds.width * ratio)); canvas.height = Math.max(1, Math.round(bounds.height * ratio)); canvas.style.width = `${bounds.width}px`; canvas.style.height = `${bounds.height}px`; canvas.getContext("2d")?.setTransform(ratio, 0, 0, ratio, 0, 0); }
      schedule();
    };
    const observer = new ResizeObserver(resize); observer.observe(host.current); resize();
    return () => { observer.disconnect(); if (frameRequest.current !== 0) cancelAnimationFrame(frameRequest.current); frameRequest.current = 0; requestDraw.current = () => {}; };
  }, []);
  useEffect(() => { requestDraw.current(); }, [dataset, fieldId, gridMode, view]);
  useEffect(() => {
    const element = host.current; if (!element) return undefined;
    let dragging = false; let moved = false; let last: [number, number] = [0, 0]; let settle = 0;
    // Interaction repaints the raster at half resolution and restores full detail once the
    // gesture settles, so dragging stays responsive at N=48/96.
    const interactive = () => {
      scene.current.quality = 0.5;
      window.clearTimeout(settle);
      settle = window.setTimeout(() => { scene.current.quality = 1; requestDraw.current(); }, 160);
    };
    const down = (event: PointerEvent) => { dragging = true; moved = false; last = [event.clientX, event.clientY]; element.setPointerCapture(event.pointerId); };
    const move = (event: PointerEvent) => { if (!dragging) return; const dx = event.clientX - last[0]; const dy = event.clientY - last[1]; if (Math.abs(dx) + Math.abs(dy) > 2) moved = true; interactive(); setView((current) => ({ ...current, pan: [current.pan[0] + dx, current.pan[1] + dy] })); last = [event.clientX, event.clientY]; };
    const up = (event: PointerEvent) => {
      if (!moved) {
        const { dataset: current, view: currentView } = scene.current;
        const bounds = element.getBoundingClientRect(); const x = event.clientX - bounds.left; const y = event.clientY - bounds.top;
        if (isPointInsideProjectedMap(x, y, bounds.width, bounds.height, currentView.zoom, currentView.pan)) {
          const origin = inverseProject(x, y, bounds.width, bounds.height, currentView.zoom, currentView.pan);
          const picked = hitTest(origin, current.grid.cellsPerPanel); const panel = panelNames.indexOf(picked.panel);
          onPickRef.current?.(panel * current.grid.cellsPerPanel ** 2 + picked.j * current.grid.cellsPerPanel + picked.i, origin);
        }
      }
      if (element.hasPointerCapture(event.pointerId)) element.releasePointerCapture(event.pointerId);
      dragging = false;
    };
    const wheel = (event: WheelEvent) => { interactive(); setView((current) => ({ ...current, zoom: Math.max(1, Math.min(8, current.zoom * (event.deltaY < 0 ? 1.1 : 0.9))) })); };
    const cancel = () => { dragging = false; moved = false; };
    element.addEventListener("pointerdown", down); element.addEventListener("pointermove", move); element.addEventListener("pointerup", up); element.addEventListener("pointercancel", cancel); element.addEventListener("wheel", wheel, { passive: true });
    return () => { window.clearTimeout(settle); element.removeEventListener("pointerdown", down); element.removeEventListener("pointermove", move); element.removeEventListener("pointerup", up); element.removeEventListener("pointercancel", cancel); element.removeEventListener("wheel", wheel); };
  }, []);
  return <div ref={host} className="map-2d" aria-label="2D global map"><canvas ref={base} /><canvas ref={grid} /><button type="button" onClick={() => setView({ zoom: 1, pan: [0, 0] })}>Reset map</button></div>;
}

let rasterScratch: HTMLCanvasElement | null = null;

function drawRasterField(context: CanvasRenderingContext2D, width: number, height: number, dataset: VisualDatasetV1, values: Float64Array, minimum: number, range: number, zoom: number, pan: readonly [number, number], quality = 1) {
  const rasterWidth = Math.max(1, Math.ceil(width * quality)); const rasterHeight = Math.max(1, Math.ceil(height * quality));
  const scratch = rasterScratch ?? (rasterScratch = document.createElement("canvas"));
  if (scratch.width !== rasterWidth || scratch.height !== rasterHeight) { scratch.width = rasterWidth; scratch.height = rasterHeight; }
  const scratchContext = scratch.getContext("2d"); if (!scratchContext) return;
  scratchContext.clearRect(0, 0, rasterWidth, rasterHeight);
  const image = scratchContext.createImageData(rasterWidth, rasterHeight);
  const pixels = new Uint32Array(image.data.buffer);
  const littleEndian = new Uint8Array(new Uint32Array([1]).buffer)[0] === 1;
  const colors = new Uint32Array(values.length);
  for (let cell = 0; cell < values.length; cell += 1) {
    const [red, green, blue] = hslToRgb(237 - 237 * (values[cell] - minimum) / range, 0.8, 0.5);
    colors[cell] = littleEndian ? (255 << 24) | (blue << 16) | (green << 8) | red
      : (red << 24) | (green << 16) | (blue << 8) | 255;
  }
  // Pixel coordinates are inverted analytically instead of through d3/hitTest object
  // allocations: this loop runs once per visible pixel on every repaint.
  const scale = Math.min(width / (2 * Math.PI), height / Math.PI) * zoom * quality;
  const centerX = (width / 2 + pan[0]) * quality; const centerY = (height / 2 + pan[1]) * quality;
  const firstX = Math.max(0, Math.ceil(centerX - Math.PI * scale - 0.5)); const lastX = Math.min(rasterWidth - 1, Math.floor(centerX + Math.PI * scale - 0.5));
  const firstY = Math.max(0, Math.ceil(centerY - Math.PI * scale / 2 - 0.5)); const lastY = Math.min(rasterHeight - 1, Math.floor(centerY + Math.PI * scale / 2 - 0.5));
  const cellsPerPanel = dataset.grid.cellsPerPanel;
  for (let y = firstY; y <= lastY; y += 1) {
    const latitude = (centerY - y - 0.5) / scale;
    const cosLatitude = Math.cos(latitude); const z = Math.sin(latitude);
    for (let x = firstX; x <= lastX; x += 1) {
      const longitude = (x + 0.5 - centerX) / scale;
      pixels[y * rasterWidth + x] =
        colors[unitToCellIndex(cosLatitude * Math.cos(longitude), cosLatitude * Math.sin(longitude), z, cellsPerPanel)];
    }
  }
  scratchContext.putImageData(image, 0, 0);
  context.save(); context.imageSmoothingEnabled = false; context.drawImage(scratch, 0, 0, width, height); context.restore();
}

function hslToRgb(hueDegrees: number, saturation: number, lightness: number): readonly [number, number, number] {
  const hue = ((hueDegrees % 360) + 360) % 360 / 60;
  const chroma = (1 - Math.abs(2 * lightness - 1)) * saturation;
  const second = chroma * (1 - Math.abs(hue % 2 - 1));
  const [red, green, blue] = hue < 1 ? [chroma, second, 0] : hue < 2 ? [second, chroma, 0] : hue < 3 ? [0, chroma, second] : hue < 4 ? [0, second, chroma] : hue < 5 ? [second, 0, chroma] : [chroma, 0, second];
  const match = lightness - chroma / 2;
  return [Math.round((red + match) * 255), Math.round((green + match) * 255), Math.round((blue + match) * 255)];
}
