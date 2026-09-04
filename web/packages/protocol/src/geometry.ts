import { ProtocolError } from "./index";

export const panelNames = ["PX", "PY", "NX", "NY", "PZ", "NZ"] as const;
export type PanelName = (typeof panelNames)[number];
export type UnitVector = readonly [number, number, number];
export type GridMode = "off" | "panel_seams" | "all_cells";

interface Basis { readonly normal: UnitVector; readonly alpha: UnitVector; readonly beta: UnitVector }
const bases: Record<PanelName, Basis> = {
  PX: { normal: [1, 0, 0], alpha: [0, 1, 0], beta: [0, 0, 1] },
  PY: { normal: [0, 1, 0], alpha: [-1, 0, 0], beta: [0, 0, 1] },
  NX: { normal: [-1, 0, 0], alpha: [0, -1, 0], beta: [0, 0, 1] },
  NY: { normal: [0, -1, 0], alpha: [1, 0, 0], beta: [0, 0, 1] },
  PZ: { normal: [0, 0, 1], alpha: [1, 0, 0], beta: [0, 1, 0] },
  NZ: { normal: [0, 0, -1], alpha: [1, 0, 0], beta: [0, -1, 0] },
};
const pi = Math.PI;

function unit(value: readonly [number, number, number]): UnitVector {
  const length = Math.hypot(...value);
  if (!Number.isFinite(length) || length === 0) throw new ProtocolError("invalid_geometry", "zero vector");
  return [value[0] / length, value[1] / length, value[2] / length];
}

function dot(a: readonly number[], b: readonly number[]): number { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

export function mapToUnitSphere(panel: PanelName, alpha: number, beta: number): UnitVector {
  if (!Number.isFinite(alpha) || !Number.isFinite(beta)) throw new ProtocolError("invalid_geometry", "panel coordinates are not finite");
  const basis = bases[panel];
  return unit([basis.normal[0] + Math.tan(alpha) * basis.alpha[0] + Math.tan(beta) * basis.beta[0],
    basis.normal[1] + Math.tan(alpha) * basis.alpha[1] + Math.tan(beta) * basis.beta[1],
    basis.normal[2] + Math.tan(alpha) * basis.alpha[2] + Math.tan(beta) * basis.beta[2]]);
}

export function inverseMap(position: UnitVector): { panel: PanelName; alpha: number; beta: number } {
  const point = unit(position);
  const axis = [Math.abs(point[0]), Math.abs(point[1]), Math.abs(point[2])];
  const selected = axis[0] >= axis[1] && axis[0] >= axis[2] ? 0 : axis[1] >= axis[2] ? 1 : 2;
  const panel: PanelName = selected === 0 ? point[0] >= 0 ? "PX" : "NX" :
    selected === 1 ? point[1] >= 0 ? "PY" : "NY" : point[2] >= 0 ? "PZ" : "NZ";
  const basis = bases[panel];
  const denominator = dot(point, basis.normal);
  return { panel, alpha: Math.atan2(dot(point, basis.alpha), denominator), beta: Math.atan2(dot(point, basis.beta), denominator) };
}

export interface CellGeometryV1 {
  readonly panel: PanelName;
  readonly i: number;
  readonly j: number;
  readonly center: UnitVector;
  readonly corners: readonly [UnitVector, UnitVector, UnitVector, UnitVector];
}

export interface EdgeGeometryV1 {
  readonly id: number;
  readonly first: UnitVector;
  readonly second: UnitVector;
  readonly panelSeam: boolean;
}

export function cellGeometry(panel: PanelName, i: number, j: number, cellsPerPanel: number): CellGeometryV1 {
  if (!Number.isInteger(cellsPerPanel) || cellsPerPanel < 1 || !Number.isInteger(i) || !Number.isInteger(j) ||
      i < 0 || i >= cellsPerPanel || j < 0 || j >= cellsPerPanel) {
    throw new ProtocolError("invalid_geometry", "cell index is invalid");
  }
  const delta = pi / (2 * cellsPerPanel);
  const coordinate = (index: number) => -pi / 4 + index * delta;
  const corners = [mapToUnitSphere(panel, coordinate(i), coordinate(j)),
    mapToUnitSphere(panel, coordinate(i + 1), coordinate(j)),
    mapToUnitSphere(panel, coordinate(i + 1), coordinate(j + 1)),
    mapToUnitSphere(panel, coordinate(i), coordinate(j + 1))] as const;
  return { panel, i, j, center: mapToUnitSphere(panel, coordinate(i + 0.5), coordinate(j + 0.5)), corners };
}

function vertexKey(value: UnitVector): string {
  return value.map((component) => Math.round(component * 1e12)).join(":");
}

// Panel geometry depends only on cellsPerPanel, so the two derived structures are cached.
// Rebuilding them per animation frame or per map pan dominates rendering cost at N=48/96.
function memoizeByCellsPerPanel<T>(build: (cellsPerPanel: number) => T): (cellsPerPanel: number) => T {
  const cache = new Map<number, T>();
  return (cellsPerPanel: number) => {
    const cached = cache.get(cellsPerPanel);
    if (cached !== undefined) return cached;
    const value = build(cellsPerPanel);
    if (cache.size >= 2) cache.delete(cache.keys().next().value as number);
    cache.set(cellsPerPanel, value);
    return value;
  };
}

export const cubedSphereCells = memoizeByCellsPerPanel((cellsPerPanel: number): readonly CellGeometryV1[] =>
  panelNames.flatMap((panel) => Array.from({ length: cellsPerPanel * cellsPerPanel }, (_, index) =>
    cellGeometry(panel, index % cellsPerPanel, Math.floor(index / cellsPerPanel), cellsPerPanel))));

export const gridEdges = memoizeByCellsPerPanel(buildGridEdges);

function buildGridEdges(cellsPerPanel: number): readonly EdgeGeometryV1[] {
  const cells = cubedSphereCells(cellsPerPanel);
  const edges = new Map<string, { first: UnitVector; second: UnitVector; panels: Set<PanelName> }>();
  for (const cell of cells) for (let side = 0; side < 4; side += 1) {
    const first = cell.corners[side]; const second = cell.corners[(side + 1) % 4];
    const key = [vertexKey(first), vertexKey(second)].sort().join("|");
    const current = edges.get(key);
    if (current) current.panels.add(cell.panel);
    else edges.set(key, { first, second, panels: new Set([cell.panel]) });
  }
  return [...edges.values()].map((edge, id) => ({ id, first: edge.first, second: edge.second,
    panelSeam: edge.panels.size > 1 }));
}

export function hitTest(position: UnitVector, cellsPerPanel: number): { panel: PanelName; i: number; j: number } {
  const mapped = inverseMap(position);
  const scale = cellsPerPanel / (pi / 2);
  return { panel: mapped.panel,
    i: Math.max(0, Math.min(cellsPerPanel - 1, Math.floor((mapped.alpha + pi / 4) * scale))),
    j: Math.max(0, Math.min(cellsPerPanel - 1, Math.floor((mapped.beta + pi / 4) * scale))) };
}
