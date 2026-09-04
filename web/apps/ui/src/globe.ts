import * as THREE from "three";
import { cubedSphereCells, fieldRange, gridEdges, mapToUnitSphere, UnitVector, VisualDatasetV1, VisualFieldV1 } from "@myplanetsim/protocol";

const minimumPanelSegments = 24;

export function globeSubdivisionsPerCell(cellsPerPanel: number): number {
  return Math.max(1, Math.ceil(minimumPanelSegments / cellsPerPanel));
}

export function globeTrianglesPerCell(cellsPerPanel: number): number {
  const subdivisions = globeSubdivisionsPerCell(cellsPerPanel);
  return 2 * subdivisions * subdivisions;
}

function requireField(dataset: VisualDatasetV1, fieldId: string): VisualFieldV1 {
  const field = dataset.fields.find((candidate) => candidate.id === fieldId);
  if (!field) throw new Error(`unknown field: ${fieldId}`);
  return field;
}

function writeCellColors(colors: Float32Array, field: VisualFieldV1, verticesPerCell: number): void {
  const [minimum, maximum] = fieldRange(field);
  const range = maximum - minimum || 1;
  const color = new THREE.Color();
  for (let cell = 0; cell < field.values.length; cell += 1) {
    color.setHSL(0.66 - 0.66 * (field.values[cell] - minimum) / range, 0.8, 0.5);
    for (let vertex = 0; vertex < verticesPerCell; vertex += 1) {
      const offset = (cell * verticesPerCell + vertex) * 3;
      colors[offset] = color.r;
      colors[offset + 1] = color.g;
      colors[offset + 2] = color.b;
    }
  }
}

// Rewrites only the colour attribute when the grid resolution is unchanged, so playback and
// field switches do not rebuild positions/indices for every frame.
export function updateGlobeColors(geometry: THREE.BufferGeometry, dataset: VisualDatasetV1, fieldId: string): boolean {
  const attribute = geometry.getAttribute("color") as THREE.BufferAttribute | undefined;
  const verticesPerCell = (globeSubdivisionsPerCell(dataset.grid.cellsPerPanel) + 1) ** 2;
  const expected = 6 * dataset.grid.cellsPerPanel ** 2 * verticesPerCell;
  if (!attribute || attribute.count !== expected || !(attribute.array instanceof Float32Array)) return false;
  writeCellColors(attribute.array, requireField(dataset, fieldId), verticesPerCell);
  attribute.needsUpdate = true;
  return true;
}

export function createGlobeGeometry(dataset: VisualDatasetV1, fieldId: string): THREE.BufferGeometry {
  const field = requireField(dataset, fieldId);
  const cells = cubedSphereCells(dataset.grid.cellsPerPanel);
  const subdivisions = globeSubdivisionsPerCell(dataset.grid.cellsPerPanel);
  const verticesPerSide = subdivisions + 1;
  const verticesPerCell = verticesPerSide * verticesPerSide;
  const positions = new Float32Array(cells.length * verticesPerCell * 3);
  const colors = new Float32Array(cells.length * verticesPerCell * 3);
  const cellIds = new Uint32Array(cells.length * verticesPerCell);
  const indices = new Uint32Array(cells.length * subdivisions * subdivisions * 6);
  const panelDelta = Math.PI / (2 * dataset.grid.cellsPerPanel);
  writeCellColors(colors, field, verticesPerCell);
  for (let cell = 0; cell < cells.length; cell += 1) {
    const geometry = cells[cell];
    for (let row = 0; row <= subdivisions; row += 1) for (let column = 0; column <= subdivisions; column += 1) {
      const vertex = cell * verticesPerCell + row * verticesPerSide + column;
      const positionOffset = vertex * 3;
      const alpha = -Math.PI / 4 + (geometry.i + column / subdivisions) * panelDelta;
      const beta = -Math.PI / 4 + (geometry.j + row / subdivisions) * panelDelta;
      positions.set(mapToUnitSphere(geometry.panel, alpha, beta), positionOffset);
      cellIds[vertex] = cell;
    }
    for (let row = 0; row < subdivisions; row += 1) for (let column = 0; column < subdivisions; column += 1) {
      const topLeft = cell * verticesPerCell + row * verticesPerSide + column;
      const topRight = topLeft + 1;
      const bottomLeft = topLeft + verticesPerSide;
      const bottomRight = bottomLeft + 1;
      const indexOffset = (cell * subdivisions * subdivisions + row * subdivisions + column) * 6;
      indices.set([topLeft, topRight, bottomRight, topLeft, bottomRight, bottomLeft], indexOffset);
    }
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geometry.setAttribute("cellId", new THREE.Uint32BufferAttribute(cellIds, 1));
  geometry.setIndex(new THREE.BufferAttribute(indices, 1));
  return geometry;
}

export function createGridGeometry(cellsPerPanel: number, mode: "off" | "panel_seams" | "all_cells"): THREE.BufferGeometry {
  const geometry = new THREE.BufferGeometry();
  if (mode === "off") return geometry;
  const positions: number[] = [];
  for (const edge of gridEdges(cellsPerPanel)) {
    if (mode === "panel_seams" && !edge.panelSeam) continue;
    const angle = Math.acos(THREE.MathUtils.clamp(edge.first[0] * edge.second[0] + edge.first[1] * edge.second[1] + edge.first[2] * edge.second[2], -1, 1));
    const subdivisions = Math.max(1, Math.ceil(angle / (Math.PI / 60)));
    for (let segment = 0; segment < subdivisions; segment += 1) {
      positions.push(...sphericalInterpolate(edge.first, edge.second, segment / subdivisions),
        ...sphericalInterpolate(edge.first, edge.second, (segment + 1) / subdivisions));
    }
  }
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
  return geometry;
}

function sphericalInterpolate(first: UnitVector, second: UnitVector, amount: number): UnitVector {
  const angle = Math.acos(THREE.MathUtils.clamp(first[0] * second[0] + first[1] * second[1] + first[2] * second[2], -1, 1));
  if (angle < 1e-12) return first;
  const denominator = Math.sin(angle);
  const firstWeight = Math.sin((1 - amount) * angle) / denominator;
  const secondWeight = Math.sin(amount * angle) / denominator;
  return [first[0] * firstWeight + second[0] * secondWeight,
    first[1] * firstWeight + second[1] * secondWeight,
    first[2] * firstWeight + second[2] * secondWeight];
}

export function pickedCellFromFace(faceIndex: number | null, trianglesPerCell = 2): number | null {
  return faceIndex === null ? null : Math.floor(faceIndex / trianglesPerCell);
}
