import * as THREE from "three";
import { cubedSphereCells, gridEdges, VisualDatasetV1 } from "@myplanetsim/protocol";

export function createGlobeGeometry(dataset: VisualDatasetV1, fieldId: string): THREE.BufferGeometry {
  const field = dataset.fields.find((candidate) => candidate.id === fieldId);
  if (!field) throw new Error(`unknown field: ${fieldId}`);
  const cells = cubedSphereCells(dataset.grid.cellsPerPanel);
  const positions = new Float32Array(cells.length * 4 * 3);
  const colors = new Float32Array(cells.length * 4 * 3);
  const cellIds = new Uint32Array(cells.length * 4);
  const minimum = Math.min(...field.values);
  const maximum = Math.max(...field.values);
  const range = maximum - minimum || 1;
  for (let cell = 0; cell < cells.length; cell += 1) {
    const color = new THREE.Color().setHSL(0.66 - 0.66 * (field.values[cell] - minimum) / range, 0.8, 0.5);
    for (let corner = 0; corner < 4; corner += 1) {
      const positionOffset = (cell * 4 + corner) * 3;
      positions.set(cells[cell].corners[corner], positionOffset);
      colors[positionOffset] = color.r;
      colors[positionOffset + 1] = color.g;
      colors[positionOffset + 2] = color.b;
      cellIds[cell * 4 + corner] = cell;
    }
  }
  const indices = new Uint32Array(cells.length * 6);
  for (let cell = 0; cell < cells.length; cell += 1) {
    const vertex = cell * 4;
    indices.set([vertex, vertex + 1, vertex + 2, vertex, vertex + 2, vertex + 3], cell * 6);
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
    positions.push(...edge.first, ...edge.second);
  }
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
  return geometry;
}

export function pickedCellFromFace(faceIndex: number | null): number | null {
  return faceIndex === null ? null : Math.floor(faceIndex / 2);
}
