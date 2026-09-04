import { geoArea, geoEquirectangular } from "d3-geo";
import { hitTest, UnitVector } from "@myplanetsim/protocol";

export function equirectangularProjection(width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]) {
  const fittedScale = Math.min(width / (2 * Math.PI), height / Math.PI);
  return geoEquirectangular().scale(fittedScale * zoom).translate([width / 2 + pan[0], height / 2 + pan[1]]);
}

export function isPointInsideProjectedMap(x: number, y: number, width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]): boolean {
  const scale = Math.min(width / (2 * Math.PI), height / Math.PI) * zoom;
  const centerX = width / 2 + pan[0]; const centerY = height / 2 + pan[1];
  return x >= centerX - Math.PI * scale && x <= centerX + Math.PI * scale &&
    y >= centerY - Math.PI * scale / 2 && y <= centerY + Math.PI * scale / 2;
}

export function projectUnit(value: UnitVector, width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]): [number, number] {
  const [longitude, latitude] = unitToLongitudeLatitude(value);
  const projected = equirectangularProjection(width, height, zoom, pan)([longitude, latitude]);
  if (!projected) throw new Error("projection failed");
  return [projected[0], projected[1]];
}

export function unitToLongitudeLatitude(value: UnitVector): [number, number] {
  return [Math.atan2(value[1], value[0]) * 180 / Math.PI,
    Math.asin(Math.max(-1, Math.min(1, value[2]))) * 180 / Math.PI];
}

export function geoCellPolygon(corners: readonly UnitVector[]) {
  let ring = corners.map(unitToLongitudeLatitude);
  let polygon = { type: "Polygon" as const, coordinates: [[...ring, ring[0]]] };
  if (geoArea(polygon) > 2 * Math.PI) {
    ring = [...ring].reverse();
    polygon = { type: "Polygon", coordinates: [[...ring, ring[0]]] };
  }
  return polygon;
}

export function inverseProject(x: number, y: number, width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]): UnitVector {
  const scale = Math.min(width / (2 * Math.PI), height / Math.PI) * zoom;
  const longitude = (x - width / 2 - pan[0]) / scale;
  const latitude = (height / 2 + pan[1] - y) / scale;
  return [Math.cos(latitude) * Math.cos(longitude), Math.cos(latitude) * Math.sin(longitude), Math.sin(latitude)];
}

export function pickMapCell(x: number, y: number, width: number, height: number, cellsPerPanel: number, zoom = 1, pan: readonly [number, number] = [0, 0]) {
  return hitTest(inverseProject(x, y, width, height, zoom, pan), cellsPerPanel);
}

// Allocation-free equivalent of hitTest plus the panel-major flat index, used by the raster
// field renderer where it runs once per visible pixel on every repaint.
export function unitToCellIndex(x: number, y: number, z: number, cellsPerPanel: number): number {
  const absX = Math.abs(x); const absY = Math.abs(y); const absZ = Math.abs(z);
  let panel: number; let alpha: number; let beta: number;
  if (absX >= absY && absX >= absZ) {
    panel = x >= 0 ? 0 : 2;
    alpha = Math.atan2(x >= 0 ? y : -y, absX); beta = Math.atan2(z, absX);
  } else if (absY >= absZ) {
    panel = y >= 0 ? 1 : 3;
    alpha = Math.atan2(y >= 0 ? -x : x, absY); beta = Math.atan2(z, absY);
  } else {
    panel = z >= 0 ? 4 : 5;
    alpha = Math.atan2(x, absZ); beta = Math.atan2(z >= 0 ? y : -y, absZ);
  }
  const scale = cellsPerPanel / (Math.PI / 2); const quarter = Math.PI / 4;
  const column = Math.max(0, Math.min(cellsPerPanel - 1, Math.floor((alpha + quarter) * scale)));
  const row = Math.max(0, Math.min(cellsPerPanel - 1, Math.floor((beta + quarter) * scale)));
  return (panel * cellsPerPanel + row) * cellsPerPanel + column;
}
