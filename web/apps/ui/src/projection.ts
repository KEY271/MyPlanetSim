import { geoEquirectangular } from "d3-geo";
import { hitTest, UnitVector } from "@myplanetsim/protocol";

export function equirectangularProjection(width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]) {
  return geoEquirectangular().scale(width * zoom / (2 * Math.PI)).translate([width / 2 + pan[0], height / 2 + pan[1]]);
}

export function projectUnit(value: UnitVector, width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]): [number, number] {
  const longitude = Math.atan2(value[1], value[0]) * 180 / Math.PI;
  const latitude = Math.asin(Math.max(-1, Math.min(1, value[2]))) * 180 / Math.PI;
  const projected = equirectangularProjection(width, height, zoom, pan)([longitude, latitude]);
  if (!projected) throw new Error("projection failed");
  return [projected[0], projected[1]];
}

export function inverseProject(x: number, y: number, width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]): UnitVector {
  const projection = equirectangularProjection(width, height, zoom, pan);
  if (!projection.invert) throw new Error("inverse projection is unavailable");
  const inverse = projection.invert([x, y]);
  if (!inverse) throw new Error("inverse projection failed");
  const longitude = inverse[0] * Math.PI / 180; const latitude = inverse[1] * Math.PI / 180;
  return [Math.cos(latitude) * Math.cos(longitude), Math.cos(latitude) * Math.sin(longitude), Math.sin(latitude)];
}

export function pickMapCell(x: number, y: number, width: number, height: number, cellsPerPanel: number, zoom = 1, pan: readonly [number, number] = [0, 0]) {
  return hitTest(inverseProject(x, y, width, height, zoom, pan), cellsPerPanel);
}

export function drawWrappedSegment(context: CanvasRenderingContext2D, first: UnitVector, second: UnitVector, width: number, height: number, zoom = 1, pan: readonly [number, number] = [0, 0]) {
  const firstPoint = projectUnit(first, width, height, zoom, pan); const secondPoint = projectUnit(second, width, height, zoom, pan);
  if (Math.abs(firstPoint[0] - secondPoint[0]) <= width / 2) { context.moveTo(firstPoint[0], firstPoint[1]); context.lineTo(secondPoint[0], secondPoint[1]); return; }
  const adjusted = secondPoint[0] < firstPoint[0] ? secondPoint[0] + width : secondPoint[0] - width;
  context.moveTo(firstPoint[0], firstPoint[1]); context.lineTo(adjusted, secondPoint[1]);
}
