import { FrameV2Column } from "@myplanetsim/protocol";

export const profileSize = { width: 300, height: 320 } as const;
export const profilePlot = { left: 62, right: 286, top: 20, bottom: 272 } as const;

export interface ProfilePoint {
  readonly index: number;
  readonly x: number;
  readonly y: number;
}

export interface ProfileTick {
  readonly index: number;
  readonly pressurePa: number;
  readonly y: number;
}

export interface ProfileGeometry {
  readonly points: readonly ProfilePoint[];
  readonly ticks: readonly ProfileTick[];
  readonly minimum: number;
  readonly maximum: number;
}

// Pressure uses a downward-increasing logarithmic axis, which is the standard reading
// order for an atmospheric column: the model top sits at the top of the chart. The axis is
// spanned by the column's own decoded pressures, so no hydrostatic relation is recomputed
// in the browser.
export function pressureToY(pressurePa: number, topPa: number, bottomPa: number): number {
  const span = Math.log(bottomPa) - Math.log(topPa);
  const fraction = span > 0 ? (Math.log(pressurePa) - Math.log(topPa)) / span : 0.5;
  return profilePlot.top + fraction * (profilePlot.bottom - profilePlot.top);
}

export function profileGeometry(column: FrameV2Column, maximumTicks = 6): ProfileGeometry {
  const levels = column.values.length;
  if (levels === 0) throw new Error("a column profile needs at least one level");
  const minimum = column.values.reduce((low, value) => Math.min(low, value), Infinity);
  const maximum = column.values.reduce((high, value) => Math.max(high, value), -Infinity);
  // A constant column has no horizontal extent, so it is centred instead of divided by
  // zero; the axis labels still report the single value.
  const span = maximum - minimum;
  const topPa = column.pressurePa[0];
  const bottomPa = column.pressurePa[levels - 1];
  const points = Array.from(column.values, (value, index) => ({
    index,
    x: span > 0
      ? profilePlot.left + ((value - minimum) / span) * (profilePlot.right - profilePlot.left)
      : (profilePlot.left + profilePlot.right) / 2,
    y: pressureToY(column.pressurePa[index], topPa, bottomPa),
  }));
  const stride = Math.max(1, Math.ceil(levels / maximumTicks));
  const ticks = points
    .filter((point) => point.index % stride === 0 || point.index === levels - 1)
    .map((point) => ({ index: point.index, pressurePa: column.pressurePa[point.index], y: point.y }));
  return { points, ticks, minimum, maximum };
}
