import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { cubedSphereCells, decodeFrame, frameV2Field, frameV2Fields, frameV2LevelDataset, frameV2Sample, gridEdges, hitTest, levelSlice, selectedColumn, VisualFrameV2 } from "@myplanetsim/protocol";
import { profileGeometry, profilePlot } from "./profile";
import { formatValue } from "./format";

const cellsPerPanel = 4;
const levels = 8;
const cellCount = 6 * cellsPerPanel * cellsPerPanel;

// A stand-in for one published FrameV2: pressure grows downward through the column, and
// every other field is a distinct function of (cell, level) so a transposed index, a
// dropped level, or a shared buffer is visible in the assertions.
function fixtureFrame(): VisualFrameV2 {
  const volume = cellCount * levels;
  const fingerprint = new TextEncoder().encode("fixture");
  const bytes = new ArrayBuffer(60 + fingerprint.length + 8 * (cellCount + 7 * volume));
  const view = new DataView(bytes);
  new Uint8Array(bytes, 0, 8).set(new TextEncoder().encode("MPSFRAM2"));
  view.setUint32(8, 2, true);
  view.setBigUint64(16, BigInt(cellsPerPanel), true);
  view.setBigUint64(24, BigInt(levels), true);
  view.setFloat64(32, 300, true);
  view.setBigUint64(40, 30n, true);
  view.setBigUint64(48, BigInt(cellCount), true);
  view.setUint32(56, fingerprint.length, true);
  new Uint8Array(bytes, 60, fingerprint.length).set(fingerprint);
  let offset = 60 + fingerprint.length;
  const write = (value: number) => { view.setFloat64(offset, value, true); offset += 8; };
  for (let cell = 0; cell < cellCount; cell += 1) write(100000 + cell);
  const volumeField = (produce: (cell: number, level: number) => void) => {
    for (let cell = 0; cell < cellCount; cell += 1) for (let level = 0; level < levels; level += 1) produce(cell, level);
  };
  volumeField((cell, level) => write(2000 + level * 12000 + cell));       // pressure
  volumeField((cell, level) => write(300 + level + cell / 1000));         // theta
  volumeField((cell, level) => write(230 + 5 * level + cell / 1000));     // temperature
  volumeField((cell, level) => write(level - 3));                         // wind x, signed
  volumeField(() => write(0));                                           // wind y
  volumeField(() => write(0));                                           // wind z
  volumeField(() => write(0.5));                                         // tracer, constant
  const frame = decodeFrame(bytes);
  if (frame.schemaVersion !== 2) throw new Error("fixture is not a FrameV2");
  return frame;
}

describe("dry hydrostatic visualization", () => {
  it("shows the same C++ value in the map, the globe, the profile, and the inspector", () => {
    const frame = fixtureFrame();
    // The 2D map and the 3D globe are handed the identical adapted dataset, so agreement
    // between the views reduces to agreement between the dataset, the column, and the sample.
    for (const field of frameV2Fields.filter((candidate) => candidate.volume)) {
      for (const level of [0, 3, levels - 1]) {
        const values = frameV2LevelDataset(frame, field.id, level).fields[0].values;
        expect(values).toHaveLength(cellCount);
        for (let cell = 0; cell < cellCount; cell += 1) {
          const column = selectedColumn(frame, field.id as "temperature", cell);
          expect(values[cell]).toBe(column.values[level]);
          expect(values[cell]).toBe(frameV2Sample(frame, field.id, cell, level).value);
          expect(column.pressurePa[level]).toBe(frameV2Sample(frame, field.id, cell, level).pressurePa);
        }
      }
    }
  });

  it("carries the frame time, step, and fingerprint into every adapted level", () => {
    const frame = fixtureFrame();
    for (let level = 0; level < levels; level += 1) {
      expect(frameV2LevelDataset(frame, "temperature", level).frame)
        .toEqual({ timeSeconds: 300, step: 30, configFingerprint: "fixture" });
    }
  });

  it("renders the top level, the bottom level, and seam and corner cells", () => {
    const frame = fixtureFrame();
    const top = frameV2LevelDataset(frame, "temperature", 0).fields[0].values;
    const bottom = frameV2LevelDataset(frame, "temperature", levels - 1).fields[0].values;
    expect(top[0]).toBeCloseTo(230, 9);
    expect(bottom[0]).toBeCloseTo(230 + 5 * (levels - 1), 9);
    // A panel corner and the first cell across a seam must be ordinary members of the
    // slice, not holes left by the cell-major flatten order.
    const corner = 0;
    const seam = cellsPerPanel * cellsPerPanel;
    for (const cell of [corner, seam, cellCount - 1]) {
      expect(Number.isFinite(top[cell])).toBe(true);
      expect(top[cell]).toBe(levelSlice(frame, "temperature", 0)[cell]);
      const geometry = cubedSphereCells(cellsPerPanel)[cell];
      expect(hitTest(geometry.center, cellsPerPanel)).toEqual({ panel: geometry.panel, i: geometry.i, j: geometry.j });
    }
  });

  it("displays a surface field, a constant field, and a signed field", () => {
    const frame = fixtureFrame();
    const surface = frameV2LevelDataset(frame, "surface_pressure", 0).fields[0];
    expect(surface.unit).toBe("Pa");
    expect(surface.values[5]).toBe(100005);
    expect(frameV2Field("surface_pressure").volume).toBe(false);
    const tracer = frameV2LevelDataset(frame, "tracer", 4).fields[0].values;
    expect([...new Set(tracer)]).toEqual([0.5]);
    const wind = frameV2LevelDataset(frame, "wind_speed", 0).fields[0].values;
    // Wind speed is a magnitude, so a negative component still displays as positive.
    expect(wind[0]).toBe(3);
    expect(frameV2Sample(frame, "wind_speed", 0, 3).value).toBe(0);
  });

  it("keeps the memoized geometry alive across field, level, and frame changes", () => {
    const frame = fixtureFrame();
    const cells = cubedSphereCells(cellsPerPanel);
    const edges = gridEdges(cellsPerPanel);
    for (const field of frameV2Fields) {
      for (let level = 0; level < levels; level += 1) {
        const value = frameV2LevelDataset(frame, field.id, field.volume ? level : 0);
        expect(value.grid.cellsPerPanel).toBe(cellsPerPanel);
        expect(cubedSphereCells(value.grid.cellsPerPanel)).toBe(cells);
        expect(gridEdges(value.grid.cellsPerPanel)).toBe(edges);
      }
    }
  });

  it("plots pressure downward on a logarithmic axis", () => {
    const frame = fixtureFrame();
    const geometry = profileGeometry(selectedColumn(frame, "temperature", 9));
    expect(geometry.points).toHaveLength(levels);
    expect(geometry.points[0].y).toBeCloseTo(profilePlot.top, 9);
    expect(geometry.points[levels - 1].y).toBeCloseTo(profilePlot.bottom, 9);
    for (let level = 1; level < levels; level += 1) {
      expect(geometry.points[level].y).toBeGreaterThan(geometry.points[level - 1].y);
    }
    const column = selectedColumn(frame, "temperature", 9);
    const top = Math.log(column.pressurePa[0]);
    const span = Math.log(column.pressurePa[levels - 1]) - top;
    // A logarithmic axis places the mid-pressure level away from the geometric midpoint a
    // linear axis would use.
    const expected = profilePlot.top + ((Math.log(column.pressurePa[4]) - top) / span) * (profilePlot.bottom - profilePlot.top);
    expect(geometry.points[4].y).toBeCloseTo(expected, 9);
    expect(geometry.points[4].y).not.toBeCloseTo((profilePlot.top + profilePlot.bottom) / 2, 3);
    expect(geometry.ticks.map((tick) => tick.index).at(-1)).toBe(levels - 1);
  });

  it("centres a constant column instead of dividing by a zero value span", () => {
    const frame = fixtureFrame();
    const geometry = profileGeometry(selectedColumn(frame, "tracer", 2));
    expect(geometry.minimum).toBe(0.5);
    expect(geometry.maximum).toBe(0.5);
    for (const point of geometry.points) {
      expect(point.x).toBeCloseTo((profilePlot.left + profilePlot.right) / 2, 9);
      expect(Number.isFinite(point.y)).toBe(true);
    }
  });

  it("formats raw solver values without hiding their magnitude", () => {
    expect(formatValue(100000)).toBe("100000");
    expect(formatValue(288.15000001)).toBe("288.15");
    expect(formatValue(-0.0000123)).toBe("-1.2300e-5");
    expect(formatValue(0)).toBe("0");
  });

  it("reserves a third view column for the profile and keeps the reduced-motion gate", () => {
    const style = readFileSync(new URL("./style.css", import.meta.url), "utf8");
    expect(style).toContain(".view-grid.with-profile");
    expect(style).toContain("prefers-reduced-motion");
    // A narrow viewport collapses all three views into one scrolling column.
    expect(style).toContain("grid-template-columns: minmax(0, 1fr); grid-auto-rows: minmax(0, 1fr)");
  });
});
