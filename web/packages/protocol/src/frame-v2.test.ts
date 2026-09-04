import { describe, expect, it } from "vitest";
import { ProtocolError, validateRunRequest } from "./index";
import { presetDescriptors, shallowWaterPresetDescriptor } from "./client";
import { decodeFrame, decodeFrameV1, decodeFrameV2, frameV2CellCount, frameV2Field, frameV2Fields, FrameV2FieldId, frameV2LevelDataset, frameV2Sample, levelSlice, selectedColumn, VisualFrameV2 } from "./visual";

const fieldCount = 7;

interface FrameOptions {
  readonly cellsPerPanel?: number;
  readonly levels?: number;
  readonly fingerprint?: string;
  readonly magic?: string;
  readonly schemaVersion?: number;
  readonly declaredCells?: number;
  readonly value?: (field: number, cell: number, level: number) => number;
}

// The fixture writes the exact byte layout Phase 5 fixed for FrameV2 so the decoder is
// tested against the contract rather than against itself.
function encodeFrameV2(options: FrameOptions = {}): ArrayBuffer {
  const cellsPerPanel = options.cellsPerPanel ?? 2;
  const levels = options.levels ?? 3;
  const cells = 6 * cellsPerPanel * cellsPerPanel;
  const volume = cells * levels;
  const fingerprint = new TextEncoder().encode(options.fingerprint ?? "fingerprint-1");
  const bytes = new ArrayBuffer(60 + fingerprint.length + 8 * (cells + fieldCount * volume));
  const view = new DataView(bytes);
  new Uint8Array(bytes, 0, 8).set(new TextEncoder().encode(options.magic ?? "MPSFRAM2"));
  view.setUint32(8, options.schemaVersion ?? 2, true);
  view.setUint32(12, 0, true);
  view.setBigUint64(16, BigInt(cellsPerPanel), true);
  view.setBigUint64(24, BigInt(levels), true);
  view.setFloat64(32, 120, true);
  view.setBigUint64(40, 12n, true);
  view.setBigUint64(48, BigInt(options.declaredCells ?? cells), true);
  view.setUint32(56, fingerprint.length, true);
  new Uint8Array(bytes, 60, fingerprint.length).set(fingerprint);
  // Distinct per (field, cell, level) so a swapped offset or a transposed flatten order
  // cannot pass by coincidence.
  const value = options.value ?? ((field, cell, level) => field * 1e6 + cell * 1e3 + level);
  let offset = 60 + fingerprint.length;
  for (let cell = 0; cell < cells; cell += 1) { view.setFloat64(offset, value(0, cell, 0), true); offset += 8; }
  for (let field = 1; field <= fieldCount; field += 1) {
    for (let cell = 0; cell < cells; cell += 1) {
      for (let level = 0; level < levels; level += 1) { view.setFloat64(offset, value(field, cell, level), true); offset += 8; }
    }
  }
  return bytes;
}

function encodeFrameV1(cellsPerPanel = 2, fingerprint = "fingerprint-1"): ArrayBuffer {
  const cells = 6 * cellsPerPanel * cellsPerPanel;
  const label = new TextEncoder().encode(fingerprint);
  const bytes = new ArrayBuffer(52 + label.length + 32 * cells);
  const view = new DataView(bytes);
  new Uint8Array(bytes, 0, 8).set(new TextEncoder().encode("MPSFRAM1"));
  view.setUint32(8, 1, true);
  view.setBigUint64(16, BigInt(cellsPerPanel), true);
  view.setFloat64(24, 60, true);
  view.setBigUint64(32, 6n, true);
  view.setBigUint64(40, BigInt(cells), true);
  view.setUint32(48, label.length, true);
  new Uint8Array(bytes, 52, label.length).set(label);
  for (let index = 0; index < 4 * cells; index += 1) {
    view.setFloat64(52 + label.length + index * 8, index < cells ? 3000 : 1, true);
  }
  return bytes;
}

describe("FrameV2 decoding", () => {
  it("dispatches on magic and keeps the FrameV1 path intact", () => {
    const v1 = decodeFrame(encodeFrameV1());
    const v2 = decodeFrame(encodeFrameV2());
    expect(v1.schemaVersion).toBe(1);
    expect(v2.schemaVersion).toBe(2);
    expect(v1).toEqual(decodeFrameV1(encodeFrameV1()));
    const bad = encodeFrameV2({ magic: "MPSFRAM9" });
    expect(() => decodeFrame(bad)).toThrow(/unsupported frame magic/);
  });

  it("reads every payload field at its own offset in cell-major then level order", () => {
    const frame = decodeFrameV2(encodeFrameV2()) as VisualFrameV2;
    expect(frame.cellsPerPanel).toBe(2);
    expect(frame.levels).toBe(3);
    expect(frame.timeSeconds).toBe(120);
    expect(frame.step).toBe(12);
    expect(frame.configFingerprint).toBe("fingerprint-1");
    expect(frameV2CellCount(frame)).toBe(24);
    expect(frame.surfacePressurePa).toHaveLength(24);
    for (const values of [frame.pressurePa, frame.potentialTemperatureK, frame.temperatureK,
      frame.windXMetersPerSecond, frame.windYMetersPerSecond, frame.windZMetersPerSecond,
      frame.tracerMixingRatio]) expect(values).toHaveLength(72);
    // Payload order is ps, pressure, theta, T, wind x/y/z, tracer, and the volume offset
    // of a cell-major frame is cell*K + k.
    expect(frame.surfacePressurePa[7]).toBe(7000);
    expect(frame.pressurePa[7 * 3 + 2]).toBe(1e6 + 7000 + 2);
    expect(frame.potentialTemperatureK[7 * 3 + 2]).toBe(2e6 + 7000 + 2);
    expect(frame.temperatureK[7 * 3 + 2]).toBe(3e6 + 7000 + 2);
    expect(frame.windXMetersPerSecond[7 * 3 + 2]).toBe(4e6 + 7000 + 2);
    expect(frame.windYMetersPerSecond[7 * 3 + 2]).toBe(5e6 + 7000 + 2);
    expect(frame.windZMetersPerSecond[7 * 3 + 2]).toBe(6e6 + 7000 + 2);
    expect(frame.tracerMixingRatio[7 * 3 + 2]).toBe(7e6 + 7000 + 2);
  });

  it("agrees between a level slice, a selected column, and a single sample", () => {
    const frame = decodeFrameV2(encodeFrameV2({ cellsPerPanel: 3, levels: 5 })) as VisualFrameV2;
    for (const field of ["pressure", "potential_temperature", "temperature", "tracer", "wind_speed"] as const) {
      for (const cell of [0, 1, 26, 53]) {
        const column = selectedColumn(frame, field, cell);
        for (let level = 0; level < frame.levels; level += 1) {
          expect(column.values[level]).toBe(levelSlice(frame, field, level)[cell]);
          expect(frameV2Sample(frame, field, cell, level).value).toBe(column.values[level]);
          expect(frameV2Sample(frame, field, cell, level).pressurePa).toBe(column.pressurePa[level]);
        }
      }
    }
  });

  it("derives wind speed from the three decoded components", () => {
    const frame = decodeFrameV2(encodeFrameV2({ value: (field) => field === 4 ? 3 : field === 5 ? 4 : field === 6 ? 12 : 1 })) as VisualFrameV2;
    expect(levelSlice(frame, "wind_speed", 1)[4]).toBeCloseTo(Math.hypot(3, 4, 12), 12);
    expect(frameV2Field("wind_speed").kind).toBe("derived");
  });

  it("adapts a level into a renderable dataset without changing the grid identity", () => {
    const frame = decodeFrameV2(encodeFrameV2()) as VisualFrameV2;
    for (const field of frameV2Fields) {
      const value = frameV2LevelDataset(frame, field.id, 1);
      expect(value.schemaVersion).toBe(1);
      expect(value.grid).toEqual({ topology: "cubed_sphere", mapping: "equiangular_gnomonic_v1", cellsPerPanel: 2, flattenOrder: "panel_major_then_j_then_i" });
      expect(value.frame).toEqual({ timeSeconds: 120, step: 12, configFingerprint: "fingerprint-1" });
      expect(value.fields).toHaveLength(1);
      expect(value.fields[0].id).toBe(field.id);
      expect(value.fields[0].unit).toBe(field.unit);
      expect(value.fields[0].provenance).toBe("FrameV2");
      expect(value.fields[0].values).toHaveLength(24);
      expect([...value.fields[0].values]).toEqual(field.volume
        ? [...levelSlice(frame, field.id as Exclude<FrameV2FieldId, "surface_pressure">, 1)]
        : [...frame.surfacePressurePa]);
    }
  });

  it("treats surface pressure as a level-free field", () => {
    const frame = decodeFrameV2(encodeFrameV2()) as VisualFrameV2;
    expect(frameV2Field("surface_pressure").volume).toBe(false);
    // The level argument is ignored rather than validated, because the viewer disables the
    // level selector for a surface field instead of copying it down the column.
    expect([...frameV2LevelDataset(frame, "surface_pressure", 0).fields[0].values])
      .toEqual([...frameV2LevelDataset(frame, "surface_pressure", 2).fields[0].values]);
    expect(frameV2Sample(frame, "surface_pressure", 5, 0)).toEqual({ value: 5000, pressurePa: null });
    expect(() => frameV2Field("depth")).toThrow(ProtocolError);
  });

  it("rejects a frame whose header, shape, or payload does not match the contract", () => {
    expect(() => decodeFrameV2(encodeFrameV2({ schemaVersion: 3 }))).toThrow(/invalid FrameV2 header/);
    expect(() => decodeFrameV2(encodeFrameV2({ cellsPerPanel: 25 }))).toThrow(/invalid FrameV2 shape/);
    expect(() => decodeFrameV2(encodeFrameV2({ levels: 31 }))).toThrow(/invalid FrameV2 shape/);
    expect(() => decodeFrameV2(encodeFrameV2({ declaredCells: 23 }))).toThrow(/invalid FrameV2 shape/);
    const truncated = encodeFrameV2().slice(0, -8);
    expect(() => decodeFrameV2(truncated)).toThrow(/invalid FrameV2 shape/);
    const complete = encodeFrameV2();
    const extended = new Uint8Array(complete.byteLength + 8);
    extended.set(new Uint8Array(complete));
    expect(() => decodeFrameV2(extended.buffer)).toThrow(/invalid FrameV2 shape/);
    expect(() => decodeFrameV2(encodeFrameV2({ value: () => Number.NaN })))
      .toThrow(/non-finite/);
    expect(() => decodeFrameV2(encodeFrameV2(), "other")).toThrow(/fingerprint mismatch/);
    expect(() => decodeFrameV2(new ArrayBuffer(16))).toThrow(/invalid FrameV2 header/);
  });

  it("rejects an out-of-range level or cell instead of reading past the payload", () => {
    const frame = decodeFrameV2(encodeFrameV2()) as VisualFrameV2;
    for (const level of [-1, 3, 1.5, Number.NaN]) {
      expect(() => levelSlice(frame, "temperature", level)).toThrow(/level is out of range/);
    }
    for (const cell of [-1, 24, 0.5]) {
      expect(() => selectedColumn(frame, "temperature", cell)).toThrow(/cell is out of range/);
      expect(() => frameV2Sample(frame, "temperature", cell, 0)).toThrow(/cell is out of range/);
    }
    expect(() => frameV2Sample(frame, "temperature", 0, 3)).toThrow(/level is out of range/);
  });
});

describe("preset descriptors", () => {
  it("keeps an undescribed preset on the shallow-water contract", () => {
    const descriptors = presetDescriptors({ protocolVersion: 1, presets: ["rest"], supportedEdits: ["gaussian_depth"] });
    expect(descriptors).toEqual([shallowWaterPresetDescriptor("rest")]);
    expect(descriptors[0].maximumCellsPerPanel).toBe(96);
    expect(descriptors[0].frameSchemaVersion).toBe(1);
  });

  it("bounds an optional vertical resolution override", () => {
    const dry = { protocolVersion: 1 as const, presetId: "dry", run: { endTimeSeconds: 10, maximumTimeStepSeconds: 1, frameIntervalSteps: 1 }, initialCondition: { edits: [] } };
    // Omitting levels keeps the preset coordinate; the browser bound is only the range the
    // FrameV2 header and the C++ control contract can carry.
    expect(validateRunRequest({ ...dry, grid: { cellsPerPanel: 4 } })).toBeTruthy();
    expect(validateRunRequest({ ...dry, grid: { cellsPerPanel: 4, levels: 1 } })).toBeTruthy();
    expect(validateRunRequest({ ...dry, grid: { cellsPerPanel: 4, levels: 30 } })).toBeTruthy();
    for (const levels of [0, 31, 2.5, -1, "8", null]) {
      expect(() => validateRunRequest({ ...dry, grid: { cellsPerPanel: 4, levels } }))
        .toThrow(/levels is outside the supported range/);
    }
  });

  it("returns described presets in the advertised order and ignores extras", () => {
    const dry = { id: "dry", modelKind: "dry_hydrostatic" as const, frameSchemaVersion: 2 as const, levels: 8, maximumLevels: 30, supportedEdits: [], maximumCellsPerPanel: 24 };
    const descriptors = presetDescriptors({
      protocolVersion: 1, presets: ["rest", "dry"], supportedEdits: ["gaussian_depth"],
      presetDetails: [dry, { ...dry, id: "unlisted" }],
    });
    expect(descriptors.map((preset) => preset.id)).toEqual(["rest", "dry"]);
    expect(descriptors[0].modelKind).toBe("shallow_water");
    expect(descriptors[1]).toEqual(dry);
  });
});
