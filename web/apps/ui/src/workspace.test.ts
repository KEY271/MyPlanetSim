import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";

describe("web workspace", () => {
  it("reserves the full viewport for top controls and side-by-side views", () => {
    const style = readFileSync(new URL("./style.css", import.meta.url), "utf8");
    expect(style).toContain("height: 100dvh");
    expect(style).toContain("grid-template-rows: auto minmax(0, 1fr)");
    expect(style).toContain("grid-template-columns: minmax(0, 1fr) minmax(0, 1fr)");
  });
});
