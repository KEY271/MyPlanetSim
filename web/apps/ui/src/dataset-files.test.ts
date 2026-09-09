import { describe, expect, it } from "vitest";
import { browserDatasetFiles } from "./dataset-files";

function file(path: string, contents: string): File {
  return { name: path.split("/").at(-1)!, webkitRelativePath: path,
    text: async () => contents,
    arrayBuffer: async () => new TextEncoder().encode(contents).buffer } as File;
}

describe("browser dataset files", () => {
  it("indexes a selected folder by safe paths", async () => {
    const files = browserDatasetFiles([
      file("dataset/manifest.json", "manifest"),
      file("dataset/means/period_000000.bin", "mean"),
    ]);
    expect(await files.readText("manifest.json")).toBe("manifest");
    expect(new TextDecoder().decode(await files.readBytes("means/period_000000.bin")))
      .toBe("mean");
    await expect(files.readText("missing.bin")).rejects.toThrow(/missing file/);
  });
});
