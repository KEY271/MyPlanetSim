import { DatasetFiles, ProtocolError, validateDatasetPath } from "@myplanetsim/protocol";

export function browserDatasetFiles(files: FileList | readonly File[]): DatasetFiles {
  const entries = Array.from(files);
  if (entries.length === 0) throw new ProtocolError("invalid_visual_dataset", "no files selected");
  const raw = entries.map((file) => (file.webkitRelativePath || file.name).replace(/^\/+/, ""));
  const firstParts = raw[0].split("/");
  const stripRoot = firstParts.length > 1 && raw.every((path) => path.split("/")[0] === firstParts[0]);
  const indexed = new Map<string, File>();
  entries.forEach((file, index) => {
    const parts = raw[index].split("/");
    const relativePath = validateDatasetPath(stripRoot ? parts.slice(1).join("/") : raw[index]);
    if (indexed.has(relativePath)) throw new ProtocolError("invalid_visual_dataset", `duplicate file ${relativePath}`);
    indexed.set(relativePath, file);
  });
  const file = (path: string) => {
    const result = indexed.get(validateDatasetPath(path));
    if (!result) throw new ProtocolError("invalid_visual_dataset", `missing file ${path}`);
    return result;
  };
  return Object.freeze({
    readText: async (path: string) => file(path).text(),
    readBytes: async (path: string) => file(path).arrayBuffer(),
  });
}
