import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  // The unit suite is vitest over src. test/ holds the node --test browser gate, which
  // drives a built page and a live gateway, so vitest must not collect it.
  test: { include: ["src/**/*.test.ts"] },
});
