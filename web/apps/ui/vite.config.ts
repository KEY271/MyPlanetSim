import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// `just dev` starts the gateway and this dev server together. Requiring the reader to copy a
// tokenised URL made the plain URL Vite prints silently fall back to the offline demo, and
// --open could not fix it because reusing an already-open tab is a same-document navigation.
// When the launcher exports a session, hand it to the client so the printed URL just works.
// This is dev-only and loopback-only; a production build always gets null.
function devSession(command: string) {
  const token = process.env.MPS_SESSION_TOKEN;
  const port = process.env.MPS_GATEWAY_PORT;
  if (command !== "serve" || !token || !port) return null;
  return { token, gateway: `http://127.0.0.1:${port}` };
}

export default defineConfig(({ command }) => ({
  plugins: [react()],
  define: { __MPS_DEV_SESSION__: JSON.stringify(devSession(command)) },
  // The unit suite is vitest over src. test/ holds the node --test browser gate, which
  // drives a built page and a live gateway, so vitest must not collect it.
  test: { include: ["src/**/*.test.ts"] },
}));
