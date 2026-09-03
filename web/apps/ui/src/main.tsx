import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "./style.css";

function App() {
  return (
    <section className="app-shell">
      <p className="eyebrow">MyPlanetSim · Phase 3</p>
      <h1>Interactive shallow-water visualizer</h1>
      <p>Offline visualizer workspace initialized. Simulation controls arrive in the protocol phase.</p>
    </section>
  );
}

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
