import { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { GaussianDepthEditV1, UnitVector, VisualDatasetV1 } from "@myplanetsim/protocol";
import { createGlobeGeometry, createGridGeometry, globeTrianglesPerCell, pickedCellFromFace, updateGlobeColors } from "./globe";

interface Globe3DProps {
  readonly dataset: VisualDatasetV1;
  readonly fieldId: string;
  readonly gridMode: "off" | "panel_seams" | "all_cells";
  readonly onPick?: (cell: number | null, origin?: UnitVector) => void;
  readonly edits: readonly GaussianDepthEditV1[];
}

interface GlobeRuntime {
  readonly surface: THREE.Mesh;
  readonly grid: THREE.LineSegments;
  readonly editGroup: THREE.Group;
}

export function Globe3D({ dataset, fieldId, gridMode, onPick, edits }: Globe3DProps) {
  const host = useRef<HTMLDivElement>(null);
  const runtime = useRef<GlobeRuntime | null>(null);
  const datasetRef = useRef(dataset); datasetRef.current = dataset;
  const onPickRef = useRef(onPick); onPickRef.current = onPick;
  const [renderError, setRenderError] = useState<string | null>(null);

  useEffect(() => {
    if (!host.current) return undefined;
    setRenderError(null);
    const scene = new THREE.Scene(); scene.background = new THREE.Color("#101820");
    const camera = new THREE.PerspectiveCamera(35, 1, 0.1, 10); camera.position.z = 3;
    let renderer: THREE.WebGLRenderer;
    try { renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false }); }
    catch (error) { setRenderError(error instanceof Error ? error.message : "WebGL initialization failed"); return undefined; }
    renderer.setClearColor("#101820", 1); renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    host.current.replaceChildren(renderer.domElement);

    const surface = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    const grid = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: "#ffffff", transparent: true, opacity: 0.55 })); grid.scale.setScalar(1.002);
    const editGroup = new THREE.Group(); const globe = new THREE.Group(); globe.add(surface, grid, editGroup); scene.add(globe);
    runtime.current = { surface, grid, editGroup };

    const raycaster = new THREE.Raycaster(); const pointer = new THREE.Vector2();
    let pointerActive = false; let dragging = false; let previousX = 0; let previousY = 0;
    const resize = () => {
      const bounds = host.current?.getBoundingClientRect(); if (!bounds) return;
      const width = Math.max(1, bounds.width); const height = Math.max(1, bounds.height);
      camera.aspect = width / height; camera.updateProjectionMatrix(); renderer.setSize(width, height, false); renderer.render(scene, camera);
    };
    const pointerDown = (event: PointerEvent) => { pointerActive = true; dragging = false; previousX = event.clientX; previousY = event.clientY; (event.currentTarget as HTMLElement).setPointerCapture(event.pointerId); };
    const pointerMove = (event: PointerEvent) => {
      if (!pointerActive) return;
      const dx = event.clientX - previousX; const dy = event.clientY - previousY;
      if (Math.abs(dx) + Math.abs(dy) > 1) dragging = true;
      globe.rotation.y += dx * 0.01; globe.rotation.x += dy * 0.01; previousX = event.clientX; previousY = event.clientY;
    };
    const wheel = (event: WheelEvent) => { camera.position.z = THREE.MathUtils.clamp(camera.position.z + event.deltaY * 0.002, 1.5, 6); };
    const pointerUp = (event: PointerEvent) => {
      if (!dragging) {
        const bounds = renderer.domElement.getBoundingClientRect();
        pointer.x = ((event.clientX - bounds.left) / bounds.width) * 2 - 1; pointer.y = -((event.clientY - bounds.top) / bounds.height) * 2 + 1;
        raycaster.setFromCamera(pointer, camera); const intersection = raycaster.intersectObject(surface)[0];
        const localPoint = intersection ? surface.worldToLocal(intersection.point.clone()).normalize() : null;
        const origin: UnitVector | undefined = localPoint ? [localPoint.x, localPoint.y, localPoint.z] : undefined;
        const currentDataset = datasetRef.current;
        onPickRef.current?.(pickedCellFromFace(intersection?.faceIndex ?? null, globeTrianglesPerCell(currentDataset.grid.cellsPerPanel)), origin);
      }
      pointerActive = false;
    };
    const pointerCancel = () => { pointerActive = false; dragging = false; };
    const element = renderer.domElement;
    element.addEventListener("pointerdown", pointerDown); element.addEventListener("pointermove", pointerMove);
    element.addEventListener("pointerup", pointerUp); element.addEventListener("pointercancel", pointerCancel); element.addEventListener("wheel", wheel, { passive: true });
    const observer = new ResizeObserver(resize); observer.observe(host.current); resize();
    let animationFrame = 0; const render = () => { renderer.render(scene, camera); animationFrame = requestAnimationFrame(render); }; render();
    return () => {
      cancelAnimationFrame(animationFrame); observer.disconnect();
      element.removeEventListener("pointerdown", pointerDown); element.removeEventListener("pointermove", pointerMove);
      element.removeEventListener("pointerup", pointerUp); element.removeEventListener("pointercancel", pointerCancel); element.removeEventListener("wheel", wheel);
      disposeObject(editGroup); surface.geometry.dispose(); (surface.material as THREE.Material).dispose(); grid.geometry.dispose(); (grid.material as THREE.Material).dispose();
      runtime.current = null; renderer.dispose(); renderer.domElement.remove();
    };
  }, []);

  useEffect(() => {
    if (!runtime.current) return;
    const surface = runtime.current.surface;
    if (updateGlobeColors(surface.geometry, dataset, fieldId)) return;
    const previous = surface.geometry; surface.geometry = createGlobeGeometry(dataset, fieldId); previous.dispose();
  }, [dataset, fieldId]);

  useEffect(() => {
    if (!runtime.current) return;
    const previous = runtime.current.grid.geometry; runtime.current.grid.geometry = createGridGeometry(dataset.grid.cellsPerPanel, gridMode); previous.dispose();
  }, [dataset.grid.cellsPerPanel, gridMode]);

  useEffect(() => {
    if (!runtime.current) return;
    const group = runtime.current.editGroup; disposeObject(group); group.clear();
    for (const edit of edits) {
      const color = edit.amplitudeMeters >= 0 ? "#ffcf56" : "#ff6fae";
      const marker = new THREE.Mesh(new THREE.SphereGeometry(0.025, 12, 8), new THREE.MeshBasicMaterial({ color })); marker.position.set(...edit.centerUnit).multiplyScalar(1.035);
      const ring = new THREE.LineLoop(editRingGeometry(edit.centerUnit, edit.sigmaRadians), new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.9 }));
      group.add(marker, ring);
    }
  }, [edits]);

  return <div aria-label="3D globe" className="globe-3d"><div ref={host} className="globe-canvas" />{renderError ? <p className="render-error">3D view unavailable: {renderError}</p> : null}</div>;
}

function disposeObject(group: THREE.Group) {
  for (const child of group.children) {
    const renderable = child as THREE.Mesh | THREE.Line;
    renderable.geometry?.dispose();
    const materials = Array.isArray(renderable.material) ? renderable.material : [renderable.material];
    for (const material of materials) material?.dispose();
  }
}

function editRingGeometry(origin: UnitVector, sigma: number): THREE.BufferGeometry {
  const normal = new THREE.Vector3(...origin).normalize();
  const reference = Math.abs(normal.z) < 0.9 ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(0, 1, 0);
  const firstTangent = reference.cross(normal).normalize(); const secondTangent = new THREE.Vector3().crossVectors(normal, firstTangent).normalize();
  const points = Array.from({ length: 96 }, (_, index) => {
    const angle = 2 * Math.PI * index / 96;
    return normal.clone().multiplyScalar(Math.cos(sigma)).addScaledVector(firstTangent, Math.sin(sigma) * Math.cos(angle)).addScaledVector(secondTangent, Math.sin(sigma) * Math.sin(angle)).multiplyScalar(1.035);
  });
  return new THREE.BufferGeometry().setFromPoints(points);
}
