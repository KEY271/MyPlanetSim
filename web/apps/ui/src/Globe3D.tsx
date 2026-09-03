import { useEffect, useRef } from "react";
import * as THREE from "three";
import { VisualDatasetV1 } from "@myplanetsim/protocol";
import { createGlobeGeometry, createGridGeometry, pickedCellFromFace } from "./globe";

interface Globe3DProps {
  readonly dataset: VisualDatasetV1;
  readonly fieldId: string;
  readonly gridMode: "off" | "panel_seams" | "all_cells";
  readonly onPick?: (cell: number | null) => void;
}

export function Globe3D({ dataset, fieldId, gridMode, onPick }: Globe3DProps) {
  const host = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!host.current) return undefined;
    const scene = new THREE.Scene();
    scene.background = new THREE.Color("#101820");
    const camera = new THREE.PerspectiveCamera(35, 1, 0.1, 10);
    camera.position.z = 3;
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    host.current.replaceChildren(renderer.domElement);
    const surface = new THREE.Mesh(createGlobeGeometry(dataset, fieldId),
      new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    const grid = new THREE.LineSegments(createGridGeometry(dataset.grid.cellsPerPanel, gridMode),
      new THREE.LineBasicMaterial({ color: "#ffffff", transparent: true, opacity: 0.55 }));
    grid.scale.setScalar(1.002);
    scene.add(surface, grid);
    const raycaster = new THREE.Raycaster();
    const pointer = new THREE.Vector2();
    let dragging = false;
    let previousX = 0;
    let previousY = 0;
    const resize = () => {
      const bounds = host.current?.getBoundingClientRect();
      if (!bounds) return;
      camera.aspect = bounds.width / Math.max(bounds.height, 1);
      camera.updateProjectionMatrix();
      renderer.setSize(bounds.width, bounds.height, false);
    };
    const pointerDown = (event: PointerEvent) => {
      dragging = false;
      previousX = event.clientX;
      previousY = event.clientY;
      (event.currentTarget as HTMLElement).setPointerCapture(event.pointerId);
    };
    const pointerMove = (event: PointerEvent) => {
      const dx = event.clientX - previousX;
      const dy = event.clientY - previousY;
      if (Math.abs(dx) + Math.abs(dy) > 1) dragging = true;
      surface.rotation.y += dx * 0.01;
      surface.rotation.x += dy * 0.01;
      grid.rotation.copy(surface.rotation);
      previousX = event.clientX;
      previousY = event.clientY;
    };
    const wheel = (event: WheelEvent) => {
      camera.position.z = THREE.MathUtils.clamp(camera.position.z + event.deltaY * 0.002, 1.5, 6);
    };
    const pointerUp = (event: PointerEvent) => {
      if (!dragging) {
        const bounds = renderer.domElement.getBoundingClientRect();
        pointer.x = ((event.clientX - bounds.left) / bounds.width) * 2 - 1;
        pointer.y = -((event.clientY - bounds.top) / bounds.height) * 2 + 1;
        raycaster.setFromCamera(pointer, camera);
        const intersection = raycaster.intersectObject(surface)[0];
        onPick?.(pickedCellFromFace(intersection?.faceIndex ?? null));
      }
    };
    const element = renderer.domElement;
    element.addEventListener("pointerdown", pointerDown);
    element.addEventListener("pointermove", pointerMove);
    element.addEventListener("pointerup", pointerUp);
    element.addEventListener("wheel", wheel, { passive: true });
    const observer = new ResizeObserver(resize);
    observer.observe(host.current);
    resize();
    let animationFrame = 0;
    const render = () => { renderer.render(scene, camera); animationFrame = requestAnimationFrame(render); };
    render();
    return () => {
      cancelAnimationFrame(animationFrame);
      observer.disconnect();
      element.removeEventListener("pointerdown", pointerDown);
      element.removeEventListener("pointermove", pointerMove);
      element.removeEventListener("pointerup", pointerUp);
      element.removeEventListener("wheel", wheel);
      surface.geometry.dispose();
      (surface.material as THREE.Material).dispose();
      grid.geometry.dispose();
      (grid.material as THREE.Material).dispose();
      renderer.dispose();
      renderer.domElement.remove();
    };
  }, [dataset, fieldId, gridMode, onPick]);
  return <div ref={host} aria-label="3D globe" className="globe-3d" />;
}
