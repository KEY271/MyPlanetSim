#!/usr/bin/env python3
"""Run and compare the registered Phase 13 moist global experiment matrix.

Every case is a separate process.  Derived configs, logs, daily samples, final
snapshots, metrics, and pairwise comparisons remain in the selected output tree.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def tracer_overrides(count):
    names = ([f"passive_{index}" for index in range(1, count)] + ["vapor"])
    result = {"tracers.names": ",".join(names)}
    for name in names[:-1]:
        result.update({
            f"tracers.{name}.role": "passive",
            f"tracers.{name}.initial_mixing_ratio": 0.001,
            f"tracers.{name}.require_nonnegative": "true",
            f"tracers.{name}.horizontal_diffusion": "true",
        })
    result.update({
        "tracers.vapor.role": "water_vapor",
        "tracers.vapor.initial_relative_humidity": 0.6,
        "tracers.vapor.require_nonnegative": "true",
        "tracers.vapor.horizontal_diffusion": "true",
    })
    return result


def passive_tracer_overrides(count):
    names = [f"passive_{index}" for index in range(1, count + 1)]
    result = {"tracers.names": ",".join(names)}
    for index, name in enumerate(names, start=1):
        result.update({
            f"tracers.{name}.role": "passive",
            f"tracers.{name}.initial_mixing_ratio": 0.001 * index,
            f"tracers.{name}.require_nonnegative": "true",
            f"tracers.{name}.horizontal_diffusion": "true",
        })
    return result


def build_cases(suite):
    # source, dt, iterations, max physics substep, overrides, removed prefixes
    moist = ("phase13_moist_pilot.cfg", 450, 5, 75, {}, ())
    cases = {}
    if suite == "one-day":
        cases.update({
            "dry-zero-water": (
                "phase12_gray_pilot.cfg", 450, 5, 75,
                {"surface.uniform_land_fraction": 0}, ()),
            "moist-lsc": moist[:4] + ({"convection.kind": "none"}, ("convection.",)),
            "moist-dry-adjustment": moist[:4] + (
                {"convection.kind": "dry_adjustment",
                 "convection.stability_tolerance_k": "1e-10"},
                ("convection.",)),
            "moist-sbm": moist,
            "land-wet": moist[:4] + (
                {"surface.uniform_land_fraction": 0.5,
                 "surface.hydrology.initial_fraction": 1}, ()),
            "land-dry": moist[:4] + (
                {"surface.uniform_land_fraction": 0.5,
                 "surface.hydrology.initial_fraction": 0}, ()),
            "all-land-closed": moist[:4] + (
                {"surface.uniform_land_fraction": 1,
                 "surface.hydrology.initial_fraction": 0.5}, ()),
            "earth-terrain": moist[:4] + (
                {"surface.geography": "earth",
                 "surface.input_file": "../../../data/earth/earth_surface_derived.csv",
                 "surface.input_fingerprint_fnv1a64": "f6c9b8886a03d401",
                 "surface.quadrature_order": 2,
                 "surface.smoothing_passes": 0},
                ("surface.geography", "surface.uniform_land_fraction")),
            "synchronous": moist[:4] + (
                {"planet.rotation_rate_rad_s": "1.991021277657232e-07",
                 "surface.uniform_land_fraction": 0.5}, ()),
        })
        for dt, iterations, substep in ((450, 5, 75), (900, 5, 150),
                                         (1800, 3, 300), (1800, 5, 300)):
            cases[f"time-{dt}-i{iterations}"] = moist[:1] + (
                dt, iterations, substep, {}, ())
        cases["explicit"] = moist[:1] + (30, 0, 30, {}, ())
        for rh in (0.7, 0.9):
            cases[f"rh-{rh}"] = moist[:4] + (
                {"convection.reference_relative_humidity": rh}, ())
        for tau in (3600, 14400):
            cases[f"tau-{tau}"] = moist[:4] + (
                {"convection.relaxation_time_s": tau}, ())
        for capacity in (50, 150, 300):
            cases[f"bucket-{capacity}"] = moist[:4] + (
                {"surface.uniform_land_fraction": 0.5,
                 "surface.hydrology.capacity_kg_m2": capacity}, ())
        for substep in (150, 300):
            cases[f"substep-{substep}"] = moist[:3] + (substep, {}, ())
        for count in (1, 2, 4):
            cases[f"components-{count}"] = moist[:4] + (
                tracer_overrides(count), ("tracers.",))
        for count in (2, 4):
            cases[f"dry-components-{count}"] = (
                "phase12_gray_pilot.cfg", 450, 5, 75,
                {"surface.uniform_land_fraction": 0,
                 **passive_tracer_overrides(count)}, ())
    else:
        # Thirty days are used for the registered timestep check and one-at-a-time
        # coefficient sensitivities. Geography remains bounded to the one-day gate.
        cases.update({
            "standard-dt450": moist,
            "standard-dt900": moist[:1] + (900, 5, 150, {}, ()),
            "standard-dt1800-i3": moist[:1] + (1800, 3, 300, {}, ()),
            "standard-dt1800-i5": moist[:1] + (1800, 5, 300, {}, ()),
        })
        for rh in (0.7, 0.9):
            cases[f"rh-{rh}"] = moist[:4] + (
                {"convection.reference_relative_humidity": rh}, ())
        for tau in (3600, 14400):
            cases[f"tau-{tau}"] = moist[:4] + (
                {"convection.relaxation_time_s": tau}, ())
        for capacity in (50, 300):
            cases[f"bucket-{capacity}"] = moist[:4] + (
                {"surface.uniform_land_fraction": 0.5,
                 "surface.hydrology.capacity_kg_m2": capacity}, ())
    return cases


def derived_config(source, overrides, removed_prefixes, destination):
    remaining = {str(key): value for key, value in overrides.items()}
    output = []
    for line in source.read_text().splitlines():
        key = line.split("=", 1)[0].strip() if "=" in line else None
        if key and any(key.startswith(prefix) for prefix in removed_prefixes):
            continue
        if key in remaining:
            output.append(f"{key} = {remaining.pop(key)}")
        else:
            output.append(line)
    output.extend(f"{key} = {value}" for key, value in remaining.items())
    destination.write_text("\n".join(output) + "\n")


def parse_metrics(text):
    result = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            result[key] = float(value)
        except ValueError:
            pass
    return result


def read_snapshot(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def temperature_rms(first, second):
    if len(first) != len(second):
        return None
    weight = 0.0
    total = 0.0
    for a, b in zip(first, second):
        if (a["cell"], a["level"]) != (b["cell"], b["level"]):
            return None
        mass = float(a["area_m2"]) * float(a["air_mass_kg_m2"])
        weight += mass
        total += mass * (float(a["temperature_k"]) -
                         float(b["temperature_k"])) ** 2
    return math.sqrt(total / weight)


def relative_difference(first, second):
    return None if second == 0 else (first - second) / abs(second)


def comparisons(results, output, suite):
    reference = "time-450-i5" if suite == "one-day" else "standard-dt450"
    names = [name for name, entry in results.items()
             if name != reference and entry["returncode"] == 0]
    answer = {}
    for name in names:
        if results.get(reference, {}).get("returncode") != 0:
            break
        first = results[name]["metrics"]
        second = results[reference]["metrics"]
        rms = temperature_rms(read_snapshot(output / f"{name}.csv"),
                              read_snapshot(output / f"{reference}.csv"))
        entry = {"temperature_rms_k": rms}
        for key in ("precipitable_water_kg_m2", "evaporation_kg_m2_s",
                    "precipitation_kg_m2_s", "toa_net_upward_w_m2"):
            if key in first and key in second:
                entry[f"{key}_difference"] = first[key] - second[key]
                entry[f"{key}_relative"] = relative_difference(first[key], second[key])
        answer[f"{name}_vs_{reference}"] = entry
    return answer


def run_suite(args, suite, days):
    output = args.output / suite
    output.mkdir(parents=True, exist_ok=True)
    results_path = output / "results.json"
    results = (json.loads(results_path.read_text())
               if args.case and results_path.exists() else {})
    failures = 0
    for name, (source_name, dt, iterations, substep, overrides,
               removed_prefixes) in build_cases(suite).items():
        if args.case and name != args.case:
            continue
        source = args.configs / source_name
        config = output / f"{name}.cfg"
        combined = dict(overrides)
        if source_name == "phase13_moist_pilot.cfg":
            combined["moisture.maximum_physics_substep_s"] = substep
        derived_config(source, combined, removed_prefixes, config)
        snapshot = output / f"{name}.csv"
        command = [str(args.binary), str(config), str(days), "6", "20", str(dt),
                   str(iterations), "20000", "5", str(snapshot)]
        print(f"[{suite}] {name}", flush=True)
        completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                                   check=False)
        (output / f"{name}.log").write_text(completed.stdout + completed.stderr)
        entry = {"command": command, "returncode": completed.returncode,
                 "days": days, "overrides": combined}
        if completed.returncode == 0:
            entry["metrics"] = parse_metrics(completed.stdout)
        else:
            entry["error"] = completed.stderr.strip()
            failures += 1
        results[name] = entry
    results_path.write_text(json.dumps(results, indent=2) + "\n")
    compared = comparisons(results, output, suite)
    (output / "comparisons.json").write_text(json.dumps(compared, indent=2) + "\n")
    print(f"[{suite}] {len(results) - failures}/{len(results)} cases completed")
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/release/tools/benchmark_moist")
    parser.add_argument("--configs", type=Path, default=ROOT / "configs")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "output/phase13-comparison")
    parser.add_argument("--suite", choices=("one-day", "thirty-day", "all"),
                        default="all")
    parser.add_argument("--case", help="run one named case and merge its result")
    args = parser.parse_args()
    failures = 0
    if args.suite in ("one-day", "all"):
        failures += run_suite(args, "one-day", 1.0)
    if args.suite in ("thirty-day", "all"):
        failures += run_suite(args, "thirty-day", 30.0)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
