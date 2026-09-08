#!/usr/bin/env python3
"""Summarize segmented Phase 13 pilot diagnostics without dependencies."""

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def config_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def diagnostic_paths(output, name):
    paths = sorted((output / "segments").glob(f"segment-*/{name}"))
    current = output / name
    if current.exists():
        paths.append(current)
    return paths


def rows(output, name):
    answer = []
    for path in diagnostic_paths(output, name):
        with path.open() as stream:
            answer.extend(dict(row, _source=str(path)) for row in csv.DictReader(stream))
    return answer


def number(row, key):
    return float(row[key])


def slope_per_day(samples, key):
    if len(samples) < 2:
        return None
    x = [number(row, "time_s") / 86400.0 for row in samples]
    y = [number(row, key) for row in samples]
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0:
        return None
    return sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y)) / denominator


def aggregate_moisture(samples, area):
    seconds = sum(number(row, "interval_seconds") for row in samples)
    evaporation = sum(number(row, "net_evaporation_kg") for row in samples)
    convective = sum(number(row, "convective_precipitation_kg") for row in samples)
    grid_scale = sum(number(row, "grid_scale_precipitation_kg") for row in samples)
    precipitation = convective + grid_scale
    weighted_pw = sum(number(row, "precipitable_water_kg_m2") *
                      number(row, "interval_seconds") for row in samples)
    return {
        "seconds": seconds,
        "mean_precipitable_water_kg_m2": weighted_pw / seconds if seconds else None,
        "evaporation_kg_m2_s": evaporation / (area * seconds) if seconds else None,
        "precipitation_kg_m2_s": precipitation / (area * seconds) if seconds else None,
        "convective_precipitation_fraction":
            convective / precipitation if precipitation else 0.0,
    }


def block_statistics(samples, area, start_time, block_days):
    blocks = {}
    block_seconds = block_days * 86400.0
    for row in samples:
        index = int((number(row, "time_s") - start_time) // block_seconds)
        blocks.setdefault(index, []).append(row)
    result = []
    for index, selected in sorted(blocks.items()):
        entry = aggregate_moisture(selected, area)
        entry.update({"start_day": index * block_days,
                      "end_day": (index + 1) * block_days})
        result.append(entry)
    return result


def radiation_mean(output, area, minimum_time):
    energy = 0.0
    seconds = 0.0
    for path in diagnostic_paths(output, "radiation_diagnostics.csv"):
        previous = None
        with path.open() as stream:
            for row in csv.DictReader(stream):
                time = number(row, "time_s")
                if previous is None:
                    previous = number(row, "segment_start_time_s")
                interval = time - previous
                previous = time
                if time >= minimum_time and interval > 0:
                    energy += number(row, "interval_toa_net_upward_j")
                    seconds += interval
    return energy / (area * seconds) if seconds else None


def config_fingerprint(output):
    checkpoint = output / "pilot.chk"
    if checkpoint.exists():
        values = config_values(checkpoint)
        if "config_fingerprint" in values:
            return values["config_fingerprint"]
    for path in diagnostic_paths(output, "run-metadata.txt"):
        values = config_values(path)
        if "metadata.config_fingerprint" in values:
            return values["metadata.config_fingerprint"]
    return None


def post_spinup_climate(output, minimum_time):
    selected = []
    for path in diagnostic_paths(output, "climate_statistics.csv"):
        moisture_path = path.parent / "moisture_diagnostics.csv"
        if not moisture_path.exists():
            continue
        with moisture_path.open() as stream:
            moisture = list(csv.DictReader(stream))
        if not moisture or number(moisture[-1], "time_s") < minimum_time:
            continue
        with path.open() as stream:
            selected.extend(row for row in csv.DictReader(stream)
                            if number(row, "accumulated_time_s") > 0.0)
    return {
        "minimum_zonal_mean_temperature_k": min(
            (number(row, "mean_temperature_k") for row in selected), default=None),
        "maximum_zonal_mean_temperature_k": max(
            (number(row, "mean_temperature_k") for row in selected), default=None),
        "minimum_zonal_mean_wind_m_s": min(
            (number(row, "mean_zonal_wind_m_s") for row in selected), default=None),
        "maximum_zonal_mean_wind_m_s": max(
            (number(row, "mean_zonal_wind_m_s") for row in selected), default=None),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "output/phase13_moist_pilot")
    parser.add_argument("--config", type=Path,
                        default=ROOT / "configs/phase13_moist_pilot.cfg")
    parser.add_argument("--spinup-days", type=float, default=200.0)
    parser.add_argument("--block-days", type=float, default=100.0)
    args = parser.parse_args()
    values = config_values(args.config)
    start_time = float(values["run.start_time_s"])
    end_time = float(values["run.end_time_s"])
    radius = float(values["planet.radius_m"])
    area = 4.0 * math.pi * radius * radius
    moisture = rows(args.output, "moisture_diagnostics.csv")
    if not moisture:
        parser.error(f"no moisture diagnostics below {args.output}")
    moisture.sort(key=lambda row: (number(row, "time_s"), number(row, "step")))
    minimum_time = start_time + args.spinup_days * 86400.0
    selected = [row for row in moisture if number(row, "time_s") >= minimum_time]
    if not selected:
        parser.error("no moisture samples remain after the spin-up cutoff")
    aggregate = aggregate_moisture(selected, area)
    first = moisture[0]
    last = moisture[-1]
    system_water = lambda row: (
        number(row, "atmospheric_water_kg") + number(row, "land_water_kg") +
        number(row, "cumulative_ocean_water_change_kg") +
        number(row, "cumulative_external_outflow_kg"))
    exchange_scale = max(
        1.0, system_water(first),
        sum(abs(number(row, "net_evaporation_kg")) +
            number(row, "convective_precipitation_kg") +
            number(row, "grid_scale_precipitation_kg") +
            number(row, "runoff_kg") for row in moisture))
    summed_water_residual = sum(
        number(row, "water_budget_residual_kg") for row in moisture)
    tracer = rows(args.output, "tracer_diagnostics.csv")
    physics = rows(args.output, "physics_diagnostics.csv")
    convection = rows(args.output, "moist_convection_diagnostics.csv")
    radiation = [row for row in rows(args.output, "radiation_diagnostics.csv")
                 if number(row, "time_s") >= minimum_time]
    result = {
        "config_fingerprint": config_fingerprint(args.output),
        "complete": number(last, "time_s") >= end_time,
        "final_day": (number(last, "time_s") - start_time) / 86400.0,
        "final_step": int(number(last, "step")),
        "spinup_days": args.spinup_days,
        "post_spinup": aggregate,
        "toa_net_upward_w_m2": radiation_mean(args.output, area, minimum_time),
        "water": {
            "system_drift_kg": system_water(last) - system_water(first),
            "system_drift_ratio": abs(system_water(last) - system_water(first)) /
                                  exchange_scale,
            "sum_budget_residual_kg": summed_water_residual,
            "sum_budget_residual_ratio": abs(summed_water_residual) /
                                         exchange_scale,
            "maximum_interval_residual_kg": max(
                abs(number(row, "water_budget_residual_kg")) for row in moisture),
            "atmospheric_water_trend_kg_per_day":
                slope_per_day(selected, "atmospheric_water_kg"),
            "land_water_trend_kg_per_day": slope_per_day(selected, "land_water_kg"),
            "precipitable_water_trend_kg_m2_per_day":
                slope_per_day(selected, "precipitable_water_kg_m2"),
        },
        "health": {
            "minimum_tracer_mixing_ratio": min(
                (number(row, "minimum_mixing_ratio") for row in tracer), default=None),
            "maximum_vapor_mixing_ratio": max(
                number(row, "maximum_vapor_mixing_ratio") for row in moisture),
            "maximum_relative_humidity": max(
                number(row, "maximum_relative_humidity") for row in moisture),
            "maximum_dilute_exceedance_area_fraction": max(
                number(row, "dilute_limit_exceedance_area_fraction")
                for row in moisture),
            "minimum_temperature_k": min(
                (number(row, "minimum_temperature_k") for row in physics),
                default=None),
            "maximum_wind_m_s": max(
                (number(row, "maximum_wind_m_s") for row in physics), default=None),
            "maximum_non_finite_count": max(
                (number(row, "non_finite_count") for row in physics), default=None),
            "maximum_temperature_increment_k": max(
                (number(row, "maximum_temperature_increment_k")
                 for row in convection), default=None),
            "maximum_moist_enthalpy_residual_j": max(
                (abs(number(row, "moist_enthalpy_budget_residual_j"))
                 for row in convection), default=None),
        },
        "work": {
            "physics_substeps": sum(int(number(row, "physics_substep_count"))
                                    for row in moisture),
            "physics_retries": sum(int(number(row, "physics_retry_count"))
                                   for row in moisture),
            "cfl_retries_post_spinup": sum(
                int(number(row, "cfl_retry_count")) for row in radiation),
            "invariant_retries_post_spinup": sum(
                int(number(row, "invariant_retry_count")) for row in radiation),
            "solver_retries_post_spinup": sum(
                int(number(row, "solver_retry_count")) for row in radiation),
        },
        "post_spinup_climate_bins": post_spinup_climate(args.output, minimum_time),
        "blocks": block_statistics(selected, area, start_time, args.block_days),
    }
    destination = args.output / "pilot-summary.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(destination)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
