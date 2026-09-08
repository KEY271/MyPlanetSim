#!/usr/bin/env python3
"""Render Phase 13 column, water-cycle, and bucket diagnostics as an SVG."""

import argparse
import csv
from pathlib import Path


def scale(value, lower, upper, first, last):
    if upper == lower:
        return 0.5 * (first + last)
    return first + (value - lower) * (last - first) / (upper - lower)


def polyline(points, colour, width=2):
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return (f'<polyline points="{coordinates}" fill="none" stroke="{colour}" '
            f'stroke-width="{width}"/>')


def read_rows(path):
    if path is None or not path.exists():
        return []
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def panel_axes(svg, left, right, top, bottom, title, x_label, y_label):
    svg.extend([
        f'<text x="{left}" y="{top - 14}" font-family="sans-serif" '
        f'font-size="16">{title}</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" '
        'stroke="#333"/>',
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" '
        'stroke="#333"/>',
        f'<text x="{0.5 * (left + right) - 70}" y="{bottom + 32}" '
        f'font-family="sans-serif" font-size="12">{x_label}</text>',
        f'<text x="{left - 50}" y="{0.5 * (top + bottom)}" '
        f'transform="rotate(-90 {left - 50} {0.5 * (top + bottom)})" '
        f'font-family="sans-serif" font-size="12">{y_label}</text>',
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--moisture", type=Path,
                        help="moisture_diagnostics.csv (defaults beside input)")
    parser.add_argument("--surface", type=Path,
                        help="moist_surface_state.csv (defaults beside input)")
    args = parser.parse_args()
    destination = args.output or args.input.with_suffix(".svg")
    rows = read_rows(args.input)
    if not rows:
        raise SystemExit("moist column CSV contains no levels")
    moisture_path = args.moisture or args.input.with_name("moisture_diagnostics.csv")
    surface_path = args.surface or args.input.with_name("moist_surface_state.csv")
    moisture_rows = read_rows(moisture_path)
    surface_rows = read_rows(surface_path)
    pressure = [float(row["pressure_pa"]) / 100.0 for row in rows]
    temperature = [float(row["temperature_k"]) for row in rows]
    humidity = [100.0 * float(row["relative_humidity"]) for row in rows]
    left, right, top, bottom = 75.0, 845.0, 55.0, 350.0
    p_min, p_max = min(pressure), max(pressure)
    t_min, t_max = min(temperature), max(temperature)
    if t_min == t_max:
        t_min, t_max = t_min - 1.0, t_max + 1.0
    temperature_points = [
        (scale(t, t_min, t_max, left, right),
         scale(p, p_min, p_max, top, bottom))
        for t, p in zip(temperature, pressure)
    ]
    humidity_points = [
        (scale(rh, 0.0, 100.0, left, right),
         scale(p, p_min, p_max, top, bottom))
        for rh, p in zip(humidity, pressure)
    ]
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="920" height="920" '
        'viewBox="0 0 920 920">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    panel_axes(svg, left, right, top, bottom,
               "Column: temperature (red), RH (blue)",
               "independent temperature and RH scales", "pressure [hPa]")
    svg.extend([polyline(temperature_points, "#d1495b"),
                polyline(humidity_points, "#00798c")])

    water_top, water_bottom = 445.0, 645.0
    panel_axes(svg, left, right, water_top, water_bottom,
               "Water cycle: precipitable water (green), precipitation (purple)",
               "time [day]", "independent scales")
    if moisture_rows:
        time = [float(row["time_s"]) / 86400.0 for row in moisture_rows]
        precipitable = [float(row["precipitable_water_kg_m2"])
                        for row in moisture_rows]
        rain = [float(row["precipitation_mm_day_display"])
                for row in moisture_rows]
        time_min, time_max = min(time), max(time)
        pw_min, pw_max = min(precipitable), max(precipitable)
        rain_min, rain_max = min(rain), max(rain)
        pw_points = [(scale(t, time_min, time_max, left, right),
                      scale(value, pw_min, pw_max, water_bottom, water_top))
                     for t, value in zip(time, precipitable)]
        rain_points = [(scale(t, time_min, time_max, left, right),
                        scale(value, rain_min, rain_max, water_bottom, water_top))
                       for t, value in zip(time, rain)]
        svg.extend([polyline(pw_points, "#2a9d55"),
                    polyline(rain_points, "#7b2cbf")])
    else:
        svg.append('<text x="90" y="545" font-family="sans-serif" '
                   'font-size="13">moisture_diagnostics.csv not found</text>')

    bucket_top, bucket_bottom = 740.0, 860.0
    panel_axes(svg, left, right, bucket_top, bucket_bottom,
               "Final land bucket distribution", "cell id",
               "water [kg m-2 land]")
    if surface_rows:
        cell = [float(row["cell_id"]) for row in surface_rows]
        bucket = [float(row["land_water_kg_m2"]) for row in surface_rows]
        cell_min, cell_max = min(cell), max(cell)
        bucket_min, bucket_max = min(bucket), max(bucket)
        bucket_points = [(scale(index, cell_min, cell_max, left, right),
                          scale(value, bucket_min, bucket_max,
                                bucket_bottom, bucket_top))
                         for index, value in zip(cell, bucket)]
        svg.append(polyline(bucket_points, "#d97706", 1.5))
        svg.extend(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2" fill="#d97706"/>'
                   for x, y in bucket_points)
    else:
        svg.append('<text x="90" y="805" font-family="sans-serif" '
                   'font-size="13">moist_surface_state.csv not found</text>')
    svg.append('</svg>')
    destination.write_text("\n".join(svg) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
