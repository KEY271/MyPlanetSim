#!/usr/bin/env python3
"""Build the bounded 1-degree surface fixture used by Phase 8 CI.

This repository fixture exercises the complete fractional remap contract without
shipping the multi-gigabyte NOAA inputs.  It is deliberately not a production
ETOPO/GSHHG derivative.  A production product must replace it only after its raw
hashes and generator are recorded in data/earth/manifest.txt.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def is_land(longitude_deg: float, latitude_deg: float) -> bool:
    shapes = (
        (-105.0, 48.0, 45.0, 30.0),  # North America
        (-60.0, -15.0, 22.0, 38.0),  # South America
        (20.0, 7.0, 28.0, 38.0),  # Africa
        (75.0, 48.0, 78.0, 32.0),  # Eurasia
        (134.0, -25.0, 20.0, 14.0),  # Australia
        (-42.0, 73.0, 12.0, 9.0),  # Greenland
    )
    if latitude_deg < -68.0:
        return True
    for centre_lon, centre_lat, width, height in shapes:
        delta_lon = math.remainder(longitude_deg - centre_lon, 360.0)
        if (delta_lon / width) ** 2 + ((latitude_deg - centre_lat) / height) ** 2 < 1:
            return True
    return False


def elevation(longitude_deg: float, latitude_deg: float) -> float:
    ridge = 2600.0 * math.exp(
        -((math.remainder(longitude_deg + 70.0, 360.0) / 9.0) ** 2)
        - (((latitude_deg + 15.0) / 35.0) ** 2)
    )
    plateau = 1800.0 * math.exp(
        -((math.remainder(longitude_deg - 85.0, 360.0) / 18.0) ** 2)
        - (((latitude_deg - 32.0) / 10.0) ** 2)
    )
    return 150.0 + ridge + plateau


def build(output_path: Path) -> None:
    samples = 8
    with output_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(
            "longitude_deg,latitude_deg,mean_surface_height_m,land_fraction\n"
        )
        for longitude_index in range(360):
            longitude = -179.5 + longitude_index
            for latitude_index in range(180):
                latitude = -89.5 + latitude_index
                land_count = 0
                height_sum = 0.0
                for sample_i in range(samples):
                    for sample_j in range(samples):
                        sample_lon = longitude + (sample_i + 0.5) / samples - 0.5
                        sample_lat = latitude + (sample_j + 0.5) / samples - 0.5
                        if is_land(sample_lon, sample_lat):
                            land_count += 1
                            height_sum += elevation(sample_lon, sample_lat)
                count = samples * samples
                fraction = land_count / count
                mean_height = height_sum / count
                output.write(
                    f"{longitude:.1f},{latitude:.1f},{mean_height:.8f},{fraction:.8f}\n"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ci-fixture-output", required=True, type=Path)
    arguments = parser.parse_args()
    arguments.ci_fixture_output.parent.mkdir(parents=True, exist_ok=True)
    build(arguments.ci_fixture_output)


if __name__ == "__main__":
    main()
