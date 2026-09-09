#!/usr/bin/env python3
"""Collect compact dry-mixing summaries beside the generated run data."""
import argparse
import csv
import json
import math
from pathlib import Path


def snapshot(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def temperature_rms(first, second):
    weight = 0.0
    total = 0.0
    for a, b in zip(first, second):
        w = float(a['area_m2']) * float(a['air_mass_kg_m2'])
        weight += w
        total += w * (float(a['temperature_k']) - float(b['temperature_k'])) ** 2
    return math.sqrt(total / weight)


def read_log(path):
    return {key: float(value) for key, value in
            (line.split('=', 1) for line in path.read_text().splitlines()
             if '=' in line)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path,
                        default=Path('output/phase12-comparison'))
    parser.add_argument('--thirty-day', type=Path, default=Path('output/phase12-30day'))
    parser.add_argument('--pilot', type=Path,
                        default=Path('output/phase12-pilot/pilot.csv'))
    parser.add_argument('--destination', type=Path,
                        default=Path('output/reports/dry-mixing-results.json'))
    args = parser.parse_args()
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    result = {
        'one_day': {name: entry.get('metrics', entry)
                    for name, entry in
                    json.loads((args.comparison / 'results.json').read_text()).items()},
        'one_day_comparisons':
            json.loads((args.comparison / 'comparisons.json').read_text()),
    }
    coarse = args.thirty_day / 'dt1800.csv'
    fine = args.thirty_day / 'dt450.csv'
    if coarse.exists() and fine.exists():
        result['thirty_day'] = {
            'dt1800': read_log(coarse.with_suffix('.log')),
            'dt450': read_log(fine.with_suffix('.log')),
            'temperature_rms_k': temperature_rms(snapshot(coarse), snapshot(fine)),
        }
    summary = Path(str(args.pilot) + '.summary.json')
    if summary.exists():
        result['pilot'] = json.loads(summary.read_text())
        result['pilot']['run'] = read_log(args.pilot.with_suffix('.log'))
    args.destination.write_text(json.dumps(result, indent=2) + '\n')
    print(args.destination)


if __name__ == '__main__':
    raise SystemExit(main())
