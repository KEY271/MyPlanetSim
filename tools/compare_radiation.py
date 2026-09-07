#!/usr/bin/env python3
"""Run the registered one-day matrix without third-party dependencies.

All commands, metrics, errors, and snapshots remain in the output directory.
A failed experiment is recorded and makes the runner exit nonzero.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=Path('build/release/tools/benchmark_radiation'))
    parser.add_argument('--output', type=Path, default=Path('output/phase11-comparison'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = {}
    for dt in (1800, 900, 450):
        for iterations in (3, 4, 5):
            cases[f'time-{dt}-i{iterations}'] = ('gray_uniform', 6, 10, dt, iterations, 1000)
    cases['explicit'] = ('gray_uniform', 6, 10, 30, 0, 1000)
    for n, k in ((6, 20), (6, 40), (12, 10), (12, 20)):
        cases[f'space-n{n}-k{k}'] = ('gray_uniform', n, k, 450, 5, 1000)
    for name in ('transparent_surface', 'gray_shortwave', 'gray_land_ocean', 'gray_tidally_locked'):
        cases[name] = (name, 6, 10, 450, 5, 1000)
    for top in (500, 250):
        cases[f'top-{top}'] = ('gray_uniform', 6, 10, 450, 5, top)
    results = {}
    for name, (preset, n, k, dt, iterations, top) in cases.items():
        snapshot = args.output / f'{name}.csv'
        command = [str(args.binary), f'configs/phase11_{preset}.cfg', '1',
                   str(n), str(k), str(dt), str(iterations), str(top), str(snapshot)]
        print(name, flush=True)
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        (args.output / f'{name}.log').write_text(completed.stdout + completed.stderr)
        result = {'command': command, 'returncode': completed.returncode}
        if completed.returncode == 0:
            result['metrics'] = {key: float(value) for key, value in
                                 (line.split('=', 1) for line in completed.stdout.splitlines())}
            with snapshot.open() as stream:
                rows = list(csv.DictReader(stream))
            total_area = sum(float(row['area_m2']) for row in rows)
            result['metrics']['mean_temperature_k'] = sum(float(row['area_m2']) * float(row['temperature_k']) for row in rows) / total_area
        else:
            result['error'] = completed.stderr
        results[name] = result
        (args.output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    comparisons = {}
    pairs = [(key, 'time-450-i5') for key in cases if key.startswith('time-')]
    pairs += [('time-1800-i5', 'explicit'), ('time-1800-i3', 'time-1800-i5'),
              ('time-1800-i4', 'time-1800-i5')]
    for name, reference in pairs:
        if results[name]['returncode'] or results[reference]['returncode']:
            continue
        with (args.output / f'{name}.csv').open() as first_file, (args.output / f'{reference}.csv').open() as second_file:
            first = list(csv.DictReader(first_file))
            second = list(csv.DictReader(second_file))
        assert len(first) == len(second)
        area = sum(float(row['area_m2']) for row in first)
        rms = math.sqrt(sum(float(a['area_m2']) * (float(a['temperature_k']) - float(b['temperature_k']))**2
                            for a, b in zip(first, second)) / area)
        toa = abs(results[name]['metrics']['toa_net_upward_w_m2'] - results[reference]['metrics']['toa_net_upward_w_m2'])
        comparisons[f'{name}_vs_{reference}'] = {'reference': reference, 'temperature_rms_k': rms, 'toa_difference_w_m2': toa,
                             'passes_one_day_accuracy': rms <= 0.5 and toa <= 1.0}
    (args.output / 'comparisons.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    return int(any(item['returncode'] for item in results.values()))


if __name__ == '__main__':
    raise SystemExit(main())
