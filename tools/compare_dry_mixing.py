#!/usr/bin/env python3
"""Run the registered Phase 12 one-day matrix without third-party dependencies.

Every command, metric, error, and snapshot stays in the output directory. Coefficient
variants are written out as derived configuration files so the exact input of each case
is recoverable. A failed experiment is recorded and makes the runner exit nonzero.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess

# name -> (preset, N, K, dt, iterations, top_pa, refinement, leg, overrides)
def build_cases():
    cases = {}
    base = ('gray_pilot', 6, 20, 450, 5, 1000, 5.0, 'all', {})
    # The four-way physics matrix of the plan on one grid and one initial state.
    for leg in ('radiation', 'convection', 'boundary_layer', 'all'):
        cases[f'matrix-{leg}'] = base[:7] + (leg, {})
    # Time-step sensitivity, plus an explicit small-step reference.
    for dt, iterations in ((1800, 3), (1800, 5), (900, 5), (450, 5)):
        cases[f'time-{dt}-i{iterations}'] = ('gray_pilot', 6, 20, dt, iterations, 1000,
                                             5.0, 'all', {})
    cases['explicit'] = ('gray_pilot', 6, 20, 30, 0, 1000, 5.0, 'all', {})
    # Grid sensitivity: horizontal, vertical, near-surface refinement, and model top.
    for n, k in ((6, 20), (6, 40), (12, 20), (12, 40)):
        cases[f'space-n{n}-k{k}'] = ('gray_pilot', n, k, 450, 5, 1000, 5.0, 'all', {})
    cases['uniform-k20'] = ('gray_pilot', 6, 20, 450, 5, 1000, 1.0, 'all', {})
    for top in (500, 250):
        cases[f'top-{top}'] = ('gray_pilot', 6, 20, 450, 5, top, 5.0, 'all', {})
    # Coefficient sensitivity, changed one family at a time.
    for critical in (0.25, 0.5):
        cases[f'richardson-{critical}'] = ('gray_pilot', 6, 20, 450, 5, 1000, 5.0, 'all',
                                           {'boundary_layer.critical_richardson': critical})
    for gustiness in (0, 0.5):
        cases[f'gustiness-{gustiness}'] = ('gray_pilot', 6, 20, 450, 5, 1000, 5.0, 'all',
                                           {'boundary_layer.gustiness_m_s': gustiness})
    cases['prandtl-1.5'] = ('gray_pilot', 6, 20, 450, 5, 1000, 5.0, 'all',
                            {'boundary_layer.turbulent_prandtl': 1.5})
    cases['rough-smooth'] = ('gray_pilot', 6, 20, 450, 5, 1000, 5.0, 'all',
                             {'surface.land_roughness_momentum_m': 0.01,
                              'surface.land_roughness_heat_m': 0.001})
    # Geography and rotation.
    cases['land_ocean'] = ('gray_land_ocean', 6, 20, 450, 5, 1000, 5.0, 'all', {})
    cases['tidally_locked'] = ('gray_tidally_locked', 6, 20, 450, 5, 1000, 5.0, 'all', {})
    return cases


def derived_config(source, overrides, destination):
    """Write `source` with `overrides` applied, keeping every other line as shipped."""
    lines = source.read_text().splitlines()
    remaining = dict(overrides)
    output = []
    for line in lines:
        key = line.split('=', 1)[0].strip() if '=' in line else None
        if key in remaining:
            output.append(f'{key} = {remaining.pop(key)}')
        else:
            output.append(line)
    for key, value in remaining.items():
        output.append(f'{key} = {value}')
    destination.write_text('\n'.join(output) + '\n')
    return destination


def read_snapshot(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def weighted_mean_temperature(rows):
    weight = sum(float(r['area_m2']) * float(r['air_mass_kg_m2']) for r in rows)
    return sum(float(r['area_m2']) * float(r['air_mass_kg_m2']) * float(r['temperature_k'])
               for r in rows) / weight


def temperature_rms(first, second):
    """Mass- and area-weighted temperature RMS; only defined on identical grids."""
    if len(first) != len(second):
        return None
    weight = 0.0
    total = 0.0
    for a, b in zip(first, second):
        if a['cell'] != b['cell'] or a['level'] != b['level']:
            return None
        w = float(a['area_m2']) * float(a['air_mass_kg_m2'])
        weight += w
        total += w * (float(a['temperature_k']) - float(b['temperature_k'])) ** 2
    return math.sqrt(total / weight)


def interval_means(path, area):
    """Accumulated-flux means over the accepted steps, not a mean of sampled rates."""
    seconds = 0.0
    toa = 0.0
    height = 0.0
    with path.open() as stream:
        for row in csv.DictReader(stream):
            dt = float(row['dt_s'])
            seconds += dt
            toa += float(row['toa_net_j'])
            height += dt * float(row['mean_h_m'])
    if seconds <= 0.0:
        return {}
    return {'toa_net_upward_w_m2': toa / (seconds * area),
            'mean_boundary_layer_height_m': height / seconds,
            'accepted_seconds': seconds}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path,
                        default=Path('build/release/tools/benchmark_dry_mixing'))
    parser.add_argument('--configs', type=Path, default=Path('configs'))
    parser.add_argument('--output', type=Path, default=Path('output/phase12-comparison'))
    parser.add_argument('--days', type=float, default=1.0)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    results = {}
    failures = 0
    for name, (preset, n, k, dt, iterations, top, refinement, leg, overrides) in \
            build_cases().items():
        source = args.configs / f'phase12_{preset}.cfg'
        config = source
        if overrides:
            config = derived_config(source, overrides, args.output / f'{name}.cfg')
        snapshot = args.output / f'{name}.csv'
        command = [str(args.binary), str(config), str(args.days), str(n), str(k),
                   str(dt), str(iterations), str(top), str(refinement), leg,
                   str(snapshot)]
        print(name, flush=True)
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        (args.output / f'{name}.log').write_text(completed.stdout + completed.stderr)
        result = {'command': command, 'returncode': completed.returncode,
                  'overrides': overrides}
        if completed.returncode == 0:
            result['metrics'] = {key: float(value) for key, value in
                                 (line.split('=', 1)
                                  for line in completed.stdout.splitlines())}
            rows = read_snapshot(snapshot)
            area = sum(float(r['area_m2']) for r in rows) / k
            result['metrics']['mean_temperature_k'] = weighted_mean_temperature(rows)
            result['metrics'].update(
                interval_means(Path(str(snapshot) + '.steps.csv'), area))
        else:
            result['error'] = completed.stderr
            failures += 1
        results[name] = result

    # RMS comparisons are only meaningful between cases that share a grid.
    comparisons = {}
    pairs = [('time-1800-i3', 'time-450-i5'), ('time-1800-i5', 'time-450-i5'),
             ('time-900-i5', 'time-450-i5'), ('time-450-i5', 'explicit'),
             ('matrix-radiation', 'matrix-all'), ('matrix-convection', 'matrix-all'),
             ('matrix-boundary_layer', 'matrix-all'),
             ('uniform-k20', 'matrix-all'), ('richardson-0.25', 'matrix-all'),
             ('richardson-0.5', 'matrix-all'), ('gustiness-0', 'matrix-all'),
             ('gustiness-0.5', 'matrix-all'), ('prandtl-1.5', 'matrix-all'),
             ('rough-smooth', 'matrix-all'), ('top-500', 'matrix-all'),
             ('top-250', 'matrix-all')]
    for first, second in pairs:
        if results.get(first, {}).get('returncode') != 0:
            continue
        if results.get(second, {}).get('returncode') != 0:
            continue
        rms = temperature_rms(read_snapshot(args.output / f'{first}.csv'),
                              read_snapshot(args.output / f'{second}.csv'))
        entry = {'temperature_rms_k': rms}
        for key in ('toa_net_upward_w_m2', 'mean_boundary_layer_height_m'):
            a = results[first]['metrics'].get(key)
            b = results[second]['metrics'].get(key)
            if a is not None and b is not None:
                entry[f'{key}_difference'] = a - b
                if b != 0.0:
                    entry[f'{key}_relative'] = (a - b) / abs(b)
        comparisons[f'{first}_vs_{second}'] = entry

    (args.output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    (args.output / 'comparisons.json').write_text(
        json.dumps(comparisons, indent=2) + '\n')
    print(f'{len(results) - failures}/{len(results)} cases completed')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
