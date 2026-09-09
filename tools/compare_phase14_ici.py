#!/usr/bin/env python3
"""Measure the fixed ICI 2/3/5 and ARK2 IMEX by dt=450/900/1800 matrix, no tuning.

Each point runs warm-up, five timings, five profiles, and allocation accounting
through benchmark_phase14.py. Optional checkpoint import tests developed states.
The reference keeps the legacy physics maximum at 300 s (450 s dynamics uses
two 225 s substeps). Same-dt comparisons isolate the time integration method;
the 450 s reference is NOT the finer 75 s physics or long-term climate gate.
The ARK2 comparator solves every vertical mode and adds fast operator
applications, so its cost is compared against ICI 2 at the same dt.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def temperature_rms(candidate, reference):
    with Path(candidate).open() as stream:
        rows = list(csv.DictReader(stream))
    with Path(reference).open() as stream:
        refs = list(csv.DictReader(stream))
    if not rows or len(rows) != len(refs):
        raise ValueError('snapshot dimensions differ')
    weighted = weight = 0.0
    for row, ref in zip(rows, refs):
        if (row['cell'], row['level'], row['area_m2']) != (
                ref['cell'], ref['level'], ref['area_m2']):
            raise ValueError('snapshot grids differ')
        mass = float(ref['area_m2']) * float(ref['air_mass_kg_m2'])
        delta = float(row['temperature_k']) - float(ref['temperature_k'])
        if not math.isfinite(mass) or not math.isfinite(delta) or mass <= 0:
            raise ValueError('invalid snapshot value')
        weighted += mass * delta**2
        weight += mass
    return math.sqrt(weighted / weight)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'build/release/tools/benchmark_moist')
    parser.add_argument('--config', type=Path, default=ROOT / 'configs/phase13_moist_pilot.cfg')
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--source-config', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--summarize-only', action='store_true',
                        help='recheck and aggregate existing measurements without rerunning')
    parser.add_argument('--measure-dt', type=int, nargs='+', choices=(450, 900, 1800),
                        help='measure only these dynamics time steps and reuse the '
                             'stored measurements of the other rows')
    args = parser.parse_args()
    if bool(args.checkpoint) != bool(args.source_config):
        parser.error('checkpoint and source-config must be supplied together')
    config_values = dict(line.split('#', 1)[0].strip().split('=', 1)
                         for line in args.config.read_text().splitlines()
                         if '=' in line.split('#', 1)[0])
    config_values = {k.strip(): v.strip() for k, v in config_values.items()}
    if float(config_values.get('moisture.maximum_physics_substep_s', '300')) != 300:
        parser.error('this matrix requires the legacy physics maximum of 300 s')
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / 'docs/validation/phase-14-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    points = {}
    provenance = {}
    for dt in (450, 900, 1800):
        for iterations in (2, 3, 5, -1):
            key = f'dt{dt}-' + ('ark2' if iterations < 0 else f'ici{iterations}')
            destination = args.output / key
            command = [str(args.binary.resolve()), str(args.config.resolve()),
                       '1', '6', '20', str(dt), str(iterations), '20000', '5',
                       str(destination.with_suffix('.csv').resolve())]
            if args.checkpoint:
                command += [str(args.checkpoint.resolve()), str(args.source_config.resolve())]
            measure = not args.summarize_only and (
                args.measure_dt is None or dt in args.measure_dt)
            if measure:
                print(f'Measuring {key}', flush=True)
                subprocess.run([sys.executable, str(ROOT / 'tools/benchmark_phase14.py'),
                                '--output', str(destination.with_suffix('.json').resolve()),
                                '--', *command], cwd=ROOT, check=True)
            if not measure and not destination.with_suffix('.json').is_file():
                # A staged partial measurement; aggregate once every row exists.
                print(f'Skipping unmeasured {key}', flush=True)
                continue
            result = json.loads(destination.with_suffix('.json').read_text())
            if result['command'] != command:
                raise ValueError(f'{key}: stored measurement command differs')
            points[key] = {k: v['median'] for k, v in result['summary'].items()}
            provenance[key] = {k: result[k] for k in (
                'command', 'commit', 'dirty', 'platform', 'cpu_model', 'conditions',
                'binary_and_config_sha256')}
            provenance[key]['elapsed_s'] = result['summary']['elapsed_s']
            provenance[key]['snapshot_sha256'] = hashlib.sha256(
                destination.with_suffix('.csv').read_bytes()).hexdigest()
            # Timers and allocation instrumentation must not alter the solution.
            for row in result['runs'] + result['profile_runs'] + [result['allocation_run']]:
                for field in ('mean_temperature_k', 'precipitable_water_kg_m2',
                              'evaporation_kg_m2_s', 'precipitation_kg_m2_s',
                              'full_rhs_per_step', 'accepted_steps', 'retries'):
                    if row[field] != points[key][field]:
                        raise ValueError(f'{key}: nondeterministic {field}')
    if len(points) != 12:
        print(f'{len(points)}/12 rows measured; comparison.json not written', flush=True)
        return
    for path in (args.binary, args.config, args.checkpoint, args.source_config):
        if path is not None and len({p['binary_and_config_sha256'][str(path.resolve())]
                                     for p in provenance.values()}) != 1:
            raise ValueError(f'matrix input changed during measurement: {path}')
    reference = points['dt450-ici5']
    comparisons = {}
    floor = manifest['conservation']['rain_evaporation_absolute_floor_kg_m2_s']
    for key, values in points.items():
        comparison = {'temperature_rms_k': temperature_rms(
            args.output / f'{key}.csv', args.output / 'dt450-ici5.csv')}
        comparison['same_dt_ici5_temperature_rms_k'] = temperature_rms(
            args.output / f'{key}.csv', args.output / (key.split('-')[0] + '-ici5.csv'))
        comparison['mean_physics_substep_s'] = (
            values['accepted_seconds'] / values['physics_substeps']
            if values['retries'] == 0 and values['physics_retries'] == 0 else None)
        for field in ('precipitable_water_kg_m2', 'evaporation_kg_m2_s',
                      'precipitation_kg_m2_s'):
            comparison[field + '_relative_difference'] = abs(values[field] - reference[field]) / max(
                abs(reference[field]), floor if field.endswith('_s') else 1e-12)
        comparison['toa_difference_w_m2'] = abs(values['toa_net_upward_w_m2'] - reference['toa_net_upward_w_m2'])
        comparison['speedup_over_ici3_same_dt'] = points[key.split('-')[0] + '-ici3']['elapsed_s'] / values['elapsed_s']
        # The IMEX adoption gate of the plan is measured against fixed two-iteration ICI.
        comparison['speedup_over_ici2_same_dt'] = points[key.split('-')[0] + '-ici2']['elapsed_s'] / values['elapsed_s']
        gate = manifest['one_day']
        comparison['measured_one_day_gates_pass'] = (
            comparison['temperature_rms_k'] <= gate['mass_weighted_temperature_rms_k']
            and comparison['precipitable_water_kg_m2_relative_difference'] <= gate['pw_relative']
            and all(abs(values[f] - reference[f]) <= max(
                        floor, gate['integrated_evaporation_precipitation_relative'] * abs(reference[f]))
                    for f in ('evaporation_kg_m2_s', 'precipitation_kg_m2_s'))
            and comparison['toa_difference_w_m2'] <= gate['toa_w_m2']
            and values['retries'] == 0 and values['physics_retries'] == 0
            and values['minimum_vapor_mixing_ratio'] >= 0
            and abs(values['water_budget_residual_ratio']) <= manifest['conservation']['water_and_enthalpy_scaled_residual']
            and abs(values['dry_mass_relative_drift']) <= manifest['conservation']['passive_mass_relative']
            and values['maximum_passive_tracer_relative_drift'] <= manifest['conservation']['passive_mass_relative'])
        comparisons[key] = comparison
    result = dict(schema_version=1, reference='dt450-ici5; legacy physics maximum 300 s, actual 225 s at dt450',
                  manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                  scope='N6/K20 one day; measured gates only, not final adoption',
                  missing=['moist enthalpy normalized residual', 'CAPE/rain tails',
                           'wave/terrain/baroclinic comparisons', '30/1200 day climate'],
                  provenance=provenance, points=points, comparisons=comparisons)
    (args.output / 'comparison.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
