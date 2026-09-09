#!/usr/bin/env python3
"""Compare the Phase 14 integrated candidate against the Phase 13 legacy baseline.

Both points use the same binary, grid, vertical coordinate, tracer registry and
model interval; only the physics schedule (and optionally the ICI iteration count)
differs.  The registered gates of `docs/validation/phase-14-manifest.json` are
evaluated against the baseline of the same run length, so a one-day comparison uses
the one-day row and a thirty-day comparison uses the thirty-day row.

Short comparisons repeat every point through benchmark_phase14.py.  Long
comparisons (30 or 1200 model days) should use `--single`, which runs each point
once in the foreground; this measures the solution difference, not the run-to-run
timing spread.  Nothing here tunes coefficients, tolerances or physical parameters.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BASELINE_CONFIG = ROOT / 'configs/phase13_moist_pilot.cfg'
CANDIDATE_CONFIG = ROOT / 'configs/phase14_moist_candidate.cfg'


def metrics(text):
    result = {}
    for line in text.splitlines():
        if '=' not in line:
            continue
        key, value = (part.strip() for part in line.split('=', 1))
        try:
            number = float(value)
        except ValueError:
            continue
        result[key] = number if math.isfinite(number) else None
    if not isinstance(result.get('elapsed_s'), float):
        raise ValueError('benchmark did not emit finite elapsed_s')
    return result


def capture(command):
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    return result.stdout.strip() if result.returncode == 0 else None


def snapshot_rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def check_same_grid(rows, references):
    if not rows or len(rows) != len(references):
        raise ValueError('snapshot dimensions differ')
    for row, reference in zip(rows, references):
        if (row['cell'], row['level'], row['area_m2']) != (
                reference['cell'], reference['level'], reference['area_m2']):
            raise ValueError('snapshot grids differ')


def mass_weighted_temperature_rms(rows, references):
    check_same_grid(rows, references)
    weighted = weight = 0.0
    for row, reference in zip(rows, references):
        mass = float(reference['area_m2']) * float(reference['air_mass_kg_m2'])
        delta = float(row['temperature_k']) - float(reference['temperature_k'])
        if not math.isfinite(mass) or not math.isfinite(delta) or mass <= 0:
            raise ValueError('invalid snapshot value')
        weighted += mass * delta**2
        weight += mass
    return math.sqrt(weighted / weight)


def zonal_means(rows, band_deg, field):
    """Area-weighted band means per (latitude band, level)."""
    totals = {}
    for row in rows:
        band = math.floor(float(row['latitude_deg']) / band_deg)
        key = (band, row['level'])
        area = float(row['area_m2'])
        value = float(row[field])
        if not math.isfinite(area) or not math.isfinite(value) or area <= 0:
            raise ValueError('invalid snapshot value')
        total = totals.setdefault(key, [0.0, 0.0])
        total[0] += area * value
        total[1] += area
    return {key: value / area for key, (value, area) in totals.items()}


def zonal_rms(rows, references, band_deg, field):
    check_same_grid(rows, references)
    candidate = zonal_means(rows, band_deg, field)
    baseline = zonal_means(references, band_deg, field)
    if candidate.keys() != baseline.keys():
        raise ValueError('zonal bands differ')
    return math.sqrt(sum((candidate[key] - baseline[key])**2
                         for key in candidate) / len(candidate))


def block_statistics(samples_path, first_day, last_day):
    """Time means and interval-integrated exchanges over (first_day, last_day].

    The day bounds are offsets from the first sample, so a run that starts from a
    developed checkpoint uses the same block as a run that starts at day zero.
    Instantaneous fields of a chaotic run separate even when the climate agrees, so
    the adoption comparison also needs averages over a block of days. Cumulative
    ledgers are differenced across the block, never sampled at one instant.
    """
    rows = snapshot_rows(samples_path)
    if not rows:
        raise ValueError(f'{samples_path}: no samples')
    origin = float(rows[0]['day'])
    first_day += origin
    last_day += origin
    block = [row for row in rows if first_day < float(row['day']) <= last_day]
    start = [row for row in rows if float(row['day']) <= first_day]
    if not block or not start:
        raise ValueError(f'{samples_path}: block ({first_day}, {last_day}] is empty')
    opening = start[-1]
    closing = block[-1]
    seconds = (float(closing['day']) - float(opening['day'])) * 86400.0
    if seconds <= 0:
        raise ValueError(f'{samples_path}: block has no duration')
    result = {f'mean_{field}': sum(float(row[field]) for row in block) / len(block)
              for field in ('mean_temperature_k', 'mean_surface_temperature_k',
                            'precipitable_water_kg_m2')}
    for field, name in (('cumulative_evaporation_kg', 'evaporation_kg_s'),
                        ('cumulative_runoff_kg', 'runoff_kg_s')):
        result[name] = (float(closing[field]) - float(opening[field])) / seconds
    result['precipitation_kg_s'] = (
        (float(closing['cumulative_convective_precipitation_kg']) +
         float(closing['cumulative_grid_scale_precipitation_kg'])) -
        (float(opening['cumulative_convective_precipitation_kg']) +
         float(opening['cumulative_grid_scale_precipitation_kg']))) / seconds
    result['convective_fraction'] = (
        (float(closing['cumulative_convective_precipitation_kg']) -
         float(opening['cumulative_convective_precipitation_kg'])) /
        max(result['precipitation_kg_s'] * seconds, 1e-30))
    result['block_days'] = float(closing['day']) - float(opening['day'])
    result['block_samples'] = float(len(block))
    return result


def global_mean_temperature(rows):
    weighted = weight = 0.0
    for row in rows:
        mass = float(row['area_m2']) * float(row['air_mass_kg_m2'])
        weighted += mass * float(row['temperature_k'])
        weight += mass
    return weighted / weight


def measure(command, destination, repetitions, environment):
    if repetitions is None:
        completed = subprocess.run(command, cwd=ROOT, env=environment,
                                   text=True, capture_output=True, check=True)
        values = metrics(completed.stdout)
        destination.write_text(json.dumps(dict(
            schema_version=1, mode='single', command=command,
            time_utc=datetime.now(timezone.utc).isoformat(),
            commit=capture(['git', 'rev-parse', 'HEAD']),
            dirty=bool(capture(['git', 'status', '--porcelain'])),
            platform=platform.platform(),
            cpu_model=capture(['sysctl', '-n', 'machdep.cpu.brand_string']),
            logical_cpu_count=os.cpu_count(),
            conditions='single foreground run; power/CPU affinity not controlled',
            binary_and_config_sha256={
                argument: hashlib.sha256(Path(argument).read_bytes()).hexdigest()
                for argument in command
                if Path(argument).is_file() and (Path(argument).suffix in
                                                 ('.cfg', '.chk') or
                                                 argument == command[0])},
            summary={key: dict(median=value, minimum=value, maximum=value)
                     for key, value in values.items()}), indent=2) + '\n')
        return
    subprocess.run([sys.executable, str(ROOT / 'tools/benchmark_phase14.py'),
                    '--repetitions', str(repetitions), '--output', str(destination),
                    '--', *command], cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path,
                        default=ROOT / 'build/release/tools/benchmark_moist')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--candidate-config', type=Path, default=CANDIDATE_CONFIG,
                        help='process-interval preset to compare against the baseline')
    parser.add_argument('--days', type=float, default=1.0)
    parser.add_argument('--n', type=int, default=6)
    parser.add_argument('--levels', type=int, default=20)
    parser.add_argument('--dt', type=float, default=1800.0)
    parser.add_argument('--iterations', type=int, default=3)
    parser.add_argument('--top-pa', type=float, default=20000.0)
    parser.add_argument('--refinement', type=float, default=5.0)
    parser.add_argument('--band-deg', type=float, default=10.0,
                        help='latitude band width of the zonal-mean gates')
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--single', action='store_true',
                        help='run each point once instead of repeating it')
    parser.add_argument('--with-ici2', action='store_true',
                        help='also measure the candidate with two ICI iterations')
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--source-config', type=Path)
    parser.add_argument('--members', type=int, default=0,
                        help='number of micro-perturbed baseline members used to '
                             'measure how much of a difference is internal variability')
    parser.add_argument('--perturbation-k', type=float, default=1e-10,
                        help='amplitude of the deterministic initial temperature '
                             'perturbation of each member')
    parser.add_argument('--statistics-from-day', type=float,
                        help='first day excluded from the block statistics; the '
                             'default discards the first half of the run')
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    if bool(args.checkpoint) != bool(args.source_config):
        parser.error('checkpoint and source-config must be supplied together')
    if not args.single and args.repetitions < 5:
        parser.error('repeated measurements need at least five repetitions')
    repetitions = None if args.single else args.repetitions
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / 'docs/validation/phase-14-manifest.json'
    manifest = json.loads(manifest_path.read_text())

    points = {'baseline': (BASELINE_CONFIG, args.iterations),
              'candidate': (args.candidate_config, args.iterations)}
    if args.with_ici2:
        points['candidate-ici2'] = (args.candidate_config, 2)
    # Members repeat the baseline method; only the starting state differs, by an
    # amplitude far below any physical or discretization difference.
    for member in range(1, args.members + 1):
        points[f'member{member}'] = (BASELINE_CONFIG, args.iterations)
    measurements = {}
    provenance = {}
    environment = {key: value for key, value in os.environ.items()
                   if key not in ('MPS_COUNT_ALLOCATIONS', 'MPS_PROFILE_RHS')}
    for key, (config, iterations) in points.items():
        destination = args.output / key
        command = [str(args.binary.resolve()), str(config.resolve()),
                   repr(args.days), str(args.n), str(args.levels), repr(args.dt),
                   str(iterations), repr(args.top_pa), repr(args.refinement),
                   str(destination.with_suffix('.csv').resolve())]
        if args.checkpoint:
            command += [str(args.checkpoint.resolve()),
                        str(args.source_config.resolve())]
        member_environment = dict(environment)
        if key.startswith('member'):
            member_environment['MPS_INITIAL_PERTURBATION_K'] = repr(
                args.perturbation_k * int(key.removeprefix('member')))
        else:
            member_environment.pop('MPS_INITIAL_PERTURBATION_K', None)
        if not args.summarize_only:
            print(f'Measuring {key}', flush=True)
            measure(command, destination.with_suffix('.json'), repetitions,
                    member_environment)
        result = json.loads(destination.with_suffix('.json').read_text())
        if result['command'] != command:
            raise ValueError(f'{key}: stored measurement command differs')
        measurements[key] = {name: value['median']
                             for name, value in result['summary'].items()}
        provenance[key] = {name: result.get(name) for name in (
            'command', 'mode', 'commit', 'dirty', 'platform', 'cpu_model',
            'conditions', 'binary_and_config_sha256')}
        provenance[key]['snapshot_sha256'] = hashlib.sha256(
            destination.with_suffix('.csv').read_bytes()).hexdigest()

    reference = measurements['baseline']
    reference_rows = snapshot_rows(args.output / 'baseline.csv')
    floor = manifest['conservation']['rain_evaporation_absolute_floor_kg_m2_s']
    gate = manifest['one_day'] if args.days <= 1 else manifest['thirty_day']
    gate_name = 'one_day' if args.days <= 1 else 'thirty_day'
    if 1 < args.days < 30:
        gate_name = 'thirty_day (applied to a shorter interval)'
    first_day = (args.statistics_from_day if args.statistics_from_day is not None
                 else args.days / 2.0)
    reference_block = block_statistics(args.output / 'baseline.csv.samples.csv',
                                       first_day, args.days)
    comparisons = {}
    for key, values in measurements.items():
        rows = snapshot_rows(args.output / f'{key}.csv')
        comparison = {
            'temperature_rms_k': mass_weighted_temperature_rms(rows, reference_rows),
            'mean_temperature_difference_k': abs(global_mean_temperature(rows) -
                                                 global_mean_temperature(reference_rows)),
            'zonal_temperature_rms_k': zonal_rms(rows, reference_rows, args.band_deg,
                                                 'temperature_k'),
            'zonal_wind_rms_m_s': zonal_rms(rows, reference_rows, args.band_deg,
                                            'zonal_wind_m_s'),
            'toa_difference_w_m2': abs(values['toa_net_upward_w_m2'] -
                                       reference['toa_net_upward_w_m2']),
            'convective_fraction_difference': abs(
                values['convective_precipitation_fraction'] -
                reference['convective_precipitation_fraction']),
            'speedup_over_baseline': reference['seconds_per_model_day'] /
                                     values['seconds_per_model_day'],
        }
        for name in ('precipitable_water_kg_m2', 'evaporation_kg_m2_s',
                     'precipitation_kg_m2_s'):
            comparison[name + '_relative_difference'] = (
                abs(values[name] - reference[name]) /
                max(abs(reference[name]), floor if name.endswith('_s') else 1e-12))
        for name in ('boundary_layer_column_calls', 'convection_column_calls',
                     'sbm_diagnostic_column_calls', 'full_rhs_per_step'):
            if reference.get(name):
                comparison[name + '_ratio'] = reference[name] / values[name] \
                    if values.get(name) else None
        conservation = manifest['conservation']
        comparison['conservation_pass'] = (
            values['retries'] == 0 and values['physics_retries'] == 0 and
            values['minimum_vapor_mixing_ratio'] >= 0 and
            abs(values['water_budget_residual_ratio']) <=
            conservation['water_and_enthalpy_scaled_residual'] and
            abs(values['dry_mass_relative_drift']) <=
            conservation['passive_mass_relative'] and
            values['maximum_passive_tracer_relative_drift'] <=
            conservation['passive_mass_relative'])
        checks = {
            'pw_relative': comparison['precipitable_water_kg_m2_relative_difference']
                           <= gate['pw_relative'],
            'evaporation_precipitation': all(
                abs(values[name] - reference[name]) <= max(
                    floor, gate['integrated_evaporation_precipitation_relative'] *
                    abs(reference[name]))
                for name in ('evaporation_kg_m2_s', 'precipitation_kg_m2_s')),
            'toa': comparison['toa_difference_w_m2'] <= gate['toa_w_m2'],
            'convective_fraction': comparison['convective_fraction_difference'] <=
                                   manifest['moist_distribution']['convective_fraction_absolute'],
            'conservation': comparison['conservation_pass'],
        }
        if gate_name == 'one_day':
            checks['temperature_rms'] = (comparison['temperature_rms_k'] <=
                                         gate['mass_weighted_temperature_rms_k'])
        else:
            checks['mean_temperature'] = (comparison['mean_temperature_difference_k']
                                          <= gate['mean_temperature_k'])
            checks['zonal_temperature_rms'] = (comparison['zonal_temperature_rms_k'] <=
                                               gate['zonal_temperature_rms_k'])
            checks['zonal_wind_rms'] = (comparison['zonal_wind_rms_m_s'] <=
                                        gate['zonal_wind_rms_m_s'])
        block = block_statistics(args.output / f'{key}.csv.samples.csv',
                                 first_day, args.days)
        comparison['block'] = block
        comparison['block_mean_temperature_difference_k'] = abs(
            block['mean_mean_temperature_k'] - reference_block['mean_mean_temperature_k'])
        comparison['block_surface_temperature_difference_k'] = abs(
            block['mean_mean_surface_temperature_k'] -
            reference_block['mean_mean_surface_temperature_k'])
        for field in ('mean_precipitable_water_kg_m2', 'evaporation_kg_s',
                      'precipitation_kg_s'):
            comparison['block_' + field + '_relative_difference'] = abs(
                block[field] - reference_block[field]) / max(
                    abs(reference_block[field]), 1e-30)
        checks['block_mean_temperature'] = (
            comparison['block_mean_temperature_difference_k'] <=
            gate['mean_temperature_k'])
        checks['block_pw'] = (
            comparison['block_mean_precipitable_water_kg_m2_relative_difference'] <=
            gate['pw_relative'])
        checks['block_evaporation_precipitation'] = all(
            comparison['block_' + field + '_relative_difference'] <=
            gate['integrated_evaporation_precipitation_relative']
            for field in ('evaporation_kg_s', 'precipitation_kg_s'))
        comparison['gate_checks'] = checks
        comparison['measured_gates_pass'] = all(checks.values())
        comparisons[key] = comparison

    # The noise floor of this comparison: the largest difference that the unchanged
    # baseline method produces from a perturbation of 1e-10 K.
    internal_spread = None
    member_keys = [key for key in comparisons if key.startswith('member')]
    if member_keys:
        internal_spread = {
            field: max(comparisons[key][field] for key in member_keys)
            for field in ('temperature_rms_k', 'mean_temperature_difference_k',
                          'zonal_temperature_rms_k', 'zonal_wind_rms_m_s',
                          'toa_difference_w_m2', 'convective_fraction_difference',
                          'precipitable_water_kg_m2_relative_difference',
                          'evaporation_kg_m2_s_relative_difference',
                          'precipitation_kg_m2_s_relative_difference',
                          'block_mean_temperature_difference_k',
                          'block_surface_temperature_difference_k',
                          'block_mean_precipitable_water_kg_m2_relative_difference',
                          'block_evaporation_kg_s_relative_difference',
                          'block_precipitation_kg_s_relative_difference')}
        internal_spread['members'] = len(member_keys)
        internal_spread['perturbation_k'] = args.perturbation_k
        internal_spread['candidate_exceeds_spread'] = {
            field: comparisons['candidate'][field] > value
            for field, value in internal_spread.items()
            if isinstance(value, float) and field != 'perturbation_k'}

    result = dict(
        schema_version=1,
        reference='baseline = configs/phase13_moist_pilot.cfg legacy 300 s schedule',
        candidate=str(args.candidate_config.relative_to(ROOT)
                      if args.candidate_config.is_absolute() else args.candidate_config),
        days=args.days, grid=dict(cells_per_panel=args.n, levels=args.levels,
                                  dt_s=args.dt, iterations=args.iterations),
        applied_gate=gate_name, band_deg=args.band_deg,
        statistics_block=f'({first_day}, {args.days}] days',
        mode='single' if args.single else f'{args.repetitions} repetitions',
        manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        missing=['CAPE distribution', 'rain p95/p99 tails', 'top-reaching fraction',
                 'terrain and finite land fraction cases'],
        internal_spread=internal_spread,
        provenance=provenance, points=measurements, comparisons=comparisons)
    (args.output / 'comparison.json').write_text(
        json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(args.output / 'comparison.json')
    for key, comparison in comparisons.items():
        print(f"{key}: speedup={comparison['speedup_over_baseline']:.4f} "
              f"gates={'pass' if comparison['measured_gates_pass'] else 'FAIL'}")
    if internal_spread is not None:
        print('internal spread of ' + str(internal_spread['members']) +
              f" perturbed baseline members (+-{args.perturbation_k} K):")
        for field, value in internal_spread.items():
            if isinstance(value, float) and field != 'perturbation_k':
                print(f"  {field}: candidate={comparisons['candidate'][field]:.6g} "
                      f"spread={value:.6g}")


if __name__ == '__main__':
    main()
