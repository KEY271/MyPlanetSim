#!/usr/bin/env python3
"""Summarize benchmark_dry_mixing's accepted-step and sampled climate diagnostics.

Accepted-step columns are integrated over the actual accepted seconds; the sampled
climate columns are integrated with the trapezoid rule over their own sample times. The
two are kept apart because the daily samples alias the diurnal cycle.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics

STEP_SUMS = ('cfl_retries', 'invariant_retries', 'solver_retries', 'toa_net_j',
             'surface_storage_j', 'legacy_sensible_j', 'legacy_drag_j',
             'boundary_layer_sensible_j', 'boundary_layer_surface_j',
             'returned_dissipation_j', 'convection_enthalpy_j',
             'convection_dry_energy_j')
STEP_MAXIMA = ('heat_residual_j', 'kinetic_identity_residual_j', 'momentum_residual_n_s',
               'tracer_change_kg', 'unstable_before', 'unstable_after', 'mean_h_m',
               'maximum_k_h_m2_s', 'shallow_fraction', 'model_top_fraction')
CLIMATE_KEYS = ('mean_temperature_k', 'mean_surface_temperature_k', 'maximum_wind_m_s',
                'toa_net_w_m2', 'surface_storage_w_m2', 'unstable_fraction',
                'minimum_dtheta_dz_k_m', 'mean_h_m')


def read_rows(path):
    with path.open() as stream:
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(stream)]
    if not rows or not all(math.isfinite(value) for row in rows
                           for value in row.values()):
        raise ValueError(f'empty or non-finite diagnostics: {path}')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('--radius-m', type=float, default=6371220.0,
                        help='planet radius used for flux normalization')
    args = parser.parse_args()
    if not math.isfinite(args.radius_m) or args.radius_m <= 0:
        parser.error('radius must be positive and finite')
    area = 4.0 * math.pi * args.radius_m**2
    steps = read_rows(Path(str(args.snapshot) + '.steps.csv'))
    climate = read_rows(Path(str(args.snapshot) + '.climate.csv'))
    final_day = steps[-1]['day']
    result = {'final_day': final_day, 'accepted_steps': len(steps), 'windows': {}}
    for name, start, end in [('first_100_days', 0, 100),
                             ('last_100_days', final_day - 100, final_day),
                             ('after_200_days', 200, final_day)]:
        selected = [row for row in steps if start < row['day'] <= end]
        samples = [row for row in climate if start <= row['day'] <= end]
        if not selected or not samples:
            continue
        window = {'accepted_steps': len(selected),
                  'median_dt_s': statistics.median(row['dt_s'] for row in selected),
                  'minimum_dt_s': min(row['dt_s'] for row in selected),
                  'maximum_dt_s': max(row['dt_s'] for row in selected),
                  'linear_iterations_per_step':
                      statistics.mean(row['linear_iterations'] for row in selected)}
        for key in STEP_SUMS:
            window[key] = sum(row[key] for row in selected)
        for key in STEP_MAXIMA:
            window[f'maximum_{key}'] = max(abs(row[key]) for row in selected)
        duration_s = sum(row['dt_s'] for row in selected)
        window['accepted_seconds'] = duration_s
        for key in ('toa_net_j', 'surface_storage_j', 'boundary_layer_sensible_j',
                    'returned_dissipation_j'):
            window[f'integrated_{key[:-2]}_w_m2'] = window[key] / (duration_s * area)
        for key in CLIMATE_KEYS:
            duration = samples[-1]['day'] - samples[0]['day']
            mean = (sum((b['day'] - a['day']) * (a[key] + b[key]) / 2
                        for a, b in zip(samples, samples[1:])) / duration
                    if duration > 0 else samples[0][key])
            window[key] = {'time_mean': mean,
                           'minimum': min(row[key] for row in samples),
                           'maximum': max(row[key] for row in samples)}
        result['windows'][name] = window
    result['maximum_sampled_mass_relative_drift'] = max(
        abs(row['mass_relative_drift']) for row in climate)
    result['final_unattributed_energy_j'] = climate[-1]['unattributed_energy_j']
    # The local process gates of the plan, normalized by the exchange they belong to.
    result['maximum_heat_residual_ratio'] = max(
        abs(row['heat_residual_j']) /
        max(area, abs(row['boundary_layer_sensible_j']) +
            abs(row['boundary_layer_surface_j']) + abs(row['returned_dissipation_j']))
        for row in steps)
    result['maximum_tracer_change_kg'] = max(abs(row['tracer_change_kg'])
                                             for row in steps)
    result['maximum_unstable_after'] = max(row['unstable_after'] for row in steps)
    destination = Path(str(args.snapshot) + '.summary.json')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
