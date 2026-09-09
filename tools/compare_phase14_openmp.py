#!/usr/bin/env python3
"""Run the Phase 14 serial/OpenMP scaling matrix in the foreground."""

import argparse
import json
import os
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('output/phase14/openmp'))
    parser.add_argument('--config', type=Path,
                        default=Path('configs/phase13_moist_pilot.cfg'))
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 2, 4])
    parser.add_argument('--days', type=float, default=1.0)
    parser.add_argument('--n', type=int, default=12)
    parser.add_argument('--levels', type=int, default=20)
    parser.add_argument('--dt', type=float, default=1800.0)
    parser.add_argument('--iterations', type=int, default=2)
    parser.add_argument('--top-pa', type=float, default=20000.0)
    parser.add_argument('--refinement', type=float, default=5.0)
    args = parser.parse_args()
    if any(thread < 1 for thread in args.threads):
        parser.error('thread counts must be positive')

    args.output.mkdir(parents=True, exist_ok=True)
    cases = [('serial', ROOT / 'build/release/tools/benchmark_moist', None)]
    cases.extend((f'openmp-{thread}',
                  ROOT / 'build/release-openmp/tools/benchmark_moist', thread)
                 for thread in args.threads)
    results = {}
    for name, binary, threads in cases:
        snapshot = args.output / f'{name}.csv'
        output = args.output / f'{name}.json'
        benchmark = [str(binary), str(args.config), str(args.days), str(args.n),
                     str(args.levels), str(args.dt), str(args.iterations),
                     str(args.top_pa), str(args.refinement), str(snapshot)]
        command = ['python3', str(ROOT / 'tools/benchmark_phase14.py'),
                   '--output', str(output), '--', *benchmark]
        environment = dict(os.environ)
        if threads is None:
            environment.pop('OMP_NUM_THREADS', None)
        else:
            environment['OMP_NUM_THREADS'] = str(threads)
            environment.setdefault('OMP_PROC_BIND', 'close')
            environment.setdefault('OMP_PLACES', 'cores')
        print(f'[{name}]', flush=True)
        subprocess.run(command, cwd=ROOT, env=environment, check=True)
        result = json.loads(output.read_text())
        results[name] = result['summary']['seconds_per_model_day']['median']

    serial = results['serial']
    summary = {
        name: {'seconds_per_model_day': elapsed,
               'speedup_over_serial': serial / elapsed}
        for name, elapsed in results.items()
    }
    destination = args.output / 'summary.json'
    destination.write_text(json.dumps(summary, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
