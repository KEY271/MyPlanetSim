#!/usr/bin/env python3
"""Repeat a benchmark with separate profiling/allocation runs and provenance.

Example: benchmark_phase14.py --output output/p14/baseline.json -- \
 build/release/tools/benchmark_moist configs/phase13_moist_pilot.cfg \
 1 6 20 1800 3 20000 5 output/p14/state.csv
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def capture(command):
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    return result.stdout.strip() if result.returncode == 0 else None


def metrics(text):
    result = {}
    for line in text.splitlines():
        if '=' not in line:
            continue
        key, value = (part.strip() for part in line.split('=', 1))
        try:
            number = float(value)
            result[key] = number if math.isfinite(number) else None
        except ValueError:
            result[key] = value
    if not isinstance(result.get('elapsed_s'), (int, float)):
        raise ValueError('benchmark did not emit finite elapsed_s')
    return result


def run(command, environment):
    completed = subprocess.run(command, cwd=ROOT, env=environment,
                               text=True, capture_output=True, check=True)
    return metrics(completed.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--conditions', default='power/CPU affinity not controlled')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or args.repetitions < 5:
        parser.error('supply a benchmark command and at least five repetitions')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    environment.pop('MPS_COUNT_ALLOCATIONS', None)
    environment.pop('MPS_PROFILE_RHS', None)
    # Process warm-up; the program also initializes its workspaces internally.
    run(command, environment)
    runs = [run(command, environment) for _ in range(args.repetitions)]
    profile_runs = [run(command, dict(environment, MPS_PROFILE_RHS='1'))
                    for _ in range(args.repetitions)]
    allocation = run(command, dict(environment, MPS_COUNT_ALLOCATIONS='1'))
    numeric = {key for key in runs[0]
               if all(isinstance(row.get(key), (int, float)) for row in runs)}
    summary = {key: dict(median=statistics.median(row[key] for row in runs),
                         minimum=min(row[key] for row in runs),
                         maximum=max(row[key] for row in runs))
               for key in sorted(numeric)}
    fingerprints = {}
    for arg in command:
        path = Path(arg)
        if path.is_file() and (path.suffix == '.cfg' or arg == command[0]):
            fingerprints[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    result = dict(schema_version=1, time_utc=datetime.now(timezone.utc).isoformat(),
                  command=command, commit=capture(['git', 'rev-parse', 'HEAD']),
                  dirty=bool(capture(['git', 'status', '--porcelain'])),
                  platform=platform.platform(), cpu=platform.processor(),
                  cpu_model=capture(['sysctl', '-n', 'machdep.cpu.brand_string']),
                  logical_cpu_count=os.cpu_count(), conditions=args.conditions,
                  environment={k: v for k, v in environment.items()
                               if k.startswith(('OMP_', 'KMP_'))},
                  binary_and_config_sha256=fingerprints,
                  compiler_cache=(Path(command[0]).parent.parent / 'CMakeCache.txt').read_text()
                    if (Path(command[0]).parent.parent / 'CMakeCache.txt').is_file() else None,
                  missing=['hardware bandwidth', 'cache misses', 'vectorization report',
                           'scratch bytes/thread', 'estimated read/write bytes'],
                  summary=summary, runs=runs, profile_runs=profile_runs,
                  profile_overhead_ratio=statistics.median(r['elapsed_s'] for r in profile_runs)
                    / summary['elapsed_s']['median'] - 1,
                  allocation_run=allocation)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(args.output)
    print(json.dumps(summary.get('seconds_per_model_day', summary['elapsed_s'])))


if __name__ == '__main__':
    main()
