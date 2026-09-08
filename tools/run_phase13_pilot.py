#!/usr/bin/env python3
"""Run the Phase 13 pilot in foreground segments with restart-safe outputs.

This command intentionally stays attached to the terminal.  Ctrl-C is forwarded to
the simulator, which writes a checkpoint; use --resume to continue it later.
"""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
import signal
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def config_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def checkpoint_position(path):
    values = config_values(path)
    return float(values["time_s"]), int(values["step"])


def next_segment(output):
    segments = output / "segments"
    segments.mkdir(parents=True, exist_ok=True)
    indexes = []
    for path in segments.glob("segment-*"):
        try:
            indexes.append(int(path.name.split("-", 1)[1]))
        except ValueError:
            pass
    index = max(indexes, default=0) + 1
    destination = segments / f"segment-{index:04d}"
    destination.mkdir()
    return destination


def archive_outputs(output, destination):
    for path in output.iterdir():
        if not path.is_file():
            continue
        if path.suffix == ".csv" or path.name == "run-metadata.txt":
            shutil.move(str(path), destination / path.name)


def run_process(command, log_path):
    print("+ " + " ".join(str(item) for item in command), flush=True)
    interrupted = False
    with log_path.open("w") as log:
        child = subprocess.Popen(
            command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, start_new_session=True)
        try:
            for line in child.stdout:
                print(line, end="", flush=True)
                log.write(line)
        except KeyboardInterrupt:
            interrupted = True
            print("\nCtrl-C: requesting a graceful checkpoint...", flush=True)
            child.send_signal(signal.SIGINT)
            for line in child.stdout:
                print(line, end="", flush=True)
                log.write(line)
        return child.wait(), interrupted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/release/my_planet_sim")
    parser.add_argument("--config", type=Path,
                        default=ROOT / "configs/phase13_moist_pilot.cfg")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--segment-days", type=float, default=100.0)
    parser.add_argument("--progress-interval-s", type=float, default=10.0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the first simulator command without running it")
    args = parser.parse_args()

    config = args.config.resolve()
    binary = args.binary.resolve()
    values = config_values(config)
    output_value = Path(values["output.directory"])
    output = output_value if output_value.is_absolute() else ROOT / output_value
    checkpoint = (args.checkpoint.resolve() if args.checkpoint else
                  output / "pilot.chk")
    end_time = float(values["run.end_time_s"])
    start_time = float(values["run.start_time_s"])
    time_step = float(values["run.time_step_s"])
    segment_steps = round(args.segment_days * 86400.0 / time_step)
    if args.segment_days <= 0 or segment_steps < 1:
        parser.error("--segment-days must contain at least one configured step")
    if not binary.is_file():
        parser.error(f"binary does not exist: {binary}")
    if args.dry_run:
        command = [str(binary), "--config", str(config),
                   "--progress-interval-s", str(args.progress_interval_s),
                   "--stop-after-step", str(segment_steps),
                   "--checkpoint", str(checkpoint)]
        if args.resume:
            command.extend(("--restart", str(checkpoint)))
        print(" ".join(command))
        return 0

    output.mkdir(parents=True, exist_ok=True)

    existing_segments = list((output / "segments").glob("segment-*")) \
        if (output / "segments").exists() else []
    if args.resume:
        if not checkpoint.is_file():
            parser.error(f"--resume requires checkpoint: {checkpoint}")
        current_time, current_step = checkpoint_position(checkpoint)
    else:
        if checkpoint.exists() or existing_segments:
            parser.error("pilot output already exists; use --resume or a new output")
        current_time, current_step = start_time, 0

    started = datetime.now(timezone.utc).isoformat()
    print(f"Phase 13 pilot started {started}; target={end_time / 86400:g} days",
          flush=True)
    interrupted = False
    while current_time < end_time:
        destination = next_segment(output)
        archive_outputs(output, destination)
        stop_step = current_step + segment_steps
        command = [str(binary), "--config", str(config),
                   "--progress-interval-s", str(args.progress_interval_s),
                   "--stop-after-step", str(stop_step),
                   "--checkpoint", str(checkpoint)]
        if checkpoint.is_file():
            command.extend(("--restart", str(checkpoint)))
        returncode, interrupted = run_process(command, destination / "run.log")
        archive_outputs(output, destination)
        if returncode != 0:
            print(f"simulator exited with status {returncode}", file=sys.stderr)
            return returncode
        if not checkpoint.is_file():
            print("simulator did not write the requested checkpoint", file=sys.stderr)
            return 1
        current_time, current_step = checkpoint_position(checkpoint)
        print(f"checkpoint: day={current_time / 86400:.6f} step={current_step}",
              flush=True)
        if interrupted:
            break

    summary_command = [sys.executable, str(ROOT / "tools/summarize_moist_pilot.py"),
                       "--output", str(output), "--config", str(config)]
    summary = subprocess.run(summary_command, cwd=ROOT, check=False)
    if summary.returncode != 0:
        return summary.returncode
    if interrupted:
        print("Pilot stopped cleanly. Resume with the same command plus --resume.")
        return 130
    print("Phase 13 pilot completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
