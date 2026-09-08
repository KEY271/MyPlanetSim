#!/usr/bin/env python3
"""Collect the small, reproducible Phase 13 matrix summaries for the report."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def successful_metrics(path):
    data = json.loads(path.read_text())
    return {name: entry["metrics"] for name, entry in data.items()
            if entry["returncode"] == 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path,
                        default=ROOT / "output/phase13-comparison")
    parser.add_argument("--destination", type=Path,
                        default=ROOT / "docs/validation/phase-13-results.json")
    args = parser.parse_args()
    result = {}
    for suite in ("one-day", "thirty-day"):
        directory = args.comparison / suite
        result[suite.replace("-", "_")] = successful_metrics(
            directory / "results.json")
        result[f"{suite.replace('-', '_')}_comparisons"] = json.loads(
            (directory / "comparisons.json").read_text())
    args.destination.write_text(json.dumps(result, indent=2) + "\n")
    print(args.destination)


if __name__ == "__main__":
    raise SystemExit(main())
