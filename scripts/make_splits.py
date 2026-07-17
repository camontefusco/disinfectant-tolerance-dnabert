#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.io import read_manifest
from disinfectant_tolerance.splits import SplitFractions, assign_grouped_splits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--output", default="data/processed/isolate_splits.csv")
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    rows = read_manifest(args.manifest)
    assignments = assign_grouped_splits(
        rows,
        fractions=SplitFractions(**config["split_fractions"]),
        seed=config["random_seed"],
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) + ["split"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "split": assignments[row["group"]]})
    print(f"Wrote leakage-aware isolate splits to {output}")


if __name__ == "__main__":
    main()

