#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.io import read_manifest
from disinfectant_tolerance.windows import iter_windows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--output", default="data/processed/windows.csv")
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    rows = read_manifest(args.manifest)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "isolate_id",
        "contig_id",
        "start",
        "end",
        "sequence",
        "label",
        "group",
    ]
    count = 0
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            for window in iter_windows(
                isolate_id=row["isolate_id"],
                assembly_path=row["assembly_path"],
                label=int(row["label"]),
                group=row["group"],
                window_size=config["window_size"],
                step=config["window_step"],
                minimum_tail_size=config["minimum_tail_size"],
            ):
                writer.writerow(window.__dict__)
                count += 1
    print(f"Wrote {count} windows to {output}")


if __name__ == "__main__":
    main()

