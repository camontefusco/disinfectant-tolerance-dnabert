#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import TypeVar

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.io import read_manifest
from disinfectant_tolerance.windows import iter_windows


T = TypeVar("T")


def evenly_spaced(items: list[T], maximum_items: int) -> list[T]:
    if maximum_items <= 0:
        raise ValueError("maximum_items must be positive")
    if len(items) <= maximum_items:
        return items
    if maximum_items == 1:
        return [items[0]]
    return [
        items[round(index * (len(items) - 1) / (maximum_items - 1))]
        for index in range(maximum_items)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--maximum-windows-per-isolate", type=int, default=16)
    parser.add_argument("--window-size", type=int, default=999)
    parser.add_argument("--window-step", type=int, default=999)
    parser.add_argument("--minimum-tail-size", type=int, default=500)
    args = parser.parse_args()

    rows = read_manifest(args.manifest)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["isolate_id", "contig_id", "start", "end", "sequence", "label", "group"]
    count = 0
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            windows = list(
                iter_windows(
                    isolate_id=row["isolate_id"],
                    assembly_path=row["assembly_path"],
                    label=int(row["label"]),
                    group=row["group"],
                    window_size=args.window_size,
                    step=args.window_step,
                    minimum_tail_size=args.minimum_tail_size,
                )
            )
            for window in evenly_spaced(windows, args.maximum_windows_per_isolate):
                writer.writerow(window.__dict__)
                count += 1
    print(f"Wrote {count} sampled windows to {output}")


if __name__ == "__main__":
    main()
