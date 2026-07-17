#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shutil
from urllib.request import urlopen


def download(url: str, output_dir: Path) -> None:
    output = output_dir / Path(url).name
    if output.exists():
        print(f"Exists: {output}")
        return
    print(f"Downloading: {url}")
    with urlopen(url) as response, output.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True)
    parser.add_argument("--output-dir", default="data/raw/fastq")
    parser.add_argument("--limit-isolates", type=int)
    parser.add_argument("--minimum-free-gb", type=float, default=15.0)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with Path(args.plan).open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    isolate_ids: list[str] = []
    for row in rows:
        isolate_id = row["isolate_id"]
        if isolate_id not in isolate_ids:
            isolate_ids.append(isolate_id)
    if args.limit_isolates is not None:
        allowed = set(isolate_ids[: args.limit_isolates])
        rows = [row for row in rows if row["isolate_id"] in allowed]

    planned_bytes = sum(int(row.get("fastq_bytes", 0) or 0) for row in rows)
    free_bytes = shutil.disk_usage(output_dir).free
    minimum_free_bytes = int(args.minimum_free_gb * 1_000_000_000)
    if planned_bytes + minimum_free_bytes > free_bytes:
        raise RuntimeError(
            "Insufficient disk space: planned FASTQ files require "
            f"{planned_bytes / 1_000_000_000:.2f} GB and the safety margin is "
            f"{args.minimum_free_gb:.2f} GB, but only "
            f"{free_bytes / 1_000_000_000:.2f} GB are free."
        )
    print(
        f"Planned FASTQ download: {planned_bytes / 1_000_000_000:.2f} GB; "
        f"free disk space: {free_bytes / 1_000_000_000:.2f} GB"
    )

    for row in rows:
        for url in row.get("fastq_urls", "").split(";"):
            if url:
                download(url, output_dir)


if __name__ == "__main__":
    main()
