#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
import shutil
from urllib.request import urlopen


def download_fasta(url: str, output: Path) -> None:
    if output.exists():
        print(f"Exists: {output}")
        return
    print(f"Downloading: {url}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url) as response, gzip.GzipFile(fileobj=response) as archive:
        with output.open("wb") as handle:
            shutil.copyfileobj(archive, handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True)
    parser.add_argument("--output-dir", default="data/raw/assemblies")
    parser.add_argument("--limit-isolates", type=int)
    parser.add_argument("--minimum-free-gb", type=float, default=15.0)
    args = parser.parse_args()

    with Path(args.plan).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if args.limit_isolates is not None:
        rows = rows[: args.limit_isolates]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    compressed_bytes = sum(int(row.get("submitted_bytes", 0) or 0) for row in rows)
    estimated_output_bytes = compressed_bytes * 5
    free_bytes = shutil.disk_usage(output_dir).free
    minimum_free_bytes = int(args.minimum_free_gb * 1_000_000_000)
    if estimated_output_bytes + minimum_free_bytes > free_bytes:
        raise RuntimeError(
            "Insufficient disk space: estimated FASTA output requires "
            f"{estimated_output_bytes / 1_000_000_000:.2f} GB and the safety "
            f"margin is {args.minimum_free_gb:.2f} GB, but only "
            f"{free_bytes / 1_000_000_000:.2f} GB are free."
        )
    print(
        f"Estimated FASTA output: {estimated_output_bytes / 1_000_000_000:.2f} GB; "
        f"free disk space: {free_bytes / 1_000_000_000:.2f} GB"
    )

    for row in rows:
        output = output_dir / f"{row['isolate_id']}.fasta"
        download_fasta(row["assembly_url"], output)


if __name__ == "__main__":
    main()

