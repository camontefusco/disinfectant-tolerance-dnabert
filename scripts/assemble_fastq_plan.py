#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shlex
import shutil
import subprocess


def command_text(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def require_tool(*names: str) -> str:
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError(f"Missing required tool: {' or '.join(names)}")


def run(parts: list[str], execute: bool) -> None:
    print(command_text(parts))
    if execute:
        subprocess.run(parts, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True)
    parser.add_argument("--fastq-dir", default="data/raw/fastq")
    parser.add_argument("--assembly-dir", default="data/raw/assemblies")
    parser.add_argument("--work-dir", default="data/raw/assembly_work")
    parser.add_argument("--qc-dir", default="reports/assembly_qc")
    parser.add_argument("--limit-isolates", type=int)
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument("--cleanup-fastq", action="store_true")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    with Path(args.plan).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if args.limit_isolates is not None:
        allowed: list[str] = []
        for row in rows:
            if row["isolate_id"] not in allowed:
                allowed.append(row["isolate_id"])
        rows = [row for row in rows if row["isolate_id"] in allowed[: args.limit_isolates]]

    fastp = require_tool("fastp")
    spades = require_tool("spades.py", "spades")
    quast = require_tool("quast.py", "quast")
    fastq_dir = Path(args.fastq_dir)
    assembly_dir = Path(args.assembly_dir)
    work_dir = Path(args.work_dir)
    qc_dir = Path(args.qc_dir)

    for row in rows:
        isolate_id = row["isolate_id"]
        urls = [url for url in row["fastq_urls"].split(";") if url]
        if len(urls) != 2:
            raise ValueError(f"Expected paired FASTQ URLs for {isolate_id}, found {len(urls)}")
        raw_1, raw_2 = (fastq_dir / Path(url).name for url in urls)
        trimmed_dir = work_dir / isolate_id / "trimmed"
        spades_dir = work_dir / isolate_id / "spades"
        trimmed_1 = trimmed_dir / f"{isolate_id}_R1.fastq.gz"
        trimmed_2 = trimmed_dir / f"{isolate_id}_R2.fastq.gz"
        assembly_path = assembly_dir / f"{isolate_id}.fasta"

        if args.execute:
            missing = [path for path in (raw_1, raw_2) if not path.exists()]
            if missing:
                raise FileNotFoundError(f"Download FASTQ files before assembly: {missing}")
            trimmed_dir.mkdir(parents=True, exist_ok=True)
            assembly_dir.mkdir(parents=True, exist_ok=True)
            qc_dir.mkdir(parents=True, exist_ok=True)

        run(
            [
                fastp,
                "--in1",
                str(raw_1),
                "--in2",
                str(raw_2),
                "--out1",
                str(trimmed_1),
                "--out2",
                str(trimmed_2),
                "--json",
                str(trimmed_dir / "fastp.json"),
                "--html",
                str(trimmed_dir / "fastp.html"),
            ],
            args.execute,
        )
        run(
            [
                spades,
                "--careful",
                "-1",
                str(trimmed_1),
                "-2",
                str(trimmed_2),
                "-o",
                str(spades_dir),
            ],
            args.execute,
        )
        if args.execute:
            shutil.copyfile(spades_dir / "scaffolds.fasta", assembly_path)
        run([quast, str(assembly_path), "-o", str(qc_dir / isolate_id)], args.execute)
        if args.execute and not args.keep_work:
            shutil.rmtree(work_dir / isolate_id)
        if args.execute and args.cleanup_fastq:
            raw_1.unlink()
            raw_2.unlink()


if __name__ == "__main__":
    main()
