#!/usr/bin/env python3
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import subprocess

import pandas as pd


def gff_coordinates(path: str) -> pd.DataFrame:
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            match = re.search(r"(?:^|;)protein_id=([^;]+)", fields[8])
            if match:
                rows.append(
                    {
                        "Protein id": match.group(1), "Contig id": fields[0],
                        "Start": int(fields[3]), "Stop": int(fields[4]), "Strand": fields[6],
                    }
                )
    return pd.DataFrame(rows).drop_duplicates("Protein id")


def annotate(row: object, args: argparse.Namespace, work_dir: Path) -> pd.DataFrame:
    output = work_dir / f"{row.assembly_accession}.tsv"
    if not output.exists():
        command = [
            args.amrfinder, "-p", row.protein_path, "--plus", "--name", row.assembly_accession,
            "-d", args.database, "--threads", str(args.threads_per_assembly), "-o", str(output),
        ]
        print("Running:", " ".join(command))
        subprocess.run(command, check=True)
    else:
        print(f"Using cached AMRFinderPlus report: {row.assembly_accession}")
    report = pd.read_csv(output, sep="\t")
    report = report.merge(gff_coordinates(row.gff_path), on="Protein id", how="left", validate="many_to_one")
    report["assembly_accession"] = row.assembly_accession
    report["taxon_query"] = row.taxon_query
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", default="data/manifests/pharma_water_pilot_annotation_files.csv")
    parser.add_argument("--amrfinder", default=".tools/amrfinderplus/bin/amrfinder")
    parser.add_argument("--database", default="data/external/amrfinderplus_db/latest")
    parser.add_argument("--work-dir", default="data/processed/pharma_water_amrfinderplus_protein")
    parser.add_argument("--output", default="data/manifests/pharma_water_pilot_amrfinderplus_hits.csv")
    parser.add_argument("--report", default="reports/pharma_water_pilot_amrfinderplus.json")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads-per-assembly", type=int, default=2)
    args = parser.parse_args()

    database = Path(args.database)
    database_version = database.resolve().name
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    annotations = pd.read_csv(args.annotations)
    rows = list(annotations.itertuples())
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        outputs = list(pool.map(lambda row: annotate(row, args, work_dir), rows))
    hits = pd.concat(outputs, ignore_index=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    hits.to_csv(output_path, index=False)
    summary = {
        "software": "NCBI AMRFinderPlus",
        "database_version": database_version,
        "mode": "PGAP protein inputs with --plus; coordinates restored from NCBI PGAP GFF. Generic cross-organism layer without organism-specific mutation calling.",
        "assemblies": len(annotations),
        "hits": len(hits),
        "scope_counts": hits["Scope"].value_counts().to_dict(),
        "type_counts": hits["Type"].value_counts().to_dict(),
        "class_counts": hits["Class"].value_counts().to_dict(),
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
