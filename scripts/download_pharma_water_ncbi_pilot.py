#!/usr/bin/env python3
from __future__ import annotations

import argparse
from io import BytesIO
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

import pandas as pd


API_ROOT = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"


def package_url(accession: str) -> str:
    return f"{API_ROOT}/{accession}/download?include_annotation_type=GENOME_FASTA"


def download_fasta(accession: str, output: Path) -> None:
    if output.exists():
        print(f"Using cached assembly: {output}")
        return
    print(f"Downloading compact NCBI package: {accession}")
    archive = ZipFile(BytesIO(urlopen(package_url(accession), timeout=180).read()))
    fasta_names = [name for name in archive.namelist() if name.endswith("_genomic.fna")]
    if len(fasta_names) != 1:
        raise ValueError(f"Expected one genomic FASTA for {accession}, found {fasta_names}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(archive.read(fasta_names[0]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", default="data/manifests/pharma_water_ncbi_pilot_plan.csv")
    parser.add_argument("--output-dir", default="data/raw/pharma_water_pilot_assemblies")
    parser.add_argument("--manifest", default="data/manifests/pharma_water_downloaded_assemblies.csv")
    parser.add_argument("--max-estimated-mb", type=float, default=200.0)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    plan = pd.read_csv(args.plan)
    estimated_mb = plan["total_sequence_length"].sum() / 1_000_000 * 1.1
    if estimated_mb > args.max_estimated_mb:
        raise ValueError(f"Estimated FASTA footprint {estimated_mb:.2f} MB exceeds cap {args.max_estimated_mb:.2f} MB")
    if not args.execute:
        raise RuntimeError("Download remains gated. Re-run with --execute after reviewing the pilot plan.")
    output_dir = Path(args.output_dir)
    paths = []
    for accession in plan["assembly_accession"]:
        output = output_dir / f"{accession}.fna"
        download_fasta(accession, output)
        paths.append(str(output))
    plan["assembly_path"] = paths
    plan["downloaded"] = [Path(path).exists() for path in paths]
    plan["fasta_bytes"] = [Path(path).stat().st_size for path in paths]
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    plan.to_csv(manifest, index=False)
    print(f"Wrote {len(plan)} assemblies to {manifest}")
    print(f"FASTA footprint MB: {plan['fasta_bytes'].sum() / 1_000_000:.2f}")


if __name__ == "__main__":
    main()
