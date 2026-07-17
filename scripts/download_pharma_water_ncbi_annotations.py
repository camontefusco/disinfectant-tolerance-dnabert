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
    return f"{API_ROOT}/{accession}/download?include_annotation_type=GENOME_FASTA,PROT_FASTA,GENOME_GFF"


def extract_annotation(accession: str, output_dir: Path) -> tuple[Path, Path]:
    protein = output_dir / f"{accession}.faa"
    gff = output_dir / f"{accession}.gff"
    if protein.exists() and gff.exists():
        print(f"Using cached annotations: {accession}")
        return protein, gff
    print(f"Downloading annotated NCBI package: {accession}")
    archive = ZipFile(BytesIO(urlopen(package_url(accession), timeout=180).read()))
    protein_names = [name for name in archive.namelist() if name.endswith("/protein.faa")]
    gff_names = [name for name in archive.namelist() if name.endswith("/genomic.gff")]
    if len(protein_names) != 1 or len(gff_names) != 1:
        raise ValueError(f"Expected protein FASTA and GFF for {accession}; found {protein_names}, {gff_names}")
    output_dir.mkdir(parents=True, exist_ok=True)
    protein.write_bytes(archive.read(protein_names[0]))
    gff.write_bytes(archive.read(gff_names[0]))
    return protein, gff


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", default="data/manifests/pharma_water_downloaded_assemblies.csv")
    parser.add_argument("--output-dir", default="data/raw/pharma_water_pilot_annotations")
    parser.add_argument("--manifest", default="data/manifests/pharma_water_pilot_annotation_files.csv")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    if not args.execute:
        raise RuntimeError("Annotation download remains gated. Re-run with --execute.")
    assemblies = pd.read_csv(args.assemblies)
    output_dir = Path(args.output_dir)
    rows = []
    for row in assemblies.itertuples():
        protein, gff = extract_annotation(row.assembly_accession, output_dir)
        rows.append(
            {
                "assembly_accession": row.assembly_accession, "taxon_query": row.taxon_query,
                "assembly_path": row.assembly_path, "protein_path": str(protein), "gff_path": str(gff),
                "protein_bytes": protein.stat().st_size, "gff_bytes": gff.stat().st_size,
            }
        )
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    data = pd.DataFrame(rows)
    data.to_csv(manifest, index=False)
    print(f"Wrote {len(data)} annotation packages to {manifest}")
    print(f"Protein and GFF footprint MB: {(data['protein_bytes'].sum() + data['gff_bytes'].sum()) / 1_000_000:.2f}")


if __name__ == "__main__":
    main()
