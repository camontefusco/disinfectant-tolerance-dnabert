#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess

import pandas as pd


def build_queries(data: pd.DataFrame, path: Path) -> pd.DataFrame:
    selected = data[data["has_biocide_annotation"]].drop_duplicates(["BacMet_ID", "protein_variant"]).copy()
    selected["query_id"] = selected["BacMet_ID"] + "_v" + selected["protein_variant"].astype(int).astype(str)
    with path.open("w") as handle:
        for row in selected.itertuples():
            handle.write(f">{row.query_id}\n{row.protein_sequence}\n")
    return selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", default="data/manifests/pharma_water_downloaded_assemblies.csv")
    parser.add_argument("--bacmet", default="data/manifests/bacmet2_exp_annotation_manifest.csv")
    parser.add_argument("--work-dir", default="data/processed/pharma_water_bacmet")
    parser.add_argument("--output", default="data/manifests/pharma_water_pilot_bacmet_biocide_hits.csv")
    args = parser.parse_args()

    if shutil.which("makeblastdb") is None or shutil.which("tblastn") is None:
        raise RuntimeError("NCBI BLAST+ with tblastn is required")
    assemblies = pd.read_csv(args.assemblies)
    bacmet = pd.read_csv(args.bacmet).fillna("")
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    query = work_dir / "bacmet_biocide_proteins.faa"
    selected = build_queries(bacmet, query)
    metadata = selected.set_index("query_id")[["BacMet_ID", "Gene_name", "Compound", "has_qac_annotation"]]
    outputs = []
    for row in assemblies.itertuples():
        database = work_dir / row.assembly_accession
        hits_path = work_dir / f"{row.assembly_accession}.tsv"
        if not hits_path.exists():
            subprocess.run(["makeblastdb", "-in", row.assembly_path, "-dbtype", "nucl", "-out", str(database)], check=True, stdout=subprocess.DEVNULL)
            subprocess.run(
                [
                    "tblastn", "-query", str(query), "-db", str(database),
                    "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
                    "-qcov_hsp_perc", "80", "-max_hsps", "1", "-evalue", "1e-20", "-out", str(hits_path),
                ],
                check=True,
            )
        names = ["query_id", "contig", "identity", "alignment_length", "query_length", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        hits = pd.read_csv(hits_path, sep="\t", names=names) if hits_path.stat().st_size else pd.DataFrame(columns=names)
        if not hits.empty:
            hits = hits[hits["identity"] >= 80].copy()
            hits["assembly_accession"] = row.assembly_accession
            hits["taxon_query"] = row.taxon_query
            outputs.append(hits)
        print(f"{row.assembly_accession}: {len(hits)} BacMet biocide-protein hits")
    hits = pd.concat(outputs, ignore_index=True) if outputs else pd.DataFrame()
    if not hits.empty:
        hits = hits.merge(metadata, on="query_id", how="left", validate="many_to_one")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    hits.to_csv(output_path, index=False)
    print(f"Wrote {len(hits)} BacMet biocide-protein hits to {output_path}")


if __name__ == "__main__":
    main()
