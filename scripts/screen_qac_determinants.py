#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd
from sklearn.metrics import balanced_accuracy_score, roc_auc_score


OUTFMT = "6 qseqid sseqid pident length qlen evalue bitscore"
COLUMNS = ["reference", "contig", "identity", "alignment_length", "reference_length", "evalue", "bitscore"]
RESULT_COLUMNS = [*COLUMNS, "gene", "coverage"]


def best_hits(reference_fasta: Path, assembly_path: Path) -> pd.DataFrame:
    result = subprocess.run(
        [
            "blastn",
            "-query",
            str(reference_fasta),
            "-subject",
            str(assembly_path),
            "-outfmt",
            OUTFMT,
            "-dust",
            "no",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not result.stdout.strip():
        return pd.DataFrame(columns=RESULT_COLUMNS)
    rows = [line.split("\t") for line in result.stdout.splitlines()]
    frame = pd.DataFrame(rows, columns=COLUMNS)
    for column in ("identity", "alignment_length", "reference_length", "evalue", "bitscore"):
        frame[column] = pd.to_numeric(frame[column])
    frame["gene"] = frame["reference"].str.split("|", regex=False).str[0]
    frame["coverage"] = frame["alignment_length"] / frame["reference_length"]
    return frame.sort_values(["gene", "bitscore"], ascending=[True, False]).drop_duplicates("gene")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--references", default="data/references/qac_determinants.fasta")
    parser.add_argument("--output", default="data/processed/qac_determinant_screen.csv")
    parser.add_argument("--report", default="reports/qac_gene_rule_baseline.json")
    parser.add_argument("--minimum-identity", type=float, default=80.0)
    parser.add_argument("--minimum-coverage", type=float, default=0.8)
    args = parser.parse_args()

    if shutil.which("blastn") is None:
        raise RuntimeError("blastn is required. Install NCBI BLAST+ first.")
    manifest = pd.read_csv(args.manifest, dtype={"label": int}).fillna("")
    hits = []
    for index, row in manifest.iterrows():
        frame = best_hits(Path(args.references), Path(row["assembly_path"]))
        for gene in ("bcrABC", "qacH", "emrE", "emrC"):
            match = frame[frame["gene"] == gene]
            record = {"isolate_id": row["isolate_id"], "gene": gene, "identity": 0.0, "coverage": 0.0, "present": False}
            if not match.empty:
                best = match.iloc[0]
                record.update({"identity": float(best["identity"]), "coverage": float(best["coverage"])})
                record["present"] = record["identity"] >= args.minimum_identity and record["coverage"] >= args.minimum_coverage
            hits.append(record)
        if (index + 1) % 25 == 0 or index + 1 == len(manifest):
            print(f"Screened {index + 1}/{len(manifest)} assemblies")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    hit_frame = pd.DataFrame(hits)
    hit_frame.to_csv(output, index=False)
    presence = hit_frame.pivot(index="isolate_id", columns="gene", values="present").reset_index()
    data = manifest.merge(presence, on="isolate_id", validate="one_to_one")
    gene_columns = ["bcrABC", "qacH", "emrE", "emrC"]
    data["known_qac_gene"] = data[gene_columns].any(axis=1)
    predictions = data["known_qac_gene"].astype(int)
    unexplained = data[(data["label"] == 1) & ~data["known_qac_gene"]]["isolate_id"].tolist()
    report = {
        "model": "known_qac_determinant_blastn_rule",
        "minimum_identity": args.minimum_identity,
        "minimum_coverage": args.minimum_coverage,
        "balanced_accuracy": balanced_accuracy_score(data["label"], predictions),
        "roc_auc": roc_auc_score(data["label"], predictions),
        "isolate_count": len(data),
        "gene_positive_isolates": int(data["known_qac_gene"].sum()),
        "gene_counts": {gene: int(data[gene].sum()) for gene in gene_columns},
        "unexplained_tolerant_count": len(unexplained),
        "unexplained_tolerant_isolates": unexplained,
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
