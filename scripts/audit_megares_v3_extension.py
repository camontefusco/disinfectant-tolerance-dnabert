#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from urllib.request import urlopen

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.io import read_fasta


FASTA_URL = "https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta"
ANNOTATIONS_URL = "https://www.meglab.org/downloads/megares_v3.00/megares_annotations_v3.00.csv"
TARGET_TYPES = {"Biocides", "Metals", "Multi-compound"}


def download_if_missing(url: str, output: Path) -> None:
    if output.exists():
        print(f"Using cached file: {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading: {url}")
    output.write_bytes(urlopen(url, timeout=120).read())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="data/raw/megares_v3")
    parser.add_argument("--manifest", default="data/manifests/megares_biocide_metal_manifest.csv")
    parser.add_argument("--report", default="reports/megares_v3_extension_audit.json")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    fasta_path = output_dir / "megares_database_v3.00.fasta"
    annotations_path = output_dir / "megares_annotations_v3.00.csv"
    download_if_missing(FASTA_URL, fasta_path)
    download_if_missing(ANNOTATIONS_URL, annotations_path)

    annotations = pd.read_csv(annotations_path)
    sequences = pd.DataFrame(
        [
            {"header": header, "sequence": sequence, "sequence_length": len(sequence)}
            for header, sequence in read_fasta(fasta_path)
        ]
    )
    data = annotations.merge(sequences, on="header", how="left", validate="one_to_one")
    if data["sequence"].isna().any():
        raise ValueError("Every MEGARes annotation must have a matching nucleotide sequence")
    subset = data[data["type"].isin(TARGET_TYPES)].copy()
    subset["requires_snp_confirmation"] = subset["header"].str.contains("RequiresSNPConfirmation", regex=False)
    subset["target_family"] = subset["type"].where(subset["type"] != "Multi-compound", "Multi-compound")
    subset["homology_group"] = subset["group"]
    subset["fragment_id"] = subset["header"].str.split("|", regex=False).str[0]
    manifest_columns = [
        "fragment_id", "header", "type", "class", "mechanism", "group",
        "homology_group", "requires_snp_confirmation", "sequence_length", "sequence",
    ]
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    subset[manifest_columns].to_csv(manifest, index=False)

    report = {
        "source": {
            "database": "MEGARes v3.0",
            "fasta_url": FASTA_URL,
            "annotations_url": ANNOTATIONS_URL,
            "fasta_bytes": fasta_path.stat().st_size,
            "annotations_bytes": annotations_path.stat().st_size,
        },
        "all_sequences": {
            "count": len(data),
            "type_counts": data["type"].value_counts().to_dict(),
        },
        "exploratory_subset": {
            "types": sorted(TARGET_TYPES),
            "count": len(subset),
            "type_counts": subset["type"].value_counts().to_dict(),
            "class_counts": subset["class"].value_counts().to_dict(),
            "mechanism_count": int(subset["mechanism"].nunique()),
            "group_count": int(subset["group"].nunique()),
            "requires_snp_confirmation": int(subset["requires_snp_confirmation"].sum()),
            "sequence_length": {
                "min": int(subset["sequence_length"].min()),
                "median": float(subset["sequence_length"].median()),
                "max": int(subset["sequence_length"].max()),
            },
        },
        "modeling_gate": {
            "status": "go_exploratory_fragment_classifier",
            "claim_boundary": (
                "Use MEGARes labels to classify resistance-associated fragments and prioritize genomic regions. "
                "Do not interpret the output as validated sanitizer-survival prediction for an isolate."
            ),
            "recommended_split": (
                "Hold out entire MEGARes group labels during evaluation to reduce close-homolog leakage. "
                "For stronger publication work, add nucleotide-similarity clustering before splitting."
            ),
        },
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
