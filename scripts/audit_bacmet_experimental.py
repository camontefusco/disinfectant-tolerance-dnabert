#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from urllib.request import urlopen

import pandas as pd


FASTA_URL = "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta"
MAPPING_URL = "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt"
QAC_PATTERN = re.compile(r"quaternary ammonium|benzylkonium|cetrimide|cetylpyridinium|dequalinium", re.I)
BIOCIDE_PATTERN = re.compile(
    r"quaternary ammonium|benzylkonium|cetrimide|cetylpyridinium|dequalinium|"
    r"triclosan|chlorhexidine|phenol|aldehyde|peroxide|hypochlorite|acid|acriflavine",
    re.I,
)
METAL_PATTERN = re.compile(
    r"\b(?:arsenic|cadmium|chromium|cobalt|copper|gold|iron|lead|mercury|nickel|"
    r"silver|tellurium|zinc|antimony|manganese|selenium|vanadium)\b",
    re.I,
)


def download_if_missing(url: str, output: Path) -> None:
    if output.exists():
        print(f"Using cached file: {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading: {url}")
    output.write_bytes(urlopen(url, timeout=120).read())


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header:
        records.append((header, "".join(sequence)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="data/raw/bacmet2_exp")
    parser.add_argument("--manifest", default="data/manifests/bacmet2_exp_annotation_manifest.csv")
    parser.add_argument("--report", default="reports/bacmet2_experimental_annotation_audit.json")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    fasta_path = output_dir / "BacMet2_EXP_database.fasta"
    mapping_path = output_dir / "BacMet2_EXP.753.mapping.txt"
    download_if_missing(FASTA_URL, fasta_path)
    download_if_missing(MAPPING_URL, mapping_path)

    mapping = pd.read_csv(mapping_path, sep="\t").fillna("")
    proteins = pd.DataFrame(
        [
            {"BacMet_ID": header.split("|", 1)[0], "protein_header": header, "protein_sequence": sequence}
            for header, sequence in read_fasta(fasta_path)
        ]
    )
    proteins["protein_variant"] = proteins.groupby("BacMet_ID").cumcount() + 1
    mapping["has_qac_annotation"] = mapping["Compound"].str.contains(QAC_PATTERN)
    mapping["has_biocide_annotation"] = mapping["Compound"].str.contains(BIOCIDE_PATTERN)
    mapping["has_metal_annotation"] = mapping["Compound"].str.contains(METAL_PATTERN)
    data = mapping.merge(proteins, on="BacMet_ID", how="left", validate="many_to_many")
    data["protein_sequence_available"] = data["protein_sequence"].ne("")
    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(manifest_path, index=False)

    alphabet = set("".join(proteins["protein_sequence"].tolist()).upper())
    nucleotide_only = alphabet <= set("ACGTN")
    report = {
        "source": {
            "database": "BacMet 2.0 experimentally confirmed resistance genes",
            "fasta_url": FASTA_URL,
            "mapping_url": MAPPING_URL,
            "fasta_bytes": fasta_path.stat().st_size,
            "mapping_bytes": mapping_path.stat().st_size,
        },
        "inventory": {
            "protein_records": len(proteins),
            "protein_bacmet_ids": int(proteins["BacMet_ID"].nunique()),
            "protein_ids_with_multiple_variants": int((proteins["BacMet_ID"].value_counts() > 1).sum()),
            "mapping_rows": len(mapping),
            "bacmet_ids": int(mapping["BacMet_ID"].nunique()),
            "uniprot_accessions": int(mapping["Accession"].nunique()),
            "qac_annotation_rows": int(mapping["has_qac_annotation"].sum()),
            "biocide_annotation_rows": int(mapping["has_biocide_annotation"].sum()),
            "metal_annotation_rows": int(mapping["has_metal_annotation"].sum()),
            "all_mapping_rows_have_protein": bool(data["protein_sequence_available"].all()),
            "fasta_is_nucleotide_only": nucleotide_only,
        },
        "external_validation_gate": {
            "status": "annotation_only_not_independent_dnabert_nucleotide_validation",
            "reason": (
                "The BacMet experimentally confirmed release is protein FASTA with UniProt accessions. "
                "A DNABERT nucleotide challenge set requires a separately audited protein-to-CDS mapping. "
                "BacMet-derived references have also contributed to MEGARes, so BacMet should be used as "
                "an interpretable annotation layer rather than claimed as independent model validation."
            ),
            "next_use": (
                "Use BacMet gene, compound, and QAC annotations beside exploratory pharmaceutical-water "
                "genome scores. Add a reviewed UniProt-to-CDS mapping only if nucleotide-level annotation "
                "or enrichment analysis is needed."
            ),
        },
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
