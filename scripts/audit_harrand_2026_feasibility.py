#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from urllib.request import urlopen
from zipfile import ZipFile

import pandas as pd


PMC_ID = "PMC12915297"
SUPPLEMENTARY_ZIP_URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMC_ID}/supplementaryFiles"
CORNELL_ITEM_ID = "1f09ec05-8d58-4c70-a477-c55c2e23bdac"
CORNELL_ORIGINAL_BUNDLE_ID = "d659677f-cbe7-4fda-ad0a-088b66175584"
CORNELL_BITSTREAMS_URL = (
    f"https://ecommons.cornell.edu/server/api/core/bundles/{CORNELL_ORIGINAL_BUNDLE_ID}/bitstreams"
)


def download_if_missing(url: str, output: Path) -> None:
    if output.exists():
        print(f"Using cached file: {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading: {url}")
    output.write_bytes(urlopen(url, timeout=120).read())


def load_json(url: str) -> dict:
    return json.loads(urlopen(url, timeout=60).read().decode())


def content_link(bitstream: dict) -> str:
    return bitstream["_links"]["content"]["href"]


def normalize(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="data/raw/harrand_2026_audit")
    parser.add_argument("--current-manifest", default="data/manifests/isolate_manifest_public_assemblies.csv")
    parser.add_argument("--report", default="reports/harrand_2026_feasibility_audit.json")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    supplementary_zip = output_dir / f"{PMC_ID}_supplementaryFiles.zip"
    download_if_missing(SUPPLEMENTARY_ZIP_URL, supplementary_zip)
    with ZipFile(supplementary_zip) as archive:
        workbook_name = "aem.01060-25-s0001.xlsx"
        workbook_path = output_dir / workbook_name
        if not workbook_path.exists():
            workbook_path.write_bytes(archive.read(workbook_name))

    bitstreams = load_json(CORNELL_BITSTREAMS_URL)["_embedded"]["bitstreams"]
    by_name = {bitstream["name"]: bitstream for bitstream in bitstreams}
    cornell_metadata_name = "Harrand_FoodMicrobe_Metadata_20241216.csv"
    cornell_readme_name = "Harrand_FoodMicrobe_readme_20241216.txt"
    cornell_metadata_path = output_dir / cornell_metadata_name
    cornell_readme_path = output_dir / cornell_readme_name
    download_if_missing(content_link(by_name[cornell_metadata_name]), cornell_metadata_path)
    download_if_missing(content_link(by_name[cornell_readme_name]), cornell_readme_path)

    metadata = pd.read_excel(workbook_path, sheet_name="Supplemental Table 1", header=1)
    metadata["Lineage"] = metadata["Lineage"].fillna("").astype(str).str.strip()
    cornell_metadata = pd.read_csv(cornell_metadata_path).dropna(axis=1, how="all")
    cornell_metadata = cornell_metadata.dropna(subset=["Isolate_ID", "File_Name"]).copy()
    current_manifest = pd.read_csv(args.current_manifest).fillna("")
    label_tokens = ("bc", "paa", "naocl", "log_reduction", "log reduction", "survival", "susceptibility")
    label_columns = [
        column for column in metadata.columns
        if any(token in str(column).lower() for token in label_tokens)
    ]

    accessions = metadata["Isolate_Accession"].fillna("").astype(str).str.strip()
    cornell_isolates = set(cornell_metadata["Isolate_ID"].astype(str))
    metadata_isolates = set(metadata["Isolate_ID"].astype(str))
    current_ids = {normalize(value) for value in current_manifest["isolate_id"]}
    overlap = sorted(
        value for value in metadata_isolates
        if normalize(value) in current_ids
    )
    assembly_zip = by_name["Harrand_FoodMicrobe_assemblies_20241216.zip"]

    report = {
        "article": {
            "title": "Listeria sanitizer tolerance at use-level concentrations shows limited association with genetic loci",
            "doi": "10.1128/aem.01060-25",
            "pmc_id": PMC_ID,
            "published_online": "2026-01-20",
        },
        "public_sources": {
            "supplementary_zip_url": SUPPLEMENTARY_ZIP_URL,
            "cornell_item_url": f"https://ecommons.cornell.edu/items/{CORNELL_ITEM_ID}",
            "cornell_assembly_zip_bytes": int(assembly_zip["sizeBytes"]),
            "cornell_assembly_zip_mb": round(int(assembly_zip["sizeBytes"]) / 1_000_000, 2),
        },
        "supplementary_table_s1": {
            "isolates": len(metadata),
            "columns": list(metadata.columns),
            "species_counts": metadata["Species"].value_counts().to_dict(),
            "lineage_counts": metadata["Lineage"].replace("", "<missing>").value_counts().to_dict(),
            "environment_counts": metadata["Environment"].value_counts().to_dict(),
            "l_monocytogenes_isolates": int((metadata["Species"] == "L. monocytogenes").sum()),
            "ecommons_accession_rows": int(accessions.str.contains("ecommons.cornell.edu", regex=False).sum()),
            "non_ecommons_accession_rows": int((~accessions.str.contains("ecommons.cornell.edu", regex=False)).sum()),
            "phenotype_label_columns": label_columns,
        },
        "cornell_deposit": {
            "metadata_isolates": len(cornell_metadata),
            "metadata_columns": list(cornell_metadata.columns),
            "s1_isolates_matching_cornell_metadata": len(metadata_isolates & cornell_isolates),
        },
        "current_project_overlap": {
            "exact_normalized_isolate_overlap": len(overlap),
            "overlapping_isolates": overlap,
        },
        "modeling_gate": {
            "status": "blocked_missing_public_isolate_level_labels" if not label_columns else "go",
            "reason": (
                "Public Table S1 provides isolate metadata and assembly-accession routes, but no isolate-level "
                "BC, PAA, or NaOCl log-reduction labels. The Cornell deposit adds assemblies and metadata only."
            ),
            "required_next_step": (
                "Ask the corresponding author for an isolate-level phenotype table with replicate-level or "
                "averaged log reductions for BC, PAA, and NaOCl before downloading assemblies for model training."
            ),
        },
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    metadata.to_csv(output_dir / "harrand_2026_table_s1_metadata.csv", index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
