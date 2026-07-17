#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import urlopen

import pandas as pd


API_ROOT = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon"
TAXA = (
    "Burkholderia cepacia complex",
    "Ralstonia pickettii",
    "Pseudomonas aeruginosa",
    "Stenotrophomonas maltophilia",
)
LEVEL_PRIORITY = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def dataset_report_url(taxon: str, page_size: int, page_token: str = "") -> str:
    query = {"page_size": page_size}
    if page_token:
        query["page_token"] = page_token
    return f"{API_ROOT}/{quote(taxon)}/dataset_report?{urlencode(query)}"


def fetch_candidates(taxon: str, minimum: int, page_size: int = 100, max_pages: int = 3) -> list[dict]:
    candidates: list[dict] = []
    page_token = ""
    for _ in range(max_pages):
        url = dataset_report_url(taxon, page_size, page_token)
        payload = json.loads(urlopen(url, timeout=120).read())
        candidates.extend(payload.get("reports", []))
        page_token = payload.get("next_page_token", "")
        refseq_current = [
            report
            for report in candidates
            if report.get("source_database") == "SOURCE_DATABASE_REFSEQ"
            and report.get("assembly_info", {}).get("assembly_status") == "current"
        ]
        if len(refseq_current) >= minimum or not page_token:
            return refseq_current
    return refseq_current


def candidate_row(taxon: str, report: dict) -> dict:
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    checkm = report.get("checkm_info", {})
    organism = report.get("organism", {})
    accession = report["accession"]
    return {
        "taxon_query": taxon,
        "assembly_accession": accession,
        "organism_name": organism.get("organism_name", ""),
        "strain": organism.get("infraspecific_names", {}).get("strain", ""),
        "biosample_accession": info.get("biosample", {}).get("accession", ""),
        "assembly_level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category", ""),
        "release_date": info.get("release_date", ""),
        "total_sequence_length": int(stats.get("total_sequence_length", 0)),
        "contigs": int(stats.get("number_of_contigs", 0)),
        "checkm_completeness": checkm.get("completeness", ""),
        "checkm_contamination": checkm.get("contamination", ""),
        "ncbi_dataset_page": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/",
        "download_approved": False,
    }


def select_pilot(rows: pd.DataFrame, per_taxon: int) -> pd.DataFrame:
    rows = rows.copy()
    rows["assembly_level_priority"] = rows["assembly_level"].map(LEVEL_PRIORITY).fillna(9)
    rows["checkm_completeness_sort"] = pd.to_numeric(rows["checkm_completeness"], errors="coerce").fillna(-1)
    rows["checkm_contamination_sort"] = pd.to_numeric(rows["checkm_contamination"], errors="coerce").fillna(999)
    rows = rows.sort_values(
        ["taxon_query", "assembly_level_priority", "checkm_contamination_sort", "checkm_completeness_sort", "assembly_accession"],
        ascending=[True, True, True, False, True],
    )
    selected = rows.groupby("taxon_query", sort=False).head(per_taxon).copy()
    return selected.drop(columns=["assembly_level_priority", "checkm_completeness_sort", "checkm_contamination_sort"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-taxon", type=int, default=5)
    parser.add_argument("--manifest", default="data/manifests/pharma_water_ncbi_pilot_plan.csv")
    parser.add_argument("--report", default="reports/pharma_water_ncbi_pilot_plan.json")
    args = parser.parse_args()

    rows = []
    candidate_counts = {}
    for taxon in TAXA:
        candidates = fetch_candidates(taxon, args.per_taxon)
        candidate_counts[taxon] = len(candidates)
        rows.extend(candidate_row(taxon, report) for report in candidates)
        print(f"{taxon}: {len(candidates)} current RefSeq candidates")
    selected = select_pilot(pd.DataFrame(rows), args.per_taxon)
    selected.insert(0, "pilot_id", [f"PW{i:03d}" for i in range(1, len(selected) + 1)])
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(manifest, index=False)

    estimated_bases = int(selected["total_sequence_length"].sum())
    report = {
        "claim_boundary": "Metadata-only pilot plan for exploratory genomic risk prioritization, not tolerance prediction.",
        "source": {"api": "NCBI Datasets API v2alpha", "api_root": API_ROOT},
        "selection": {
            "taxa": list(TAXA),
            "requested_per_taxon": args.per_taxon,
            "candidate_counts": candidate_counts,
            "selected_count": len(selected),
            "selected_by_taxon": selected["taxon_query"].value_counts().to_dict(),
            "estimated_total_sequence_bases": estimated_bases,
            "estimated_uncompressed_fasta_mb_upper_bound": round(estimated_bases / 1_000_000 * 1.1, 2),
        },
        "download_gate": {
            "status": "planned_not_downloaded",
            "reason": "Review the selected RefSeq accessions and projected footprint before downloading assemblies.",
            "next_step": "Download only approved pilot assemblies, fragment them, and score them with both interpretable database matches and calibrated exploratory models.",
        },
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
