#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


FIGSHARE_COLLECTION_ID = "6293046"
BIOPROJECT = "PRJNA884650"
NCBI_TOOL = "disinfectant-tolerance-dnabert"
NCBI_EMAIL = "cmontefusco@example.com"


def fetch_json(url: str, timeout: int = 60) -> dict[str, Any] | list[dict[str, Any]]:
    request = Request(url, headers={"User-Agent": "disinfectant-tolerance-dnabert/0.1"})
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def ncbi_url(endpoint: str, params: dict[str, str | int]) -> str:
    base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi"
    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
    return f"{base}?{urlencode({**params, **common})}"


def is_machine_readable_table_file(file_record: dict[str, Any]) -> bool:
    name = str(file_record.get("name", "")).lower()
    mimetype = str(file_record.get("mimetype", "")).lower()
    extensions = (".csv", ".tsv", ".txt", ".xlsx", ".xls", ".json", ".xml")
    return name.endswith(extensions) or any(token in mimetype for token in ["csv", "excel", "spreadsheet", "json", "xml", "text/plain"])


def fetch_figshare_inventory(collection_id: str = FIGSHARE_COLLECTION_ID) -> pd.DataFrame:
    articles = fetch_json(f"https://api.figshare.com/v2/collections/{collection_id}/articles?page_size=1000")
    rows: list[dict[str, Any]] = []
    for article in articles:
        detail = fetch_json(str(article["url"]))
        files = detail.get("files", []) or [{}]
        for file_record in files:
            rows.append(
                {
                    "article_id": detail.get("id"),
                    "title": detail.get("title"),
                    "doi": detail.get("doi"),
                    "figshare_url": detail.get("figshare_url"),
                    "file_id": file_record.get("id"),
                    "file_name": file_record.get("name", ""),
                    "file_size_bytes": file_record.get("size"),
                    "file_mimetype": file_record.get("mimetype", ""),
                    "download_url": file_record.get("download_url", ""),
                    "machine_readable_table_file": is_machine_readable_table_file(file_record),
                    "mentions_ddac": "ddac" in (str(detail.get("title", "")) + " " + str(detail.get("description", ""))).lower(),
                }
            )
    return pd.DataFrame(rows)


def fetch_assembly_ids(bioproject: str = BIOPROJECT) -> list[str]:
    url = ncbi_url(
        "esearch",
        {"db": "assembly", "term": bioproject, "retmax": 500, "retmode": "json"},
    )
    result = fetch_json(url)
    return list(result["esearchresult"]["idlist"])


def chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[index:index + size] for index in range(0, len(values), size)]


def extract_strain(summary: dict[str, Any]) -> str:
    biosource = summary.get("biosource", {}) or {}
    for item in biosource.get("infraspecieslist", []) or []:
        if str(item.get("sub_type", "")).lower() in {"strain", "isolate"}:
            return str(item.get("sub_value", ""))
    return str(biosource.get("isolate", "") or "")


def parse_assembly_summaries(summary_payloads: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for payload in summary_payloads:
        result = payload.get("result", {})
        for uid in result.get("uids", []):
            item = result[str(uid)]
            rows.append(
                {
                    "uid": uid,
                    "assembly_accession": item.get("assemblyaccession", ""),
                    "genbank_accession": (item.get("synonym", {}) or {}).get("genbank", ""),
                    "refseq_accession": (item.get("synonym", {}) or {}).get("refseq", ""),
                    "assembly_name": item.get("assemblyname", ""),
                    "organism": item.get("organism", ""),
                    "species_name": item.get("speciesname", ""),
                    "biosample_accession": item.get("biosampleaccn", ""),
                    "strain": extract_strain(item),
                    "coverage": item.get("coverage", ""),
                    "assembly_status": item.get("assemblystatus", ""),
                    "contig_n50": item.get("contign50", ""),
                    "scaffold_n50": item.get("scaffoldn50", ""),
                    "ftp_path_refseq": item.get("ftppath_refseq", ""),
                    "ftp_path_genbank": item.get("ftppath_genbank", ""),
                }
            )
    return pd.DataFrame(rows).sort_values(["strain", "assembly_accession"]).reset_index(drop=True)


def fetch_assembly_inventory(bioproject: str = BIOPROJECT) -> pd.DataFrame:
    ids = fetch_assembly_ids(bioproject)
    payloads: list[dict[str, Any]] = []
    for batch in chunked(ids, 25):
        url = ncbi_url("esummary", {"db": "assembly", "id": ",".join(batch), "retmode": "json"})
        payloads.append(fetch_json(url))
        time.sleep(0.34)
    return parse_assembly_summaries(payloads)


def build_decision(figshare: pd.DataFrame, assemblies: pd.DataFrame) -> dict[str, Any]:
    ddac_files = figshare[figshare["mentions_ddac"]].copy()
    raw_ddac_files = ddac_files[ddac_files["machine_readable_table_file"]].copy()
    sequenced_panel_files = figshare[figshare["title"].astype(str).str.contains("sequenced strains|genomic characterization", case=False, regex=True)]
    phenotype_files = figshare[
        figshare["title"].astype(str).str.contains("antimicrobial resistance|decreased susceptibility|DDAC", case=False, regex=True)
        | figshare["mentions_ddac"]
    ]
    can_build_manifest = bool(len(assemblies) >= 1 and len(raw_ddac_files) >= 1)
    return {
        "study": "Pottier et al. 2023 Scientific Reports",
        "doi": "10.1038/s41598-023-29590-0",
        "figshare_collection": "10.6084/m9.figshare.c.6293046.v6",
        "bioproject": BIOPROJECT,
        "figshare_items": int(figshare["article_id"].nunique()),
        "figshare_files": int(len(figshare)),
        "machine_readable_table_files": int(figshare["machine_readable_table_file"].sum()),
        "ddac_mentioning_files": int(len(ddac_files)),
        "machine_readable_ddac_files": int(len(raw_ddac_files)),
        "ncbi_assembly_records": int(len(assemblies)),
        "unique_strains_in_assembly": int(assemblies["strain"].replace("", pd.NA).dropna().nunique()) if "strain" in assemblies else 0,
        "sequenced_panel_figshare_items": sequenced_panel_files[["title", "doi", "file_name", "file_mimetype"]].to_dict("records"),
        "phenotype_figshare_items": phenotype_files[["title", "doi", "file_name", "file_mimetype", "machine_readable_table_file"]].to_dict("records"),
        "can_build_supervised_manifest_without_manual_extraction": can_build_manifest,
        "decision": (
            "blocked_missing_machine_readable_row_level_ddac_labels"
            if not can_build_manifest
            else "candidate_manifest_ready_for_row_level_join"
        ),
        "recommended_next_action": (
            "Do not train yet. Assemblies are available through NCBI Assembly, but the public figshare phenotype materials are image/SVG/TIFF files rather than a row-level CSV/XLSX table. "
            "Contact authors or perform documented manual/OCR extraction before building a supervised DDAC manifest."
            if not can_build_manifest
            else "Join raw DDAC phenotype rows to NCBI strain IDs and audit label balance before modeling."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figshare-output", default="data/manifests/pottier_2023_figshare_inventory.csv")
    parser.add_argument("--assembly-output", default="data/manifests/pottier_2023_ncbi_assembly_inventory.csv")
    parser.add_argument("--decision-output", default="reports/pottier_2023_ddac_manifest_verification.json")
    args = parser.parse_args()

    figshare = fetch_figshare_inventory()
    assemblies = fetch_assembly_inventory()
    decision = build_decision(figshare, assemblies)

    Path(args.figshare_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.assembly_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.decision_output).parent.mkdir(parents=True, exist_ok=True)
    figshare.to_csv(args.figshare_output, index=False)
    assemblies.to_csv(args.assembly_output, index=False)
    Path(args.decision_output).write_text(json.dumps(decision, indent=2) + "\n")
    print(f"Wrote {args.figshare_output}")
    print(f"Wrote {args.assembly_output}")
    print(f"Wrote {args.decision_output}")
    print(decision["decision"])


if __name__ == "__main__":
    main()
