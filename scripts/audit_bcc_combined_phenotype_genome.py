#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


MOORE_2009_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC2640157/"
KIM_2015_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12381686/"
RUSHTON_2013_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3697374/"
NCBI_TOOL = "disinfectant-tolerance-dnabert"
NCBI_EMAIL = "cmontefusco@example.com"


def normalize_text(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def flatten_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    if isinstance(frame.columns, pd.MultiIndex):
        names = []
        for column in frame.columns:
            parts = [str(part) for part in column if not str(part).startswith("Unnamed")]
            names.append(" | ".join(parts))
        frame.columns = names
    else:
        frame.columns = [str(column) for column in frame.columns]
    return frame


def clean_value(value: object) -> str:
    text = str(value).strip()
    if text.lower() in {"nan", "none", ""}:
        return ""
    return text.replace("\u200a", "").replace("\xa0", " ")


def extract_moore_2009(tables: list[pd.DataFrame]) -> pd.DataFrame:
    table = flatten_columns(tables[1])
    table.columns = [
        "species_or_strain",
        "chx_mic_mg_l",
        "chx_mbc_mg_l",
        "cpc_mic_mg_l",
        "cpc_mbc_mg_l",
        "triclosan_mic_mg_l",
        "bzk_mic_mg_l",
    ]
    rows: list[dict[str, Any]] = []
    current_species = ""
    for _, row in table.iterrows():
        name = clean_value(row["species_or_strain"])
        if not name:
            continue
        if name.startswith("B.") and "(" in name and ")" in name:
            current_species = re.sub(r"\s*\([^)]*\)", "", name).strip()
            continue
        for compound, endpoint, value_column in [
            ("chlorhexidine", "MIC", "chx_mic_mg_l"),
            ("chlorhexidine", "MBC", "chx_mbc_mg_l"),
            ("cetylpyridinium chloride", "MIC", "cpc_mic_mg_l"),
            ("cetylpyridinium chloride", "MBC", "cpc_mbc_mg_l"),
            ("triclosan", "MIC", "triclosan_mic_mg_l"),
            ("benzalkonium chloride", "MIC", "bzk_mic_mg_l"),
        ]:
            value = clean_value(row[value_column])
            if value:
                rows.append(
                    {
                        "study_key": "moore_2009_bcc_biocide_susceptibility",
                        "citation_short": "Moore et al. 2009, Journal of Antimicrobial Chemotherapy",
                        "doi": "10.1093/jac/dkn540",
                        "source_url": MOORE_2009_URL,
                        "species_reported": current_species,
                        "strain_id": name,
                        "isolation_source": "",
                        "compound": compound,
                        "endpoint": endpoint,
                        "timepoint": "baseline",
                        "value": value,
                        "units": "mg/L",
                        "row_level_phenotype": True,
                    }
                )
    return pd.DataFrame(rows)


def extract_kim_2015(tables: list[pd.DataFrame]) -> pd.DataFrame:
    table = flatten_columns(tables[0])
    table.columns = [
        "species_reported",
        "strain_id",
        "isolation_source",
        "chx_initial_ug_ml",
        "chx_day40_ug_ml",
        "bzk_initial_ug_ml",
        "bzk_day40_ug_ml",
    ]
    rows: list[dict[str, Any]] = []
    for _, row in table.iterrows():
        strain = clean_value(row["strain_id"])
        if not strain:
            continue
        for compound, timepoint, value_column in [
            ("chlorhexidine", "baseline", "chx_initial_ug_ml"),
            ("chlorhexidine", "day_40", "chx_day40_ug_ml"),
            ("benzalkonium chloride", "baseline", "bzk_initial_ug_ml"),
            ("benzalkonium chloride", "day_40", "bzk_day40_ug_ml"),
        ]:
            value = clean_value(row[value_column])
            if value:
                rows.append(
                    {
                        "study_key": "kim_2015_bcc_chx_bzk",
                        "citation_short": "Kim/Wang et al. 2015, J Ind Microbiol Biotechnol",
                        "doi": "10.1007/s10295-015-1605-x",
                        "source_url": KIM_2015_URL,
                        "species_reported": clean_value(row["species_reported"]),
                        "strain_id": strain,
                        "isolation_source": clean_value(row["isolation_source"]),
                        "compound": compound,
                        "endpoint": "MIC",
                        "timepoint": timepoint,
                        "value": value,
                        "units": "ug/mL",
                        "row_level_phenotype": True,
                    }
                )
    return pd.DataFrame(rows)


def extract_public_phenotypes() -> pd.DataFrame:
    moore_tables = pd.read_html(MOORE_2009_URL)
    kim_tables = pd.read_html(KIM_2015_URL)
    rows = pd.concat([extract_moore_2009(moore_tables), extract_kim_2015(kim_tables)], ignore_index=True)
    rows["strain_norm"] = rows["strain_id"].map(normalize_text)
    return rows


def fetch_json(url: str, timeout: int = 60) -> dict[str, Any]:
    request = Request(url, headers={"User-Agent": "disinfectant-tolerance-dnabert/0.1"})
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def ncbi_url(endpoint: str, params: dict[str, str | int]) -> str:
    base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi"
    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
    return f"{base}?{urlencode({**params, **common})}"


def extract_strain(summary: dict[str, Any]) -> str:
    biosource = summary.get("biosource", {}) or {}
    for item in biosource.get("infraspecieslist", []) or []:
        if str(item.get("sub_type", "")).lower() in {"strain", "isolate"}:
            return str(item.get("sub_value", ""))
    return str(biosource.get("isolate", "") or "")


def summarize_assembly(uid: str, item: dict[str, Any], query_strain: str) -> dict[str, Any]:
    strain = extract_strain(item)
    haystack = " ".join(
        [
            strain,
            str(item.get("assemblyname", "")),
            str(item.get("organism", "")),
            str(item.get("biosampleaccn", "")),
            str(item.get("assemblyaccession", "")),
        ]
    )
    return {
        "query_strain": query_strain,
        "query_norm": normalize_text(query_strain),
        "uid": uid,
        "assembly_accession": item.get("assemblyaccession", ""),
        "genbank_accession": (item.get("synonym", {}) or {}).get("genbank", ""),
        "refseq_accession": (item.get("synonym", {}) or {}).get("refseq", ""),
        "organism": item.get("organism", ""),
        "biosample_accession": item.get("biosampleaccn", ""),
        "assembly_strain": strain,
        "assembly_status": item.get("assemblystatus", ""),
        "coverage": item.get("coverage", ""),
        "ftp_path_refseq": item.get("ftppath_refseq", ""),
        "ftp_path_genbank": item.get("ftppath_genbank", ""),
        "exact_norm_match": normalize_text(query_strain) == normalize_text(strain),
        "contains_norm_match": normalize_text(query_strain) in normalize_text(haystack),
    }


def query_assembly_for_strain(strain: str, retmax: int = 5) -> list[dict[str, Any]]:
    term = f'Burkholderia[Organism] AND "{strain}"[All Fields]'
    search = fetch_json(ncbi_url("esearch", {"db": "assembly", "term": term, "retmax": retmax, "retmode": "json"}))
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    summary = fetch_json(ncbi_url("esummary", {"db": "assembly", "id": ",".join(ids), "retmode": "json"}))
    result = summary.get("result", {})
    return [summarize_assembly(uid, result[str(uid)], strain) for uid in result.get("uids", [])]


def query_all_strains(strains: list[str], sleep_seconds: float = 0.34) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for strain in strains:
        rows.extend(query_assembly_for_strain(strain))
        time.sleep(sleep_seconds)
    if rows:
        return pd.DataFrame(rows).sort_values(["query_strain", "exact_norm_match", "contains_norm_match"], ascending=[True, False, False])
    return pd.DataFrame(
        columns=[
            "query_strain", "query_norm", "uid", "assembly_accession", "genbank_accession",
            "refseq_accession", "organism", "biosample_accession", "assembly_strain",
            "assembly_status", "coverage", "ftp_path_refseq", "ftp_path_genbank",
            "exact_norm_match", "contains_norm_match",
        ]
    )


def build_candidate_manifest(phenotypes: pd.DataFrame, matches: pd.DataFrame) -> pd.DataFrame:
    if matches.empty:
        return pd.DataFrame()
    best = (
        matches[matches["exact_norm_match"] | matches["contains_norm_match"]]
        .sort_values(["query_strain", "exact_norm_match", "contains_norm_match"], ascending=[True, False, False])
        .drop_duplicates("query_strain")
    )
    merged = phenotypes.merge(best, left_on="strain_id", right_on="query_strain", how="inner", validate="many_to_one")
    return merged


def build_report(phenotypes: pd.DataFrame, matches: pd.DataFrame, manifest: pd.DataFrame) -> dict[str, Any]:
    bzk = phenotypes[phenotypes["compound"].eq("benzalkonium chloride")]
    chx = phenotypes[phenotypes["compound"].eq("chlorhexidine")]
    matched_strains = set(manifest["strain_id"]) if not manifest.empty else set()
    return {
        "studies_used_for_row_level_phenotypes": sorted(phenotypes["study_key"].unique().tolist()),
        "context_studies_not_row_level_from_public_html": [
            {
                "study_key": "rushton_2013_bcc_preservatives",
                "doi": "10.1128/AAC.00140-13",
                "source_url": RUSHTON_2013_URL,
                "reason": "PMC HTML tables expose aggregate species/preservative counts and B. lata adaptive variants, not a row-level 83-strain MIC/MBC table.",
            }
        ],
        "phenotype_rows": int(len(phenotypes)),
        "unique_phenotyped_strains": int(phenotypes["strain_id"].nunique()),
        "unique_bzk_strains": int(bzk["strain_id"].nunique()),
        "unique_chx_strains": int(chx["strain_id"].nunique()),
        "assembly_match_candidate_rows": int(len(matches)),
        "exact_or_contains_matched_strains": int(len(matched_strains)),
        "candidate_manifest_rows": int(len(manifest)),
        "candidate_manifest_unique_strains": int(manifest["strain_id"].nunique()) if not manifest.empty else 0,
        "candidate_manifest_unique_bzk_strains": int(manifest[manifest["compound"].eq("benzalkonium chloride")]["strain_id"].nunique()) if not manifest.empty else 0,
        "decision": "candidate_supervised_bcc_manifest_created" if len(matched_strains) >= 20 else "insufficient_exact_public_genome_matches_for_supervised_model",
        "recommended_next_action": (
            "Audit candidate genome matches manually, choose a harmonized BZK/CHX endpoint, and create grouped splits."
            if len(matched_strains) >= 20
            else "Use these Bcc phenotypes for biological context unless more strain-to-genome matches are found."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotypes-output", default="data/manifests/bcc_combined_biocide_phenotype_rows.csv")
    parser.add_argument("--matches-output", default="data/manifests/bcc_ncbi_assembly_match_candidates.csv")
    parser.add_argument("--manifest-output", default="data/manifests/bcc_candidate_supervised_manifest.csv")
    parser.add_argument("--report-output", default="reports/bcc_combined_phenotype_genome_audit.json")
    parser.add_argument("--skip-ncbi", action="store_true")
    args = parser.parse_args()

    phenotypes = extract_public_phenotypes()
    strains = sorted(phenotypes["strain_id"].dropna().unique().tolist())
    matches = pd.DataFrame() if args.skip_ncbi else query_all_strains(strains)
    manifest = build_candidate_manifest(phenotypes, matches)
    report = build_report(phenotypes, matches, manifest)

    Path(args.phenotypes_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.matches_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.manifest_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report_output).parent.mkdir(parents=True, exist_ok=True)
    phenotypes.to_csv(args.phenotypes_output, index=False)
    matches.to_csv(args.matches_output, index=False)
    manifest.to_csv(args.manifest_output, index=False)
    Path(args.report_output).write_text(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {args.phenotypes_output}")
    print(f"Wrote {args.matches_output}")
    print(f"Wrote {args.manifest_output}")
    print(f"Wrote {args.report_output}")
    print(report["decision"])


if __name__ == "__main__":
    main()
