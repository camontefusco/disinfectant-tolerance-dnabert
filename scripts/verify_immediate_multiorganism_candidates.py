#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import urljoin, urlencode
from urllib.request import Request, urlopen

import pandas as pd
from bs4 import BeautifulSoup


NCBI_TOOL = "disinfectant-tolerance-dnabert"
NCBI_EMAIL = "cmontefusco@example.com"

IMMEDIATE_CANDIDATES = [
    {
        "study_key": "proteus_2025_chx_cationic",
        "citation_short": "Bennett et al. 2025, Microbiology",
        "organism_group": "Proteus mirabilis",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12304739/",
        "doi": "10.1099/mic.0.001580",
        "pmcid": "PMC12304739",
        "reported_endpoint": "chlorhexidine and other cationic biocides",
        "reported_project": "PRJNA1154625; PRJNA475751; PRJNA12624",
        "phenotype_terms": ["chlorhexidine", "CHD", "cationic biocide", "MIC", "Table S1"],
        "accession_terms": ["PRJNA1154625", "PRJNA475751", "PRJNA12624", "accession", "Table S1"],
        "expected_status": "most_promising_if_table_s1_is_machine_readable",
    },
    {
        "study_key": "enterococcus_faecium_2019_chx",
        "citation_short": "Duarte et al. 2019, Applied and Environmental Microbiology",
        "organism_group": "Enterococcus faecium",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6856339/",
        "doi": "10.1128/AEM.01589-19",
        "pmcid": "PMC6856339",
        "reported_endpoint": "chlorhexidine MIC distribution for 106 isolates",
        "reported_project": "",
        "phenotype_terms": ["chlorhexidine", "MIC", "106", "Supplemental"],
        "accession_terms": ["accession", "BioProject", "genome", "sequence"],
        "expected_status": "strong_if_supplement_has_accessions",
    },
    {
        "study_key": "enterococcus_faecalis_2022_chx",
        "citation_short": "Pereira et al. 2022, mSphere",
        "organism_group": "Enterococcus faecalis",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9430118/",
        "doi": "10.1128/msphere.00486-22",
        "pmcid": "PMC9430118",
        "reported_endpoint": "chlorhexidine MIC/MBC and ChlR-EfrEF evolution",
        "reported_project": "",
        "phenotype_terms": ["chlorhexidine", "MIC", "MBC", "Supplementary Table"],
        "accession_terms": ["accession", "BioProject", "genome", "sequence"],
        "expected_status": "second_tier_if_accessions_are_public",
    },
    {
        "study_key": "hospital_environment_chx_2024",
        "citation_short": "Hartmann lab 2024, medRxiv",
        "organism_group": "Mixed hospital environmental bacteria",
        "url": "https://www.medrxiv.org/content/10.1101/2024.10.07.24315058v1.full-text",
        "doi": "10.1101/2024.10.07.24315058",
        "pmcid": "",
        "reported_endpoint": "chlorhexidine MICs and ARG/taxonomy table",
        "reported_project": "PRJNA1169385",
        "phenotype_terms": ["chlorhexidine MIC", "Table S4", "taxonomic", "ARG"],
        "accession_terms": ["PRJNA1169385", "SRA", "GitHub"],
        "expected_status": "high_value_but_preprint_and_data_release_risk",
    },
]


def fetch_text(url: str, timeout: int = 60) -> str:
    request = Request(url, headers={"User-Agent": "disinfectant-tolerance-dnabert/0.1"})
    with urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="ignore")


def cached_text(raw_dir: Path, study_key: str, url: str) -> str:
    path = raw_dir / f"{study_key}.html"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(fetch_text(url))
    return path.read_text(errors="ignore")


def normalize_space(value: str) -> str:
    return re.sub(r"\s+", " ", value).strip()


def find_supplement_links(html: str, base_url: str) -> list[dict[str, object]]:
    soup = BeautifulSoup(html, "html.parser")
    rows: list[dict[str, object]] = []
    for link in soup.find_all("a", href=True):
        href = str(link.get("href", ""))
        text = normalize_space(link.get_text(" ", strip=True))
        combined = f"{text} {href}"
        if not re.search(r"supp|table|data|xls|xlsx|csv|tsv|docx|pdf|figshare|github", combined, re.I):
            continue
        absolute = urljoin(base_url, href)
        lower = absolute.lower()
        rows.append(
            {
                "link_text": text[:300],
                "href": absolute,
                "extension": Path(lower.split("?")[0]).suffix,
                "machine_readable_candidate": bool(re.search(r"\.(csv|tsv|xlsx|xls|json|txt)(\?|$)", lower)),
                "spreadsheet_candidate": bool(re.search(r"\.(xlsx|xls|csv|tsv)(\?|$)", lower)),
            }
        )
    # Remove exact duplicate links while preserving order.
    seen: set[str] = set()
    unique: list[dict[str, object]] = []
    for row in rows:
        href = str(row["href"])
        if href in seen:
            continue
        seen.add(href)
        unique.append(row)
    return unique


def table_inventory(html: str) -> list[dict[str, object]]:
    soup = BeautifulSoup(html, "html.parser")
    records: list[dict[str, object]] = []
    for index, table in enumerate(soup.find_all("table")):
        text = normalize_space(table.get_text(" ", strip=True))
        records.append(
            {
                "table_index": index,
                "characters": len(text),
                "mentions_mic": bool(re.search(r"\bMIC\b|minimum inhibitory", text, re.I)),
                "mentions_accession": bool(re.search(r"accession|BioProject|SRA|SAMN|PRJNA", text, re.I)),
                "mentions_biocide": bool(re.search(r"chlorhexidine|benzalkonium|biocide|triclosan|DDAC|QAC", text, re.I)),
                "preview": text[:500],
            }
        )
    return records


def term_hits(text: str, terms: list[str]) -> dict[str, bool]:
    lowered = text.lower()
    return {term: term.lower() in lowered for term in terms}


def ncbi_url(endpoint: str, params: dict[str, str | int]) -> str:
    base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi"
    common = {"tool": NCBI_TOOL, "email": NCBI_EMAIL, "retmode": "json"}
    return f"{base}?{urlencode({**params, **common})}"


def fetch_json(url: str, timeout: int = 60) -> dict[str, Any]:
    request = Request(url, headers={"User-Agent": "disinfectant-tolerance-dnabert/0.1"})
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def ncbi_count(db: str, term: str) -> int | None:
    if not term:
        return None
    try:
        payload = fetch_json(ncbi_url("esearch", {"db": db, "term": term, "retmax": 0}))
        return int(payload["esearchresult"]["count"])
    except Exception:
        return None


def extract_bioprojects(value: str) -> list[str]:
    return sorted(set(re.findall(r"PRJNA\d+", value or "")))


def build_candidate_records(raw_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates: list[dict[str, object]] = []
    supplements: list[dict[str, object]] = []
    tables: list[dict[str, object]] = []

    for candidate in IMMEDIATE_CANDIDATES:
        html = cached_text(raw_dir, candidate["study_key"], candidate["url"])
        fetch_blocked = bool(re.search(r"recaptcha|challengepage|unusual traffic", html, re.I))
        article_text = normalize_space(BeautifulSoup(html, "html.parser").get_text(" ", strip=True))
        supp_links = find_supplement_links(html, candidate["url"])
        table_rows = table_inventory(html)
        phenotype_hits = term_hits(article_text, candidate["phenotype_terms"])
        accession_hits = term_hits(article_text, candidate["accession_terms"])
        bioprojects = extract_bioprojects(str(candidate["reported_project"]) + " " + article_text)

        for link in supp_links:
            supplements.append({"study_key": candidate["study_key"], **link})
        for row in table_rows:
            tables.append({"study_key": candidate["study_key"], **row})

        sra_counts = {project: ncbi_count("sra", project) for project in bioprojects}
        assembly_counts = {project: ncbi_count("assembly", project) for project in bioprojects}
        time.sleep(0.34)

        machine_readable_supplements = sum(bool(row["machine_readable_candidate"]) for row in supp_links)
        spreadsheet_supplements = sum(bool(row["spreadsheet_candidate"]) for row in supp_links)
        phenotype_signal = any(phenotype_hits.values())
        accession_signal = any(accession_hits.values()) or bool(bioprojects)
        has_public_sequence_signal = bool(bioprojects) or any((value or 0) > 0 for value in sra_counts.values())
        can_build_without_manual_extraction = bool(
            phenotype_signal
            and accession_signal
            and has_public_sequence_signal
            and spreadsheet_supplements > 0
        )

        if fetch_blocked:
            verification_status = "fetch_blocked_by_provider"
            next_action = "Retry with manual browser/download or alternate source before making a readiness decision."
            can_build_without_manual_extraction = False
        elif can_build_without_manual_extraction:
            verification_status = "candidate_for_manifest_extraction"
            next_action = "Download machine-readable supplement and join row-level phenotype rows to public sequence accessions."
        elif machine_readable_supplements > 0 and not has_public_sequence_signal:
            verification_status = "phenotype_table_possible_but_sequence_link_unverified"
            next_action = "Inspect supplement columns, then query exact isolate names in NCBI Assembly/SRA."
        elif has_public_sequence_signal and machine_readable_supplements == 0:
            verification_status = "sequence_public_but_row_level_table_not_found"
            next_action = "Do not model yet; recover table from supplement, author data, or documented manual extraction."
        else:
            verification_status = "not_ready_from_automated_public_metadata"
            next_action = "Use as context unless row-level supplement and accession mapping are found manually."

        candidates.append(
            {
                "study_key": candidate["study_key"],
                "citation_short": candidate["citation_short"],
                "organism_group": candidate["organism_group"],
                "doi": candidate["doi"],
                "url": candidate["url"],
                "reported_endpoint": candidate["reported_endpoint"],
                "reported_project": candidate["reported_project"],
                "bioprojects_detected": ";".join(bioprojects),
                "sra_counts_by_project": json.dumps(sra_counts, sort_keys=True),
                "assembly_counts_by_project": json.dumps(assembly_counts, sort_keys=True),
                "article_tables_detected": len(table_rows),
                "supplement_links_detected": len(supp_links),
                "machine_readable_supplement_links": machine_readable_supplements,
                "spreadsheet_supplement_links": spreadsheet_supplements,
                "phenotype_terms_found": json.dumps(phenotype_hits, sort_keys=True),
                "accession_terms_found": json.dumps(accession_hits, sort_keys=True),
                "phenotype_signal": phenotype_signal,
                "accession_signal": accession_signal,
                "public_sequence_signal": has_public_sequence_signal,
                "fetch_blocked": fetch_blocked,
                "can_build_without_manual_extraction": can_build_without_manual_extraction,
                "verification_status": verification_status,
                "expected_status": candidate["expected_status"],
                "next_action": next_action,
            }
        )

    return pd.DataFrame(candidates), pd.DataFrame(supplements), pd.DataFrame(tables)


def build_report(candidates: pd.DataFrame) -> dict[str, object]:
    manifest_ready = candidates[candidates["can_build_without_manual_extraction"]].copy()
    sequence_public_no_table = candidates[candidates["verification_status"].eq("sequence_public_but_row_level_table_not_found")]
    return {
        "candidates_checked": int(len(candidates)),
        "candidate_manifest_extraction_ready": int(len(manifest_ready)),
        "sequence_public_but_row_level_table_not_found": int(len(sequence_public_no_table)),
        "top_next_candidates": candidates["study_key"].tolist(),
        "decision": (
            "extract_manifest_for_ready_candidate"
            if len(manifest_ready)
            else "no_immediate_candidate_manifest_ready_from_automated_metadata"
        ),
        "recommended_next_action": (
            "Proceed to row-level manifest extraction for: " + ", ".join(manifest_ready["study_key"].tolist())
            if len(manifest_ready)
            else "Do not start model training. The next step is manual supplement inspection/author-data recovery for the strongest near misses, especially Proteus and Enterococcus."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw/immediate_multiorganism_candidate_pages")
    parser.add_argument("--candidate-output", default="data/manifests/immediate_multiorganism_candidate_verification.csv")
    parser.add_argument("--supplement-output", default="data/manifests/immediate_multiorganism_supplement_links.csv")
    parser.add_argument("--table-output", default="data/manifests/immediate_multiorganism_article_table_inventory.csv")
    parser.add_argument("--report-output", default="reports/immediate_multiorganism_candidate_verification.json")
    args = parser.parse_args()

    candidates, supplements, tables = build_candidate_records(Path(args.raw_dir))
    report = build_report(candidates)

    Path(args.candidate_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report_output).parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.candidate_output, index=False)
    supplements.to_csv(args.supplement_output, index=False)
    tables.to_csv(args.table_output, index=False)
    Path(args.report_output).write_text(json.dumps(report, indent=2) + "\n")

    print(f"Wrote {args.candidate_output}")
    print(f"Wrote {args.supplement_output}")
    print(f"Wrote {args.table_output}")
    print(f"Wrote {args.report_output}")
    print(report["decision"])


if __name__ == "__main__":
    main()
