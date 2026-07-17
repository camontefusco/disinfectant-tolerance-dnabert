#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd
from bs4 import BeautifulSoup


STUDIES = {
    "ab_2017": {
        "pmcid": "PMC5622949",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5622949/",
        "citation_short": "Lin et al. 2017, Frontiers in Microbiology",
        "doi": "10.3389/fmicb.2017.01836",
        "organism": "Acinetobacter baumannii",
        "endpoint": "MICs for triclosan, chlorhexidine acetate, benzalkonium bromide, sodium hypochlorite, hydrogen peroxide, ethanol",
    },
    "ab_env_2024": {
        "pmcid": "PMC11097793",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11097793/",
        "citation_short": "Kheljan et al. 2024, Scientific Reports",
        "doi": "10.1038/s41598-024-61827-0",
        "organism": "Acinetobacter baumannii",
        "endpoint": "Aggregate MIC distribution for benzalkonium chloride, chlorhexidine, formaldehyde, triclosan",
    },
    "ab_source_2025": {
        "pmcid": "PMC12661710",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12661710/",
        "citation_short": "Gavini et al. 2025, Antimicrobial Resistance & Infection Control",
        "doi": "10.1186/s13756-025-01658-6",
        "organism": "Acinetobacter baumannii",
        "endpoint": "MBCs for sodium hypochlorite, ethanol, hydrogen peroxide, chlorhexidine, benzalkonium chloride",
    },
    "kp_agedcare_2024": {
        "pmcid": "PMC11051875",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11051875/",
        "citation_short": "Kleyn et al. 2024, Microorganisms",
        "doi": "10.3390/microorganisms12040751",
        "organism": "Klebsiella pneumoniae complex",
        "endpoint": "WGS BioProject PRJNA949397; aggregate CHG/TRI/BZK MIC distribution; row-level antibiotic MIC table",
    },
    "kp_crkp_2022": {
        "pmcid": "PMC9219816",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9219816/",
        "citation_short": "Li et al. 2022, Frontiers in Microbiology",
        "doi": "10.3389/fmicb.2022.931592",
        "organism": "Carbapenem-resistant Klebsiella pneumoniae",
        "endpoint": "MIC/MBC for seven disinfectants and PCR efflux genes; public WGS linkage not confirmed",
    },
}


def fetch_text(url: str, timeout: int = 60) -> str:
    request = Request(url, headers={"User-Agent": "disinfectant-tolerance-dnabert/0.1"})
    with urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="ignore")


def cached_html(raw_dir: Path, study_key: str) -> str:
    meta = STUDIES[study_key]
    path = raw_dir / f"{meta['pmcid']}.html"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(fetch_text(meta["url"]))
    return path.read_text(errors="ignore")


def table_rows(table) -> list[list[str]]:
    rows: list[list[str]] = []
    for tr in table.find_all("tr"):
        cells = [re.sub(r"\s+", " ", cell.get_text(" ", strip=True)) for cell in tr.find_all(["th", "td"])]
        if cells:
            rows.append(cells)
    return rows


def table_caption(table) -> str:
    caption = table.find_previous(class_="caption")
    if caption:
        return re.sub(r"\s+", " ", caption.get_text(" ", strip=True))
    return re.sub(r"\s+", " ", table.get_text(" ", strip=True)[:250])


def inventory_article_tables(raw_dir: Path) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for key, meta in STUDIES.items():
        soup = BeautifulSoup(cached_html(raw_dir, key), "html.parser")
        for index, table in enumerate(soup.find_all("table")):
            rows = table_rows(table)
            text = " ".join(" ".join(row) for row in rows[:6])
            records.append(
                {
                    "study_key": key,
                    "citation_short": meta["citation_short"],
                    "organism": meta["organism"],
                    "table_index": index,
                    "caption": table_caption(table),
                    "rows": len(rows),
                    "max_columns": max([len(row) for row in rows], default=0),
                    "mentions_biocide": bool(re.search(r"BZK|benzalkonium|chlorhexidine|CHG|triclosan|biocide|disinfect", text, re.I)),
                    "mentions_accession": bool(re.search(r"BioProject|SAMN|accession", text, re.I)),
                    "looks_row_level": bool(re.search(r"Isolate|Strain", text, re.I)) and len(rows) > 10,
                    "preview": text[:500],
                }
            )
    return pd.DataFrame(records)


def extract_ab_2017_biocide_rows(raw_dir: Path) -> pd.DataFrame:
    soup = BeautifulSoup(cached_html(raw_dir, "ab_2017"), "html.parser")
    tables = soup.find_all("table")
    if len(tables) < 4:
        return pd.DataFrame()
    rows = table_rows(tables[3])
    # Article Table 4: row-level isolates with triclosan and chlorhexidine acetate MICs
    # plus relative gene-expression values. BZK is only available in grouped Table 3.
    records: list[dict[str, object]] = []
    for row in rows:
        if not row or not re.match(r"AB\d+", row[0]):
            continue
        isolate = row[0]
        antibiotic_susceptibility = row[1] if len(row) > 1 else ""
        mic_values = row[2:4]
        for compound, value in zip(
            ["triclosan", "chlorhexidine acetate"],
            mic_values,
            strict=False,
        ):
            records.append(
                {
                    "study_key": "ab_2017",
                    "citation_short": STUDIES["ab_2017"]["citation_short"],
                    "organism": STUDIES["ab_2017"]["organism"],
                    "isolate_id": isolate,
                    "compound": compound,
                    "endpoint": "MIC",
                    "value": value,
                    "units": "ug/mL",
                    "antibiotic_susceptibility": antibiotic_susceptibility,
                    "genome_link_status": "not_confirmed",
                    "source_url": STUDIES["ab_2017"]["url"],
                    "doi": STUDIES["ab_2017"]["doi"],
                }
            )
    return pd.DataFrame(records)


def build_study_audit(table_inventory: pd.DataFrame, ab_rows: pd.DataFrame) -> pd.DataFrame:
    records = []
    for key, meta in STUDIES.items():
        subset = table_inventory[table_inventory["study_key"].eq(key)]
        row_level_biocide = False
        public_genomes = "not_confirmed"
        usable = "no"
        reason = ""
        if key == "ab_2017":
            row_level_biocide = not ab_rows.empty
            reason = "Row-level biocide MICs are in article Table 3, but public WGS/assembly accessions for AB01-AB47 are not reported."
        elif key == "kp_agedcare_2024":
            public_genomes = "BioProject PRJNA949397 reported"
            reason = "Public WGS accessions and row-level antibiotic MICs are reported, but CHG/TRI/BZK biocide MICs are aggregate distribution counts in Table 4."
        elif key == "ab_source_2025":
            reason = "Article reports 64 isolates and supplementary XLSX files, but downloaded supplement links require interstitial handling; row-level MBC extraction not completed in this automated audit."
        elif key == "ab_env_2024":
            reason = "Article table exposed through PMC is aggregate MIC distribution for 96 isolates, not row-level isolate labels."
        else:
            reason = "Public article evidence is phenotype/PCR-focused; row-level public WGS linkage not confirmed."
        records.append(
            {
                "study_key": key,
                "citation_short": meta["citation_short"],
                "doi": meta["doi"],
                "url": meta["url"],
                "organism": meta["organism"],
                "endpoint": meta["endpoint"],
                "article_tables": int(len(subset)),
                "row_level_biocide_phenotypes_found": row_level_biocide,
                "public_genome_linkage": public_genomes,
                "supervised_ml_ready": usable,
                "reason": reason,
            }
        )
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw/organism_audit_pmcs")
    parser.add_argument("--table-inventory", default="data/manifests/acinetobacter_klebsiella_table_inventory.csv")
    parser.add_argument("--phenotype-output", default="data/manifests/acinetobacter_2017_row_level_biocide_mic_rows.csv")
    parser.add_argument("--study-audit", default="data/manifests/acinetobacter_klebsiella_study_readiness_audit.csv")
    parser.add_argument("--report-output", default="reports/acinetobacter_klebsiella_phenotype_genome_audit.json")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    table_inventory = inventory_article_tables(raw_dir)
    ab_rows = extract_ab_2017_biocide_rows(raw_dir)
    study_audit = build_study_audit(table_inventory, ab_rows)

    report = {
        "studies_audited": int(len(study_audit)),
        "row_level_ab_2017_biocide_rows": int(len(ab_rows)),
        "row_level_ab_2017_isolates": int(ab_rows["isolate_id"].nunique()) if not ab_rows.empty else 0,
        "supervised_ml_ready_studies": int(study_audit["supervised_ml_ready"].eq("yes").sum()),
        "decision": "no_additional_supervised_organism_ready",
        "best_near_misses": [
            "A. baumannii 2017: row-level biocide MICs, but no public genome linkage.",
            "K. pneumoniae 2024 aged-care: public genomes, but biocide MICs are aggregate.",
            "A. baumannii 2025 source study: likely useful xlsx supplements, but automated download needs separate handling.",
        ],
        "recommended_next_action": "Do not start modeling. Use these papers as rationale/context unless author/raw supplemental row-level phenotype-genome data are recovered.",
    }

    Path(args.table_inventory).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report_output).parent.mkdir(parents=True, exist_ok=True)
    table_inventory.to_csv(args.table_inventory, index=False)
    ab_rows.to_csv(args.phenotype_output, index=False)
    study_audit.to_csv(args.study_audit, index=False)
    Path(args.report_output).write_text(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {args.table_inventory}")
    print(f"Wrote {args.phenotype_output}")
    print(f"Wrote {args.study_audit}")
    print(f"Wrote {args.report_output}")
    print(report["decision"])


if __name__ == "__main__":
    main()
