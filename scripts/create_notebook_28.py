#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def markdown(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)}


def code(text: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": text.splitlines(keepends=True)}


cells = [
    markdown(
        """# 28 - Extract the Proteus 2025 candidate manifest

This notebook attempts the first real second-organism manifest extraction.

Candidate: Bennett et al. 2025, *Proteus mirabilis* chlorhexidine/cationic-biocide study.

Why this candidate:

- the article states Table S1 contains clinical metadata and accession numbers;
- public sequence projects are reported: `PRJNA1154625`, `PRJNA475751`, and `PRJNA12624`;
- PMC exposes `mic-171-01580-s001.xlsx` as a supplementary spreadsheet.

Important: direct command-line download of the PMC XLSX may return a secure "Preparing to download" HTML page. If that happens, download the XLSX in a browser and place it at `data/raw/proteus_2025/mic-171-01580-s001.xlsx` or leave it in `~/Downloads/`.
"""
    ),
    code(
        """from __future__ import annotations
from pathlib import Path
import json
import subprocess
import sys
import pandas as pd

PROJECT_ROOT = Path.cwd()
if PROJECT_ROOT.name == "notebooks":
    PROJECT_ROOT = PROJECT_ROOT.parent

SHEET_INVENTORY = PROJECT_ROOT / "data" / "manifests" / "proteus_2025_table_s1_sheet_inventory.csv"
ISOLATE_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "proteus_2025_candidate_isolate_manifest.csv"
PHENOTYPE_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "proteus_2025_candidate_phenotype_manifest.csv"
REPORT = PROJECT_ROOT / "reports" / "proteus_2025_manifest_extraction.json"
"""
    ),
    markdown("## 1. Run extraction"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "extract_proteus_2025_manifest.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(REPORT.read_text())
report
"""
    ),
    markdown("## 2. Workbook sheet inventory"),
    code(
        """if SHEET_INVENTORY.exists() and SHEET_INVENTORY.stat().st_size:
    try:
        sheet_inventory = pd.read_csv(SHEET_INVENTORY)
        display(sheet_inventory)
    except pd.errors.EmptyDataError:
        print("No sheet inventory yet: valid XLSX supplement not found.")
else:
    print("No sheet inventory yet: valid XLSX supplement not found.")
"""
    ),
    markdown("## 3. Candidate isolate manifest"),
    code(
        """if ISOLATE_MANIFEST.exists() and ISOLATE_MANIFEST.stat().st_size:
    try:
        isolates = pd.read_csv(ISOLATE_MANIFEST)
        print("Rows:", len(isolates))
        display(isolates.head(30))
    except pd.errors.EmptyDataError:
        print("No isolate manifest yet.")
else:
    print("No isolate manifest yet.")
"""
    ),
    markdown("## 4. Candidate phenotype manifest"),
    code(
        """if PHENOTYPE_MANIFEST.exists() and PHENOTYPE_MANIFEST.stat().st_size:
    try:
        phenotypes = pd.read_csv(PHENOTYPE_MANIFEST)
        print("Rows:", len(phenotypes))
        if len(phenotypes):
            display(pd.crosstab(phenotypes["compound_or_endpoint"], phenotypes["phenotype_type"]))
            display(phenotypes.head(40))
    except pd.errors.EmptyDataError:
        print("No phenotype manifest yet.")
else:
    print("No phenotype manifest yet.")
"""
    ),
    markdown(
        """## 5. Decision

Current expected outcomes:

- `candidate_manifest_extracted`: proceed to audit label balance, accession mappings, and whether enough independent isolates exist for lineage-aware evaluation.
- `blocked_manual_supplement_download_required`: download the Table S1 XLSX manually from the PMC article and rerun this notebook.
- `workbook_found_but_manifest_extraction_failed`: inspect the sheet inventory and update the column-mapping rules.

Do not download reads or train a Proteus model until this notebook produces a credible row-level phenotype-genome manifest.
"""
    ),
]

notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.10"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

output = ROOT / "notebooks" / "28_extract_proteus_2025_candidate_manifest.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
