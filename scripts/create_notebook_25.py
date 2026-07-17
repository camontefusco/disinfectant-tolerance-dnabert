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
        """# 25 - *Acinetobacter* and *Klebsiella* phenotype-genome audit

This notebook checks whether another organism can become a supervised biocide-tolerance benchmark.

Strict rule:

- row-level biocide phenotype required;
- public genome/accession linkage required;
- aggregate MIC distributions are context only;
- no labels inferred from figures, distributions, or group summaries.
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

TABLES = PROJECT_ROOT / "data" / "manifests" / "acinetobacter_klebsiella_table_inventory.csv"
AB_ROWS = PROJECT_ROOT / "data" / "manifests" / "acinetobacter_2017_row_level_biocide_mic_rows.csv"
STUDIES = PROJECT_ROOT / "data" / "manifests" / "acinetobacter_klebsiella_study_readiness_audit.csv"
REPORT = PROJECT_ROOT / "reports" / "acinetobacter_klebsiella_phenotype_genome_audit.json"
"""
    ),
    markdown("## 1. Run the audit"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "audit_acinetobacter_klebsiella_phenotype_genome.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)
json.loads(REPORT.read_text())
"""
    ),
    markdown("## 2. Study readiness"),
    code(
        """study_audit = pd.read_csv(STUDIES)
display(study_audit[[
    "study_key", "organism", "row_level_biocide_phenotypes_found",
    "public_genome_linkage", "supervised_ml_ready", "reason", "doi", "url"
]])
"""
    ),
    markdown("## 3. Extracted row-level *A. baumannii* 2017 MIC rows"),
    code(
        """ab_rows = pd.read_csv(AB_ROWS)
print("Rows:", len(ab_rows))
if len(ab_rows):
    print("Unique isolates:", ab_rows["isolate_id"].nunique())
    display(pd.crosstab(ab_rows["compound"], ab_rows["endpoint"]))
    display(ab_rows.head(30))
"""
    ),
    markdown("## 4. Article table inventory"),
    code(
        """tables = pd.read_csv(TABLES)
display(tables[[
    "study_key", "table_index", "rows", "max_columns",
    "mentions_biocide", "mentions_accession", "looks_row_level", "caption"
]])
"""
    ),
    markdown(
        """## 5. Decision

No additional organism is ready for supervised modeling from these public materials.

Best near misses:

- *A. baumannii* 2017 has row-level biocide MICs for 47 isolates, but no public WGS/assembly linkage.
- *K. pneumoniae* aged-care 2024 has public WGS BioProject `PRJNA949397`, but biocide MICs are aggregate distribution counts.
- *A. baumannii* 2025 may have useful supplementary XLSX files, but row-level extraction needs a separate supplement-download workflow and genome linkage audit.

Use these papers as literature support unless row-level phenotype-genome tables are recovered.
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

output = ROOT / "notebooks" / "25_acinetobacter_klebsiella_phenotype_genome_audit.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
