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
    markdown("""# 14 - Audit the Harrand use-level sanitizer-survival dataset

This notebook is a feasibility gate for the harder sanitizer-survival extension. It downloads only small public audit files, inventories the cohort and assembly routes, checks overlap with the current project, and verifies whether isolate-level phenotype labels are publicly available. It deliberately stops before bulk genome download.
"""),
    markdown("""## 1. Run the complete audit

Use **Run All**. The cached public files total approximately 1.2 MB. No assembly ZIP or genome data are downloaded.
"""),
    code("""from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

PROJECT_ROOT = Path.cwd()
if PROJECT_ROOT.name == "notebooks":
    PROJECT_ROOT = PROJECT_ROOT.parent
REPORT = PROJECT_ROOT / "reports" / "harrand_2026_feasibility_audit.json"
command = [sys.executable, str(PROJECT_ROOT / "scripts" / "audit_harrand_2026_feasibility.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)
audit = json.loads(REPORT.read_text())
"""),
    markdown("""## 2. Cohort inventory

The paper reports measured 30-second exposure outcomes for BC at 300 ppm and PAA at 80 ppm across 501 isolates. A subset of 108 isolates was also exposed to NaOCl at 500 ppm. Table S1 should contain genome-accessibility metadata for the full cohort.
"""),
    code("""table_s1 = audit["supplementary_table_s1"]
print("Total isolates:", table_s1["isolates"])
print("L. monocytogenes:", table_s1["l_monocytogenes_isolates"])
display(pd.Series(table_s1["species_counts"], name="isolates").to_frame())
display(pd.Series(table_s1["lineage_counts"], name="isolates").to_frame())
display(pd.Series(table_s1["environment_counts"], name="isolates").to_frame())
"""),
    markdown("""## 3. Assembly availability and storage gate

The Cornell deposit contains a compact assembly ZIP for isolates not published elsewhere. This notebook records its size without downloading it.
"""),
    code("""sources = audit["public_sources"]
deposit = audit["cornell_deposit"]
print("Cornell assembly ZIP MB:", sources["cornell_assembly_zip_mb"])
print("Cornell metadata isolates:", deposit["metadata_isolates"])
print("Rows matching Cornell metadata:", deposit["s1_isolates_matching_cornell_metadata"])
print("Normalized overlap with current BC-MIC project:", audit["current_project_overlap"]["exact_normalized_isolate_overlap"])
"""),
    markdown("""## 4. Modeling gate

A public genome-accessibility table is not automatically a public training dataset. The required labels are isolate-level BC, PAA, and NaOCl log reductions, ideally with replicate-level values or clearly documented averages.
"""),
    code("""gate = audit["modeling_gate"]
print("Status:", gate["status"])
print("Reason:", gate["reason"])
print("Required next step:", gate["required_next_step"])
print("Phenotype-like columns found in public Table S1:", table_s1["phenotype_label_columns"])
"""),
    markdown("""## 5. Scientific decision

If the status is `blocked_missing_public_isolate_level_labels`, stop before genome download. Request the phenotype table from the corresponding author. The public article supports the research question, but a model cannot be trained honestly from figures, summary statistics, or synthetic labels.

The MEGARes-centered pharmaceutical-water extension remains feasible independently as **exploratory sequence-risk prioritization**, not as validated sanitizer-survival prediction.
"""),
]

notebook = {
    "cells": cells,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
output = ROOT / "notebooks" / "14_audit_harrand_use_level_sanitizer_survival_dataset.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
