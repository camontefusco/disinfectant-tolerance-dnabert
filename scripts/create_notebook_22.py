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
        """# 22 - Literature audit for pharma-relevant labeled biocide datasets

This notebook asks whether we can add a second phenotype-labeled benchmark closer to pharmaceutical-water organisms without fabricating labels.

Strict inclusion rule for supervised ML:

1. isolate-level disinfectant/biocide phenotype labels;
2. public genome accessions or recoverable assemblies;
3. row-level mapping between the phenotype table and the sequence records.

Mechanistic studies are valuable, but they are not converted into training labels unless the paper provides isolate-level phenotypes and genome access.
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

CSV = PROJECT_ROOT / "data" / "manifests" / "literature_phenotype_dataset_audit.csv"
REPORT = PROJECT_ROOT / "reports" / "literature_phenotype_dataset_audit.json"
"""
    ),
    markdown("## 1. Build the literature audit"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "audit_literature_phenotype_datasets.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)

report = json.loads(REPORT.read_text())
report["summary"]
"""
    ),
    markdown("## 2. Candidate studies ranked by readiness"),
    code(
        """audit = pd.read_csv(CSV)
display(audit[[
    "priority", "citation_short", "organism_group", "biocide_endpoint",
    "reported_n_phenotyped", "reported_n_sequenced", "supervised_ml_readiness",
    "accession_or_project", "doi"
]])
"""
    ),
    markdown("## 3. What can become the next labeled benchmark?"),
    code(
        """ready = audit[audit["supervised_ml_readiness"].eq("highest_priority_candidate")]
conditional = audit[audit["supervised_ml_readiness"].str.startswith("conditional", na=False)]
case_studies = audit[audit["supervised_ml_readiness"].isin(["case_study_only", "literature_context_only"])]

print("Highest-priority supervised candidate:")
display(ready[["citation_short", "organism_group", "biocide_endpoint", "recommended_next_action", "caution", "url"]])

print("Conditional candidates requiring strain-to-genome verification:")
display(conditional[["citation_short", "organism_group", "biocide_endpoint", "recommended_next_action", "caution", "url"]])

print("Case-study/context references, not training data:")
display(case_studies[["citation_short", "organism_group", "biocide_endpoint", "recommended_next_action", "caution", "url"]])
"""
    ),
    markdown(
        """## 4. Decision

Do not add unlabeled organisms to improve model claims. The only identified non-*Listeria* organism worth pursuing as a second supervised benchmark is the *Pseudomonas aeruginosa* DDAC study, because the literature reports both DDAC phenotypes and public sequence deposition.

The next notebook should therefore be narrow: extract the Pottier et al. 2023 phenotype table and verify row-level mapping to BioProject PRJNA884650. If that mapping fails, stop and keep the pharmaceutical-water stage exploratory.
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

output = ROOT / "notebooks" / "22_literature_phenotype_dataset_audit.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
