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
        """# 26 - Multi-organism biocide dataset readiness audit

This notebook asks a narrow question:

Can any organism beyond *Listeria monocytogenes* support a second supervised biocide-tolerance benchmark using only public data?

Strict inclusion rule:

- row-level biocide phenotype labels are required;
- isolate IDs must map to public genome/SRA/assembly accessions;
- aggregate MIC distributions, figure-only data, and experimental-evolution labels are not enough for training;
- small datasets may still be useful as external case studies.
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

AUDIT_CSV = PROJECT_ROOT / "data" / "manifests" / "multiorganism_biocide_ml_readiness_audit.csv"
REPORT = PROJECT_ROOT / "reports" / "multiorganism_biocide_ml_readiness_audit.json"
"""
    ),
    markdown("## 1. Run the audit"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "audit_multiorganism_ml_ready_biocide_datasets.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(REPORT.read_text())
report["summary"]
"""
    ),
    markdown("## 2. Ranked candidates"),
    code(
        """audit = pd.read_csv(AUDIT_CSV)
display(audit[[
    "priority", "study_key", "organism_group", "endpoint", "reported_n_phenotyped",
    "reported_n_sequenced", "row_level_labels_status", "public_genome_status",
    "accession_or_project", "ml_readiness", "audit_batch", "doi"
]])
"""
    ),
    markdown("## 3. Immediate audit batch"),
    code(
        """immediate = audit[audit["audit_batch"].eq("immediate")].copy()
display(immediate[[
    "study_key", "organism_group", "endpoint", "row_level_label_evidence",
    "genome_linkage_evidence", "next_action", "caution", "url"
]])
"""
    ),
    markdown("## 4. Small external-case candidates"),
    code(
        """external = audit[audit["audit_batch"].eq("external_case")].copy()
display(external[[
    "study_key", "organism_group", "endpoint", "reported_n_phenotyped",
    "accession_or_project", "next_action", "caution", "url"
]])
"""
    ),
    markdown("## 5. Context-only or mechanistic datasets"),
    code(
        """context = audit[audit["audit_batch"].eq("context")].copy()
display(context[[
    "study_key", "organism_group", "endpoint", "ml_readiness",
    "row_level_labels_status", "public_genome_status", "caution", "url"
]])
"""
    ),
    markdown(
        """## 6. Decision

The next bounded batch is:

1. *Proteus mirabilis* 2025 cationic-biocide dataset.
2. *Enterococcus faecium* 2019 chlorhexidine dataset.
3. *Enterococcus faecalis* 2022 chlorhexidine dataset.
4. Hartmann lab mixed hospital-environment chlorhexidine preprint.

Do not train yet. The next notebook should extract or verify supplements for this batch and produce one of two outcomes:

- a row-level phenotype-genome manifest for a second supervised benchmark; or
- a documented stop decision that no public non-*Listeria* organism is ML-ready without author contact or manual data recovery.

The strongest current paper remains the *Listeria* supervised benchmark plus a transparent, evidence-based audit showing why pharma-relevant organisms are exploratory rather than validated tolerance predictors.
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

output = ROOT / "notebooks" / "26_audit_multiorganism_ml_ready_biocide_datasets.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
