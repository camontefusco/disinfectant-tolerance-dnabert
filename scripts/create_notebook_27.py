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
        """# 27 - Verify immediate multi-organism candidates

This notebook verifies the four highest-priority non-*Listeria* candidates identified in notebook 26.

It checks only small public metadata:

- article tables;
- supplementary links;
- machine-readable spreadsheet candidates;
- BioProject/SRA/Assembly signals.

It does not download raw reads or assemblies.
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

CANDIDATES = PROJECT_ROOT / "data" / "manifests" / "immediate_multiorganism_candidate_verification.csv"
SUPPLEMENTS = PROJECT_ROOT / "data" / "manifests" / "immediate_multiorganism_supplement_links.csv"
TABLES = PROJECT_ROOT / "data" / "manifests" / "immediate_multiorganism_article_table_inventory.csv"
REPORT = PROJECT_ROOT / "reports" / "immediate_multiorganism_candidate_verification.json"
"""
    ),
    markdown("## 1. Run public metadata verifier"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "verify_immediate_multiorganism_candidates.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)
json.loads(REPORT.read_text())
"""
    ),
    markdown("## 2. Candidate readiness"),
    code(
        """candidate_audit = pd.read_csv(CANDIDATES)
display(candidate_audit[[
    "study_key", "organism_group", "bioprojects_detected",
    "supplement_links_detected", "machine_readable_supplement_links",
    "spreadsheet_supplement_links", "phenotype_signal", "accession_signal",
    "public_sequence_signal", "fetch_blocked",
    "can_build_without_manual_extraction", "verification_status", "next_action"
]])
"""
    ),
    markdown("## 3. Supplement links"),
    code(
        """supplements = pd.read_csv(SUPPLEMENTS)
display(supplements[[
    "study_key", "link_text", "href", "extension",
    "machine_readable_candidate", "spreadsheet_candidate"
]].sort_values(["study_key", "spreadsheet_candidate"], ascending=[True, False]))
"""
    ),
    markdown("## 4. Article table inventory"),
    code(
        """tables = pd.read_csv(TABLES)
display(tables[[
    "study_key", "table_index", "mentions_mic", "mentions_accession",
    "mentions_biocide", "preview"
]])
"""
    ),
    markdown(
        """## 5. Decision

The automated public-metadata check identifies **Proteus mirabilis** as the only candidate ready to attempt row-level manifest extraction:

- public supplementary spreadsheet candidate: `mic-171-01580-s001.xlsx`;
- public SRA signal: `PRJNA1154625` and `PRJNA475751`;
- article text explicitly states Table S1 contains metadata and accession numbers.

Near misses:

- *Enterococcus faecium*: public sequence signal, but only PDF supplement exposed by PMC metadata.
- *Enterococcus faecalis*: the PMC fetch was blocked by provider challenge in the automated run, so this remains unresolved, not rejected.
- Hospital-environment CHX preprint: public SRA signal and GitHub link, but no machine-readable Table S4 discovered by this metadata check.

Next step: try Proteus manifest extraction, but expect PMC's secure download page may require manual browser download of the XLSX.
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

output = ROOT / "notebooks" / "27_verify_immediate_multiorganism_candidates.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
