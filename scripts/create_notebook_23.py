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
        """# 23 - Verify Pottier 2023 *Pseudomonas aeruginosa* DDAC manifest feasibility

This notebook tests whether Pottier et al. 2023 can become a second supervised phenotype benchmark.

It does **not** create labels from plots or aggregate statements. A supervised manifest is allowed only if public row-level DDAC phenotype data can be joined to public genome assemblies by strain/sample ID.

Reference:

- Pottier et al. 2023, *Scientific Reports*, DOI `10.1038/s41598-023-29590-0`
- Figshare collection `10.6084/m9.figshare.c.6293046.v6`
- NCBI BioProject `PRJNA884650`
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

FIGSHARE = PROJECT_ROOT / "data" / "manifests" / "pottier_2023_figshare_inventory.csv"
ASSEMBLIES = PROJECT_ROOT / "data" / "manifests" / "pottier_2023_ncbi_assembly_inventory.csv"
DECISION = PROJECT_ROOT / "reports" / "pottier_2023_ddac_manifest_verification.json"
"""
    ),
    markdown("## 1. Query figshare and NCBI Assembly"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "verify_pottier_2023_ddac_manifest.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)

decision = json.loads(DECISION.read_text())
decision
"""
    ),
    markdown("## 2. Figshare file inventory"),
    code(
        """figshare = pd.read_csv(FIGSHARE)
display(figshare[[
    "title", "doi", "file_name", "file_mimetype",
    "machine_readable_table_file", "mentions_ddac", "download_url"
]])
"""
    ),
    markdown("## 3. NCBI Assembly inventory"),
    code(
        """assemblies = pd.read_csv(ASSEMBLIES)
print("Assemblies:", len(assemblies))
print("Unique strains:", assemblies["strain"].nunique())
display(assemblies[[
    "strain", "assembly_accession", "genbank_accession", "refseq_accession",
    "biosample_accession", "coverage", "assembly_status", "ftp_path_refseq"
]].head(20))
"""
    ),
    markdown(
        """## 4. Decision

If the decision is `blocked_missing_machine_readable_row_level_ddac_labels`, do not train a *P. aeruginosa* DDAC model yet.

The assemblies are usable, but the phenotype labels need a defensible extraction route:

1. author-provided table;
2. machine-readable supplement found elsewhere;
3. documented manual/OCR extraction from figure/table images with independent review.

Until then, Pottier et al. 2023 remains a high-priority candidate, not a completed second benchmark.
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

output = ROOT / "notebooks" / "23_verify_pottier_2023_ddac_manifest.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
