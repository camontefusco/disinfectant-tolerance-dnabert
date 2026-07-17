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
        """# 24 - Bcc combined biocide phenotype-genome audit

This notebook tests whether multiple Bcc papers can form a second supervised benchmark closer to pharmaceutical-water and water-based-product risk.

Strict rule:

- row-level strain phenotype required;
- public genome match required;
- aggregate paper tables are context only;
- no inferred labels from species means, figures, or narrative text.

Primary row-level sources:

- Moore et al. 2009, DOI `10.1093/jac/dkn540`
- Kim/Wang et al. 2015, DOI `10.1007/s10295-015-1605-x`

Context source:

- Rushton et al. 2013, DOI `10.1128/AAC.00140-13`
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

PHENOTYPES = PROJECT_ROOT / "data" / "manifests" / "bcc_combined_biocide_phenotype_rows.csv"
MATCHES = PROJECT_ROOT / "data" / "manifests" / "bcc_ncbi_assembly_match_candidates.csv"
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "bcc_candidate_supervised_manifest.csv"
REPORT = PROJECT_ROOT / "reports" / "bcc_combined_phenotype_genome_audit.json"
"""
    ),
    markdown("## 1. Extract row-level phenotypes and query NCBI Assembly"),
    code(
        """command = [sys.executable, str(PROJECT_ROOT / "scripts" / "audit_bcc_combined_phenotype_genome.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)

report = json.loads(REPORT.read_text())
report
"""
    ),
    markdown("## 2. Extracted phenotype rows"),
    code(
        """phenotypes = pd.read_csv(PHENOTYPES)
print("Rows:", len(phenotypes))
print("Unique strains:", phenotypes["strain_id"].nunique())
display(pd.crosstab([phenotypes["study_key"], phenotypes["compound"]], phenotypes["endpoint"]))
display(phenotypes.head(30))
"""
    ),
    markdown("## 3. NCBI assembly match candidates"),
    code(
        """matches = pd.read_csv(MATCHES)
print("Candidate rows:", len(matches))
if len(matches):
    print("Exact/contains matched query strains:", matches[matches["exact_norm_match"] | matches["contains_norm_match"]]["query_strain"].nunique())
    display(matches[[
        "query_strain", "assembly_accession", "organism", "biosample_accession",
        "assembly_strain", "exact_norm_match", "contains_norm_match", "ftp_path_refseq"
    ]].head(50))
else:
    print("No assembly matches returned.")
"""
    ),
    markdown("## 4. Candidate supervised manifest"),
    code(
        """manifest = pd.read_csv(MANIFEST)
print("Manifest rows:", len(manifest))
if len(manifest):
    print("Unique matched strains:", manifest["strain_id"].nunique())
    display(pd.crosstab(manifest["compound"], manifest["endpoint"]))
    display(manifest.head(50))
else:
    print("No supervised manifest yet.")
"""
    ),
    markdown(
        """## 5. Decision

If at least 20 strain-level phenotype records have credible public genome matches, this can become a small Bcc benchmark. Otherwise it remains biological support for the pharmaceutical-water prioritization paper.

Before modeling, manually audit every genome match. NCBI strain-name search can return related or ambiguous assemblies, especially for collection strain aliases.
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

output = ROOT / "notebooks" / "24_bcc_combined_phenotype_genome_audit.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
