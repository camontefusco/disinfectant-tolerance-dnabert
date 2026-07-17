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
    markdown("""# 21 - Publication readiness and final stop rule

This notebook turns the completed benchmark and pharmaceutical-water pilot into a bounded publication dossier. It does not add organisms, labels, or model families.
"""),
    code("""from __future__ import annotations
from pathlib import Path
import subprocess
import sys
import pandas as pd

PROJECT_ROOT = Path.cwd()
if PROJECT_ROOT.name == "notebooks":
    PROJECT_ROOT = PROJECT_ROOT.parent
DOSSIER = PROJECT_ROOT / "reports" / "publication_readiness_dossier.md"
TOP = PROJECT_ROOT / "reports" / "publication_top_candidate_regions.csv"
BENCHMARK = PROJECT_ROOT / "reports" / "publication_benchmark_table.csv"
"""),
    markdown("""## 1. Build the dossier

The dossier defines the final publishable claim, stop rule, benchmark table, evidence enrichment, and strongest candidate regions.
"""),
    code("""command = [sys.executable, str(PROJECT_ROOT / "scripts" / "build_publication_dossier.py")]
print(" ".join(command))
subprocess.run(command, cwd=PROJECT_ROOT, check=True)
print(DOSSIER.read_text())
"""),
    markdown("""## 2. Candidate regions for the Results section"""),
    code("""display(pd.read_csv(TOP).head(40))"""),
    markdown("""## 3. Benchmark table for the Results section"""),
    code("""display(pd.read_csv(BENCHMARK).round(3))"""),
]

notebook = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}}, "nbformat": 4, "nbformat_minor": 5}
output = ROOT / "notebooks" / "21_publication_readiness_and_final_dossier.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
