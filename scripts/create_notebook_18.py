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
    markdown("""# 18 - Plan a disk-capped pharmaceutical-water genome pilot

This notebook queries the official NCBI Datasets API for a small RefSeq pilot across four pharmaceutical-water-relevant taxa. It writes a reviewable metadata manifest and footprint estimate only.

It does **not** download assemblies and does **not** claim disinfectant-tolerance prediction.
"""),
    markdown("""## 1. Configuration

Use **Run All**. The default plan selects up to five current RefSeq assemblies per taxon. Change `PER_TAXON` only after reviewing disk availability.
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
PER_TAXON = 5
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "pharma_water_ncbi_pilot_plan.csv"
REPORT = PROJECT_ROOT / "reports" / "pharma_water_ncbi_pilot_plan.json"
"""),
    markdown("""## 2. Query NCBI metadata and build the pilot plan

The selection prefers complete genomes and chromosomes, then lower CheckM contamination and higher completeness when those metadata are available. The manifest leaves `download_approved` set to `False`.
"""),
    code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "plan_pharma_water_ncbi_pilot.py"),
    "--per-taxon", str(PER_TAXON), "--manifest", str(MANIFEST), "--report", str(REPORT),
]
print(" ".join(command))
if not REPORT.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(REPORT.read_text())
plan = pd.read_csv(MANIFEST)
display(pd.Series(report["selection"], name="value").to_frame())
display(plan)
"""),
    markdown("""## 3. Review gate

Review the accessions, organism names, assembly levels, and projected footprint before downloading anything.

The follow-on scoring notebook should download only approved assemblies and keep three channels side by side:

1. MEGARes, BacMet, and AMRFinderPlus matches;
2. calibrated exploratory scores for full and partial fragments;
3. organism, assembly, contig, and region context.
"""),
]

notebook = {
    "cells": cells,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
output = ROOT / "notebooks" / "18_plan_disk_capped_pharma_water_genome_pilot.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
