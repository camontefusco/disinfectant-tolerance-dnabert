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
    markdown("""# 19 - Download and screen the disk-capped pharmaceutical-water pilot

This notebook downloads only the twenty reviewed RefSeq assemblies from notebook 18, finds interpretable MEGARes nucleotide matches, adds deterministic background windows, and scores reviewable 500 bp regions with the calibrated partial-fragment model.

The output is **exploratory genomic risk prioritization**, not disinfectant-tolerance prediction.
"""),
    markdown("""## 1. Configuration

Use **Run All**. Existing assemblies and reports are reused. The downloader enforces a `200 MB` estimated FASTA ceiling and discards temporary ZIP packages after extracting genomic FASTA files.
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
PLAN = PROJECT_ROOT / "data" / "manifests" / "pharma_water_ncbi_pilot_plan.csv"
ASSEMBLIES = PROJECT_ROOT / "data" / "manifests" / "pharma_water_downloaded_assemblies.csv"
FRAGMENTS = PROJECT_ROOT / "data" / "manifests" / "pharma_water_pilot_scored_fragments.csv"
REPORT = PROJECT_ROOT / "reports" / "pharma_water_pilot_screen.json"
"""),
    markdown("""## 2. Download compact NCBI FASTA packages

The download is already bounded by the reviewed manifest. Each NCBI ZIP is held in memory only long enough to extract its genomic FASTA.
"""),
    code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "download_pharma_water_ncbi_pilot.py"),
    "--plan", str(PLAN), "--manifest", str(ASSEMBLIES), "--execute",
]
print(" ".join(command))
if not ASSEMBLIES.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
assemblies = pd.read_csv(ASSEMBLIES)
print("Assemblies:", len(assemblies))
print("FASTA MB:", round(assemblies["fasta_bytes"].sum() / 1_000_000, 2))
"""),
    markdown("""## 3. Screen MEGARes hits and calculate calibrated exploratory scores

MEGARes nucleotide references are aligned to each assembly with BLAST+. Candidate windows are centered on matches and supplemented with deterministic background windows. Scores use the notebook 17 partial-fragment winner: a five-seed calibrated k-mer random-forest ensemble.
"""),
    code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "screen_pharma_water_pilot.py"),
    "--assemblies", str(ASSEMBLIES), "--fragments", str(FRAGMENTS), "--report", str(REPORT),
]
print(" ".join(command))
if not REPORT.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(REPORT.read_text())
fragments = pd.read_csv(FRAGMENTS)
display(pd.Series(report, name="value").to_frame())
display(fragments.sort_values("exploratory_biocide_relevance_score", ascending=False).head(30))
"""),
    markdown("""## 4. Interpretation gate

Review the top-ranked regions and MEGARes matches by organism. Scores prioritize fragments for inspection; they do not predict whether an isolate survives a disinfectant exposure.

The next hardening step is AMRFinderPlus annotation and a reviewed protein-aware BacMet mapping layer beside these regions. DNABERT2 can then be applied selectively to a small set of intact candidate regions, rather than indiscriminately to every genomic window.
"""),
]
notebook = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}}, "nbformat": 4, "nbformat_minor": 5}
output = ROOT / "notebooks" / "19_download_and_screen_pharma_water_pilot.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
