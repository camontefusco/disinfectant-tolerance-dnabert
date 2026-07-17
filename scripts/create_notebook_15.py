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
    markdown("""# 15 - Audit MEGARes v3 and run a homology-aware fragment baseline

This notebook starts the exploratory pharmaceutical-water extension. It downloads the official MEGARes v3 nucleotide FASTA and annotations, isolates biocide-, metal-, and multi-compound-resistance sequences, and runs a first fragment-family baseline with entire MEGARes gene groups held out.

The output is **exploratory sequence-risk prioritization**, not validated isolate-level disinfectant-tolerance prediction.
"""),
    markdown("""## 1. Run the complete compact workflow

Use **Run All**. The official MEGARes files total approximately 11 MB. This notebook does not download pharmaceutical-water organism genomes and does not run DNABERT2 yet.
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
AUDIT_REPORT = PROJECT_ROOT / "reports" / "megares_v3_extension_audit.json"
BASELINE_REPORT = PROJECT_ROOT / "reports" / "megares_homology_aware_kmer_baseline.json"
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "megares_biocide_metal_manifest.csv"

audit_command = [sys.executable, str(PROJECT_ROOT / "scripts" / "audit_megares_v3_extension.py")]
baseline_command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "train_megares_homology_aware_baseline.py"),
    "--manifest", str(MANIFEST), "--output", str(BASELINE_REPORT),
]
print(" ".join(audit_command))
subprocess.run(audit_command, cwd=PROJECT_ROOT, check=True)
print(" ".join(baseline_command))
subprocess.run(baseline_command, cwd=PROJECT_ROOT, check=True)
audit = json.loads(AUDIT_REPORT.read_text())
baseline = json.loads(BASELINE_REPORT.read_text())
"""),
    markdown("""## 2. Inspect the curated exploratory subset

The primary labels are MEGARes hierarchy categories. They indicate resistance-associated sequence families, not measured behavior in a manufacturing-water system.
"""),
    code("""subset = audit["exploratory_subset"]
print("Subset fragments:", subset["count"])
print("Mechanisms:", subset["mechanism_count"])
print("Gene groups:", subset["group_count"])
print("Sequence lengths:", subset["sequence_length"])
display(pd.Series(subset["type_counts"], name="fragments").to_frame())
display(pd.Series(subset["class_counts"], name="fragments").sort_values(ascending=False).to_frame())
"""),
    markdown("""## 3. Inspect the homology-aware baseline

The split holds out complete MEGARes gene groups. This is stricter than a random sequence split and reduces close-homolog leakage. A publication workflow should later add nucleotide-similarity clustering before splitting.
"""),
    code("""print("Train fragments:", baseline["train_fragments"])
print("Test fragments:", baseline["test_fragments"])
print("Train groups:", baseline["train_groups"])
print("Test groups:", baseline["test_groups"])
print("Balanced accuracy:", round(baseline["balanced_accuracy"], 3))
display(pd.DataFrame(baseline["classification_report"]).T.round(3))
"""),
    markdown("""## 4. Scientific decision

Proceed to DNABERT2 only if:

- the MEGARes subset contains enough biocide and metal fragments for meaningful grouped evaluation;
- the homology-aware baseline is reproducible;
- the next claim remains fragment-level sequence-risk prioritization;
- pharmaceutical-water genomes are screened only after the classifier and interpretable database search are benchmarked.

The next bundled notebook should add similarity-cluster splitting, fragment-length/noise experiments, calibration, and DNABERT2-versus-k-mer comparison before any organism-level application.
"""),
]

notebook = {
    "cells": cells,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
output = ROOT / "notebooks" / "15_audit_megares_and_run_homology_aware_fragment_baseline.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
