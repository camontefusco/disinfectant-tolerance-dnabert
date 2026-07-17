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
    markdown("""# 17 - Narrow the MEGARes task and audit BacMet evidence

This notebook asks a more useful exploratory question: can sequence models distinguish biocide-relevant MEGARes fragments from metal-only resistance fragments while holding out nucleotide-similarity clusters?

It also audits the official BacMet experimentally confirmed release. BacMet remains an interpretable external annotation layer, not an independent DNABERT nucleotide validation dataset.
"""),
    markdown("""## 1. Configuration

Use **Run All**. Existing similarity clusters, robustness conditions, and DNABERT2 embeddings from notebook 16 are reused. BacMet adds only small reference files.
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
CLUSTERED = PROJECT_ROOT / "data" / "manifests" / "megares_similarity_clustered_manifest.csv"
CONDITIONS = PROJECT_ROOT / "data" / "processed" / "megares_conditions"
EMBEDDINGS = PROJECT_ROOT / "data" / "processed" / "megares_embeddings"
BACMET_REPORT = PROJECT_ROOT / "reports" / "bacmet2_experimental_annotation_audit.json"
BINARY_REPORT = PROJECT_ROOT / "reports" / "megares_binary_biocide_relevance_benchmark.json"
"""),
    markdown("""## 2. Audit the official BacMet experimentally confirmed release

The BacMet release is protein FASTA plus a compound-mapping table. That is useful biological evidence, but it cannot be passed directly into a nucleotide DNABERT model. Because BacMet-derived references have also informed MEGARes, this notebook records it as an annotation audit rather than independent validation.
"""),
    code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "audit_bacmet_experimental.py"),
    "--report", str(BACMET_REPORT),
]
print(" ".join(command))
if not BACMET_REPORT.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
bacmet = json.loads(BACMET_REPORT.read_text())
display(pd.Series(bacmet["inventory"], name="value").to_frame())
print(bacmet["external_validation_gate"]["status"])
print(bacmet["external_validation_gate"]["reason"])
"""),
    markdown("""## 3. Evaluate a narrower binary task across five similarity-held splits

Positive fragments are MEGARes `Biocides` plus multi-compound classes explicitly containing biocide resistance. Negative fragments are `Metals` only. The four drug-and-metal records without a biocide annotation are excluded.

The benchmark reuses the full, centered-500-bp, and deterministic-1%-substitution representations from notebook 16. It reports balanced accuracy, ROC AUC, average precision, F1, log loss, Brier score, and calibration error.
"""),
    code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "evaluate_megares_binary_biocide_relevance.py"),
    "--manifest", str(CLUSTERED), "--conditions-dir", str(CONDITIONS),
    "--embeddings-dir", str(EMBEDDINGS), "--output", str(BINARY_REPORT),
]
print(" ".join(command))
if not BINARY_REPORT.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(BINARY_REPORT.read_text())
print(report["task"])
summary = pd.DataFrame(report["summary"])
display(summary.sort_values(["condition", "calibrated_balanced_accuracy_mean"], ascending=[True, False]).round(3))
"""),
    markdown("""## 4. Scientific decision

Use this binary benchmark to decide whether the exploratory model has enough similarity-held signal to screen pharmaceutical-water genomes.

The next stage should keep three channels visible beside each other:

1. interpretable MEGARes, BacMet, and AMRFinderPlus matches;
2. calibrated exploratory fragment scores from the selected representation;
3. organism and contig context for manual review.

Do not call the resulting ranking disinfectant-tolerance prediction. BacMet nucleotide mapping can be added later as a reviewed annotation-enrichment step, but it is not required before a small organism-screening pilot.
"""),
]

notebook = {
    "cells": cells,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
output = ROOT / "notebooks" / "17_narrow_megares_biocide_task_and_audit_bacmet.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
