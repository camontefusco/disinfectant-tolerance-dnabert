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
    markdown("""# 16 - Benchmark similarity-aware MEGARes fragment models and robustness

This notebook evaluates the exploratory resistance-associated fragment classifier. It builds nucleotide-similarity clusters, generates deterministic partial-fragment and substitution-noise conditions, extracts reusable DNABERT2 embeddings, and compares calibrated k-mer, random-forest, XGBoost, and DNABERT2 strategies.

The output is **fragment-level genomic risk prioritization**, not validated isolate-level disinfectant-tolerance prediction.
"""),
    markdown("""## 1. Configuration

Use **Run All**. Completed artifacts are reused automatically. The DNABERT2 embedding stage is the slow part, but it embeds only `1,163` fragments for each of three conditions.
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
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "megares_biocide_metal_manifest.csv"
CLUSTERED = PROJECT_ROOT / "data" / "manifests" / "megares_similarity_clustered_manifest.csv"
CONDITIONS = PROJECT_ROOT / "data" / "processed" / "megares_conditions"
EMBEDDINGS = PROJECT_ROOT / "data" / "processed" / "megares_embeddings"
LOCAL_MODEL = PROJECT_ROOT / "data" / "external" / "dnabert2_pytorch_fallback"
REPORT = PROJECT_ROOT / "reports" / "megares_similarity_robustness_benchmark.json"
EMBEDDINGS.mkdir(parents=True, exist_ok=True)
"""),
    markdown("""## 2. Build similarity clusters and robustness conditions

Similarity clusters are connected components at at least 90% nucleotide identity and 80% bidirectional coverage. Entire clusters remain held out during evaluation.
"""),
    code("""cluster_command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "cluster_megares_similarity.py"),
    "--manifest", str(MANIFEST), "--output", str(CLUSTERED),
]
condition_command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "prepare_megares_fragment_conditions.py"),
    "--manifest", str(CLUSTERED), "--output-dir", str(CONDITIONS),
]
if not CLUSTERED.exists():
    subprocess.run(cluster_command, cwd=PROJECT_ROOT, check=True)
if not (CONDITIONS / "noise_1%.csv").exists():
    subprocess.run(condition_command, cwd=PROJECT_ROOT, check=True)
clustered = pd.read_csv(CLUSTERED)
print("Fragments:", len(clustered))
print("Similarity clusters:", clustered["similarity_cluster"].nunique())
"""),
    markdown("""## 3. Extract reusable DNABERT2 embeddings

The three conditions are full sequence, centered 500 bp fragment, and deterministic 1% nucleotide substitutions. Existing embedding archives are reused.
"""),
    code("""for condition in ("full", "partial_500bp", "noise_1%"):
    output = EMBEDDINGS / f"{condition}.npz"
    command = [
        sys.executable, str(PROJECT_ROOT / "scripts" / "embed_dnabert2.py"),
        "--windows", str(CONDITIONS / f"{condition}.csv"),
        "--model", str(LOCAL_MODEL), "--maximum-windows-per-isolate", "1",
        "--output", str(output),
    ]
    if not output.exists():
        print("Running:", " ".join(command))
        subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    print(condition, "embedding MB:", round(output.stat().st_size / 1_000_000, 2))
"""),
    markdown("""## 4. Evaluate calibrated models across five repeated similarity-cluster splits

Calibration is learned from development fragments only. The test set remains held out. Metrics include balanced accuracy, macro F1, log loss, Brier score, and expected calibration error.
"""),
    code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "evaluate_megares_similarity_robustness.py"),
    "--manifest", str(CLUSTERED), "--conditions-dir", str(CONDITIONS),
    "--embeddings-dir", str(EMBEDDINGS), "--output", str(REPORT),
]
print(" ".join(command))
if not REPORT.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(REPORT.read_text())
summary = pd.DataFrame(report["summary"])
display(summary.sort_values(["condition", "calibrated_balanced_accuracy_mean"], ascending=[True, False]).round(3))
"""),
    markdown("""## 5. Scientific decision

Select a model based on similarity-held-out performance, robustness, and calibration together. Do not select on clean balanced accuracy alone.

Proceed to pharmaceutical-water genome screening only after:

- the DNABERT2 comparison is interpreted against the simpler sequence baselines;
- partial-fragment degradation is understood;
- calibrated probabilities are treated as prioritization scores rather than tolerance probabilities;
- AMRFinderPlus or an equivalent interpretable annotation layer remains visible beside ML scores.
"""),
]

notebook = {
    "cells": cells,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
output = ROOT / "notebooks" / "16_benchmark_megares_similarity_robustness_and_calibration.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
