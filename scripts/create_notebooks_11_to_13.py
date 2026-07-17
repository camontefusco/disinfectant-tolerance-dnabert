#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"


def markdown(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)}


def code(text: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": text.splitlines(keepends=True)}


def write(name: str, cells: list[dict]) -> None:
    notebook = {
        "cells": cells,
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    (NOTEBOOKS / name).write_text(json.dumps(notebook, indent=1) + "\n")


SETUP = """from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

PROJECT_ROOT = Path.cwd()
if PROJECT_ROOT.name == "notebooks":
    PROJECT_ROOT = PROJECT_ROOT.parent
print("Project root:", PROJECT_ROOT)
"""


write(
    "11_screen_known_qac_determinants_and_explain_errors.ipynb",
    [
        markdown("""# 11 - Screen known QAC determinants and explain errors

This notebook bundles the interpretable baseline: retrieve a curated public reference panel, screen the 197 existing assemblies with BLAST+, evaluate a known-determinant rule, and inspect tolerant isolates not explained by that rule. It does not download additional genomes.
"""),
        markdown("""## 1. Configuration

The reference panel contains `bcrABC`, `qacH`, `emrE`, and `emrC`. The default rule requires at least 80% nucleotide identity and 80% reference coverage. This is a compact hypothesis-oriented screen, not a replacement for a full resistome workflow.
"""),
        code(SETUP + """
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "isolate_manifest_public_assemblies.csv"
REFERENCES = PROJECT_ROOT / "data" / "references" / "qac_determinants.fasta"
SCREEN = PROJECT_ROOT / "data" / "processed" / "qac_determinant_screen.csv"
REPORT = PROJECT_ROOT / "reports" / "qac_gene_rule_baseline.json"
RUN_REFERENCE_FETCH = False
RUN_QAC_SCREEN = False
"""),
        markdown("""## 2. Retrieve the small curated reference panel

Set `RUN_REFERENCE_FETCH = True` for the first run. This downloads only four short public NCBI sequence records.
"""),
        code("""command = [sys.executable, str(PROJECT_ROOT / "scripts" / "fetch_qac_references.py"), "--output", str(REFERENCES)]
print(" ".join(command))
if RUN_REFERENCE_FETCH:
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
print("Reference FASTA exists:", REFERENCES.exists())
"""),
        markdown("""## 3. Run BLAST screening and evaluate the rule baseline

Set `RUN_QAC_SCREEN = True`. The screen reuses the compact FASTA assemblies already on disk.
"""),
        code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "screen_qac_determinants.py"),
    "--manifest", str(MANIFEST), "--references", str(REFERENCES),
    "--output", str(SCREEN), "--report", str(REPORT),
]
print(" ".join(command))
if RUN_QAC_SCREEN:
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
report = json.loads(REPORT.read_text())
display(pd.Series(report))
"""),
        markdown("""## 4. Review determinant frequencies and unexplained tolerant isolates

The most interesting subset is tolerant isolates without one of the four screened determinants. Those isolates test whether ML contributes information beyond a simple known-gene rule.
"""),
        code("""hits = pd.read_csv(SCREEN)
manifest = pd.read_csv(MANIFEST)
presence = hits.pivot(index="isolate_id", columns="gene", values="present").reset_index()
review = manifest.merge(presence, on="isolate_id", validate="one_to_one")
gene_columns = ["bcrABC", "qacH", "emrE", "emrC"]
review["known_qac_gene"] = review[gene_columns].any(axis=1)
print("Gene-positive counts:")
display(review[gene_columns].sum().sort_values(ascending=False).to_frame("isolates"))
print("Tolerance by known-QAC rule:")
display(pd.crosstab(review["known_qac_gene"], review["label"], margins=True))
print("Tolerant isolates not explained by the four-reference panel:")
display(review[(review["label"] == 1) & ~review["known_qac_gene"]][["isolate_id", "group", "mic_bc_mg_l", "source_metadata"]])
"""),
        markdown("""## 5. Decision

Notebook 12 is worth running if DNABERT2 should be tested as a complement to the known-gene rule. Keep the claim narrow: the rule screen is interpretable but intentionally compact, while a manuscript workflow should later add a full annotated resistome pipeline.
"""),
    ],
)

write(
    "12_run_dnabert2_representation_sensitivity_grid.ipynb",
    [
        markdown("""# 12 - Run the DNABERT2 representation sensitivity grid

This notebook bundles the model-strategy experiments. It creates one 64-window table, extracts window-level DNABERT2 embeddings once, and evaluates 18 strategies: 3 window caps x 2 aggregations x 3 classifiers. Each strategy is evaluated across the same 10 lineage-aware partitions.
"""),
        markdown("""## 1. Configuration and storage estimate

The largest new artifact is the compressed window-level embedding archive. Expect roughly 40 MB, not gigabytes. DNABERT2 inference is the slow stage and may take approximately 15-25 minutes on Apple MPS.
"""),
        code(SETUP + """
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "isolate_manifest_public_assemblies.csv"
WINDOWS = PROJECT_ROOT / "data" / "processed" / "dnabert2_full_64_sampled_windows.csv"
MEAN_EMBEDDINGS = PROJECT_ROOT / "data" / "processed" / "dnabert2_full_64_mean_embeddings.npz"
WINDOW_EMBEDDINGS = PROJECT_ROOT / "data" / "processed" / "dnabert2_full_64_window_embeddings.npz"
REPORT = PROJECT_ROOT / "reports" / "dnabert2_strategy_evaluation.json"
LOCAL_MODEL = PROJECT_ROOT / "data" / "external" / "dnabert2_pytorch_fallback"
BUILD_64_WINDOWS = True
EMBED_64_WINDOWS = True
RUN_STRATEGY_GRID = True
"""),
        markdown("""## 2. Build evenly spaced 64-window input

This remains compact and reuses the existing assemblies. Completed artifacts are reused automatically, so interrupted runs can resume with **Run All**.
"""),
        code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "build_sampled_windows.py"),
    "--manifest", str(MANIFEST), "--output", str(WINDOWS),
    "--maximum-windows-per-isolate", "64",
]
print(" ".join(command))
if BUILD_64_WINDOWS and not WINDOWS.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
if WINDOWS.exists():
    print("Window table MB:", round(WINDOWS.stat().st_size / 1_000_000, 2))
else:
    print("Window table is missing. Set BUILD_64_WINDOWS = True and rerun this cell.")
"""),
        markdown("""## 3. Extract window-level DNABERT2 embeddings once

The local PyTorch-attention fallback prepared earlier is reused. Existing embeddings are not recomputed.
"""),
        code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "embed_dnabert2.py"),
    "--windows", str(WINDOWS), "--model", str(LOCAL_MODEL),
    "--maximum-windows-per-isolate", "64", "--output", str(MEAN_EMBEDDINGS),
    "--window-output", str(WINDOW_EMBEDDINGS),
]
print(" ".join(command))
if EMBED_64_WINDOWS and not WINDOW_EMBEDDINGS.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
if WINDOW_EMBEDDINGS.exists():
    print("Window embedding archive MB:", round(WINDOW_EMBEDDINGS.stat().st_size / 1_000_000, 2))
else:
    print("Window embeddings are missing. Build windows first, then set EMBED_64_WINDOWS = True.")
"""),
        markdown("""## 4. Evaluate 18 representation and classifier strategies

This stage should be much faster than transformer inference. Existing reports are reused.
"""),
        code("""command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "evaluate_dnabert2_strategies.py"),
    "--window-embeddings", str(WINDOW_EMBEDDINGS), "--manifest", str(MANIFEST),
    "--output", str(REPORT),
]
print(" ".join(command))
if RUN_STRATEGY_GRID and WINDOW_EMBEDDINGS.exists() and not REPORT.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
if REPORT.exists():
    report = json.loads(REPORT.read_text())
    summary = pd.DataFrame(report["summary"]).sort_values("roc_auc_mean", ascending=False)
    display(summary.round(3))
else:
    print("Strategy report is missing. Complete the window-embedding stage above, then rerun this cell.")
"""),
        markdown("""## 5. Decision

Prefer a strategy only if its improvement is reasonably consistent across lineage-aware folds. Small differences should be treated as uncertainty, not as a reason to add model complexity.
"""),
    ],
)

write(
    "13_build_final_feasibility_decision_package.ipynb",
    [
        markdown("""# 13 - Build the final feasibility decision package

This notebook gathers the completed experiments into one concise evidence package. It is the stopping point for the internal proof-of-concept phase.
"""),
        markdown("""## 1. Build the final comparison table

This includes the conventional k-mer baseline, the original frozen-DNABERT2 baseline, the interpretable QAC rule if available, and the strongest sensitivity-grid strategy if available.
"""),
        code(SETUP + """
SUMMARY_CSV = PROJECT_ROOT / "reports" / "final_feasibility_summary.csv"
SUMMARY_MD = PROJECT_ROOT / "reports" / "final_feasibility_summary.md"
RUN_SUMMARY = True
command = [
    sys.executable, str(PROJECT_ROOT / "scripts" / "build_feasibility_summary.py"),
    "--output", str(SUMMARY_CSV), "--markdown", str(SUMMARY_MD),
]
print(" ".join(command))
if RUN_SUMMARY:
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
summary = pd.read_csv(SUMMARY_CSV)
display(summary.round(3))
"""),
        markdown("""## 2. Check whether the feasibility target has been reached

The immediate target is:

> Determine whether a genome-language-model representation identifies low-level benzalkonium-chloride MIC tolerance signal beyond lineage-aware k-mer and known-QAC-determinant baselines.

Use this checklist:

| Requirement | Evidence | Interpretation |
|---|---|---|
| Public measured phenotype labels | 197 assembly-linked isolates from the published cohort | Satisfied for feasibility |
| Lineage-aware evaluation | Repeated grouped train/dev/test partitions | Satisfied for feasibility |
| Conventional baseline | Repeated capped k-mer TF-IDF evaluation | Satisfied |
| Interpretable baseline | BLAST+ screen for four known QAC determinants | Satisfied as a compact baseline; expand for manuscript |
| Representation sensitivity | 18 DNABERT2 strategies | Satisfied after notebook 12 |
| External validation | Independent compatible MIC-labeled cohort | Still missing |
| Facility-level sanitizer-survival claim | Separate measured endpoint | Out of scope for this phase |
"""),
        markdown("""## 3. Stop/go rule

Stop the internal proof-of-concept after this notebook. Proceed to an external-validation phase only if DNABERT2 remains directionally stronger than the conventional and rule baselines across grouped folds, and if the unexplained tolerant subset remains scientifically interesting.

For a portfolio or conference-ready project, notebook 13 is enough. For a manuscript-oriented claim, the next required work is external validation and a broader resistome annotation workflow, not endless tuning on this cohort.
"""),
        code("""print(SUMMARY_MD.read_text())
"""),
    ],
)

print("Wrote notebooks 11, 12, and 13")
