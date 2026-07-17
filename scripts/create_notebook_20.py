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
    markdown("""# 20 - Harden pharmaceutical-water evidence and selectively score intact regions

This notebook adds an interpretable biological evidence layer to the pharmaceutical-water pilot. It retrieves compact NCBI protein and GFF annotations, runs AMRFinderPlus with the `--plus` catalog, screens experimentally supported BacMet biocide proteins, links nearby calls to candidate regions, and selectively applies DNABERT2 to intact high-priority regions.

The output remains **exploratory genomic risk prioritization**, not validated disinfectant-tolerance prediction.
"""),
    markdown("""## 1. Configuration

Use **Run All**. Every expensive stage is cached. The project-local AMRFinderPlus toolchain and its pinned NCBI database are prepared outside the notebook so their disk footprint can be reviewed explicitly.
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
PYTHON = sys.executable
ANNOTATIONS = PROJECT_ROOT / "data" / "manifests" / "pharma_water_pilot_annotation_files.csv"
AMRFINDER_HITS = PROJECT_ROOT / "data" / "manifests" / "pharma_water_pilot_amrfinderplus_hits.csv"
BACMET_HITS = PROJECT_ROOT / "data" / "manifests" / "pharma_water_pilot_bacmet_biocide_hits.csv"
EVIDENCE = PROJECT_ROOT / "data" / "manifests" / "pharma_water_pilot_evidence_matrix.csv"
CANDIDATES = PROJECT_ROOT / "data" / "processed" / "pharma_water_selective_dnabert2_candidates.csv"
CANDIDATE_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "pharma_water_selective_dnabert2_candidates.csv"
CANDIDATE_EMBEDDINGS = PROJECT_ROOT / "data" / "processed" / "pharma_water_selective_dnabert2_embeddings.npz"
FINAL = PROJECT_ROOT / "data" / "manifests" / "pharma_water_pilot_evidence_matrix_with_dnabert2.csv"
LOCAL_MODEL = PROJECT_ROOT / "data" / "external" / "dnabert2_pytorch_fallback"
"""),
    markdown("""## 2. Retrieve compact NCBI protein and GFF annotations

The same twenty RefSeq accessions are used. ZIP packages are held in memory only long enough to extract protein FASTA and GFF files.
"""),
    code("""command = [PYTHON, str(PROJECT_ROOT / "scripts" / "download_pharma_water_ncbi_annotations.py"), "--execute"]
print(" ".join(command))
if not ANNOTATIONS.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
annotations = pd.read_csv(ANNOTATIONS)
print("Annotation packages:", len(annotations))
print("Protein and GFF MB:", round((annotations["protein_bytes"].sum() + annotations["gff_bytes"].sum()) / 1_000_000, 2))
"""),
    markdown("""## 3. Run AMRFinderPlus with the stress-response catalog

AMRFinderPlus receives PGAP protein inputs plus `--plus`. The runner restores genomic coordinates from the downloaded PGAP GFF files and caches per-assembly TSV files. This cross-organism pilot uses a comparable generic annotation layer; organism-specific mutation calling is not claimed.
"""),
    code("""command = [PYTHON, str(PROJECT_ROOT / "scripts" / "run_amrfinderplus_pilot.py")]
print(" ".join(command))
if not AMRFINDER_HITS.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
amr = pd.read_csv(AMRFINDER_HITS)
display(pd.crosstab(amr["Type"], amr["Scope"], margins=True))
display(amr[["assembly_accession", "taxon_query", "Element symbol", "Element name", "Scope", "Class", "Subclass"]].head(30))
"""),
    markdown("""## 4. Screen experimentally supported BacMet biocide proteins

BacMet is screened as a translated-protein annotation layer. Protein matches support interpretation; they are not an independent nucleotide DNABERT validation dataset.
"""),
    code("""command = [PYTHON, str(PROJECT_ROOT / "scripts" / "screen_bacmet_protein_pilot.py")]
print(" ".join(command))
if not BACMET_HITS.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
bacmet = pd.read_csv(BACMET_HITS)
print("BacMet biocide-protein hits:", len(bacmet))
display(bacmet[["assembly_accession", "taxon_query", "BacMet_ID", "Gene_name", "identity", "has_qac_annotation"]].head(30))
"""),
    markdown("""## 5. Merge nearby interpretable evidence

The evidence matrix links the nearest AMRFinderPlus and BacMet calls within `2,000 bp` of each scored region.
"""),
    code("""command = [PYTHON, str(PROJECT_ROOT / "scripts" / "build_pharma_water_evidence_matrix.py")]
print(" ".join(command))
if not EVIDENCE.exists():
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
evidence = pd.read_csv(EVIDENCE)
display(evidence.sort_values("exploratory_biocide_relevance_score", ascending=False).head(30))
"""),
    markdown("""## 6. Selectively score intact candidate regions with DNABERT2

Only the top `25` MEGARes-centered regions per organism group are embedded. The final score averages five similarity-held calibrated DNABERT2 plus XGBoost models from the notebook 17 strategy.
"""),
    code("""prepare = [PYTHON, str(PROJECT_ROOT / "scripts" / "prepare_selective_dnabert2_candidates.py")]
if not CANDIDATES.exists():
    subprocess.run(prepare, cwd=PROJECT_ROOT, check=True)
embed = [
    PYTHON, str(PROJECT_ROOT / "scripts" / "embed_dnabert2.py"),
    "--windows", str(CANDIDATES), "--model", str(LOCAL_MODEL),
    "--maximum-windows-per-isolate", "1", "--output", str(CANDIDATE_EMBEDDINGS),
]
print(" ".join(embed))
if not CANDIDATE_EMBEDDINGS.exists():
    subprocess.run(embed, cwd=PROJECT_ROOT, check=True)
score = [PYTHON, str(PROJECT_ROOT / "scripts" / "score_selective_dnabert2_candidates.py")]
if not FINAL.exists():
    subprocess.run(score, cwd=PROJECT_ROOT, check=True)
final = pd.read_csv(FINAL)
display(final[final["selective_dnabert2_biocide_relevance_score"].notna()].sort_values("selective_dnabert2_biocide_relevance_score", ascending=False).head(30))
"""),
    markdown("""## 7. Interpretation gate

Prioritize regions where multiple evidence channels agree:

1. interpretable MEGARes class;
2. nearby AMRFinderPlus `--plus` call;
3. nearby experimentally supported BacMet biocide or QAC protein;
4. elevated partial-fragment score;
5. elevated selective intact-region DNABERT2 score.

These channels identify review-worthy genomic regions. They do not establish disinfectant tolerance, pharmaceutical-water persistence, or causal mechanisms.
"""),
]
notebook = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.10"}}, "nbformat": 4, "nbformat_minor": 5}
output = ROOT / "notebooks" / "20_harden_pharma_water_evidence_and_selective_dnabert2.ipynb"
output.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {output}")
