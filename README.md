# Genomic Language Models for Disinfectant-Tolerance Risk Screening

This repository is a reproducible feasibility-study scaffold for evaluating
whether genomic language models add useful signal beyond simpler sequence
models when predicting low-level disinfectant tolerance in bacteria.

The first benchmark is deliberately narrow:

> Can sequence-only models predict low-level benzalkonium-chloride (BC)
> tolerance in *Listeria monocytogenes* isolates while respecting bacterial
> population structure?

The initial endpoint is the published BC MIC threshold from Gmeiner et al.
(2025):

```text
BC MIC >= 1.25 mg/L -> tolerant
BC MIC <  1.25 mg/L -> sensitive
```

Use-level sanitizer survival, biofilm persistence, and pharmaceutical-water
risk prioritization are separate endpoints. They must not be treated as
interchangeable labels.

## Why Start Here?

The project has two goals:

1. Reproduce a credible sequence-only baseline before using a transformer.
2. Test whether frozen DNABERT2 embeddings justify more advanced modeling.

A later exploratory extension can profile public genomes from
pharmaceutical-water-relevant microorganisms such as *Burkholderia cepacia*
complex, *Ralstonia*, *Pseudomonas aeruginosa*, *Stenotrophomonas maltophilia*,
and *Methylobacterium*. That extension is a genomic risk-prioritization study,
not a validated sanitizer-survival predictor.

## Project Layout

```text
configs/                 versioned experiment settings
data/manifests/          accession, phenotype, lineage, and provenance tables
data/raw/                downloaded assemblies or reads, ignored by git
data/processed/          windows, splits, and embeddings, ignored by git
notebooks/               exploratory work only
reports/                 metrics and figures
scripts/                 runnable pipeline entry points
src/                     reusable Python utilities
tests/                   small correctness tests
```

## Setup

Create a Python environment and install the lightweight baseline:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Install PyTorch and Transformers only when the baseline workflow is ready:

```bash
python -m pip install -e ".[dnabert]"
```

## First Feasibility Run

### 1. Build The First Manifest In Notebook 01

Open:

```text
notebooks/01_build_listeria_feasibility_manifest.ipynb
```

The notebook downloads the public supplementary workbook, inspects its sheets,
normalizes the relevant columns, applies the published BC threshold, samples a
reproducible feasibility cohort, and writes:

```text
data/manifests/isolate_manifest.csv
```

The notebook leaves `group` blank when the selected supplementary table does not
contain a lineage-aware grouping column. In that case, add clonal-complex or
PopPUNK groups before running the split script.

### 2. Review The Manifest

Start from:

```text
data/manifests/isolate_manifest.example.csv
```

Create:

```text
data/manifests/isolate_manifest.csv
```

The manifest must include an isolate-level phenotype label and a lineage-aware
group, such as clonal complex or PopPUNK cluster. See
`data/manifests/README.md`.

### 3. Resolve Public Sequence Identifiers

Open:

```text
notebooks/02_resolve_ena_accessions.ipynb
```

The supplementary workbook mixes run, experiment, and sample accessions. This
notebook resolves those identifiers through ENA, writes an isolate-level
resolved manifest, and creates a run-level retrieval plan with estimated
download sizes.

### 4. Create Splits And A Smoke-Test Plan

Open:

```text
notebooks/03_make_lineage_splits_and_smoke_plan.ipynb
```

This notebook creates lineage-aware train/dev/test splits and selects a small,
inspectable FASTQ retrieval plan for the first assembly smoke test.

### 5. Download Public ENA Assemblies

Open:

```text
notebooks/04_download_public_ena_assemblies.ipynb
```

This notebook resolves ENA `SEQUENCE_ASSEMBLY` analysis records and downloads
compact public FASTA files. Prefer this route over raw FASTQ retrieval and local
assembly whenever public assemblies are available.

### 6. Download Public FASTQ Files

For a small public subset:

```bash
python scripts/download_fastq_plan.py \
  --plan data/manifests/ena_smoke_retrieval_plan.csv \
  --output-dir data/raw/fastq \
  --limit-isolates 3
```

The current modeling scripts consume assembled FASTA files. Assemblies can be
downloaded directly when available or generated from FASTQ files using a
versioned assembly workflow such as SPAdes plus QUAST.

### 7. Generate Leakage-Aware Splits

```bash
python scripts/make_splits.py \
  --manifest data/manifests/isolate_manifest.csv
```

All isolates in one lineage group are assigned to the same split. This avoids
inflated results caused by close relatives appearing in both training and test
data.

### 8. Generate Genomic Windows

```bash
python scripts/build_windows.py \
  --manifest data/manifests/isolate_manifest.csv
```

The default window size is `999` bases to mirror the original DNABERT notebook.

### 9. Train A Lightweight Baseline

```bash
python scripts/train_kmer_baseline.py
```

This trains a character-k-mer TF-IDF logistic-regression model at isolate
level. It is intentionally simple and should be treated as the first reference
point.

### 10. Extract Frozen DNABERT2 Embeddings

```bash
python scripts/embed_dnabert2.py
python scripts/train_embedding_baseline.py
```

The embedding script creates one vector per isolate by mean-pooling genomic
window embeddings. This is a feasibility baseline. A stronger later model can
replace mean pooling with attention-based multiple-instance learning.

### 11. Review Assemblies And Run The Smoke Baseline

Open:

```text
notebooks/05_review_assemblies_and_run_smoke_baseline.ipynb
```

After the 12-isolate assembly smoke test is complete, this notebook reviews
assembly availability, creates an assembly-backed manifest, generates genomic
windows, and runs the character-k-mer baseline as a technical validation.

### 12. Run The Full Public-Assembly Baseline

Open:

```text
notebooks/06_download_full_assembly_cohort_and_run_kmer_baseline.ipynb
```

### 26. Audit Pharma-Relevant Phenotype Datasets Before Adding Organisms

Open:

```text
notebooks/22_literature_phenotype_dataset_audit.ipynb
```

This notebook checks real literature references for organism-level biocide
phenotypes plus public genome access. It does not create labels from mechanism
papers. The current decision is that *Pseudomonas aeruginosa* DDAC data from
Pottier et al. (2023) is the only high-priority non-*Listeria* supervised
candidate; Bcc studies are biologically important but conditional until exact
strain-to-genome mappings are verified.

### 27. Verify The Pottier 2023 DDAC Candidate Before Modeling

Open:

```text
notebooks/23_verify_pottier_2023_ddac_manifest.ipynb
```

This notebook checks the best second-benchmark lead: Pottier et al. (2023)
*Pseudomonas aeruginosa* DDAC susceptibility with BioProject `PRJNA884650`.
The current reproducible decision is that NCBI Assembly exposes the expected
77 genomes, but the public figshare phenotype materials are figure/image files
rather than a row-level machine-readable DDAC table. Do not train a supervised
*P. aeruginosa* DDAC model unless row-level labels are recovered from an author
table, another public machine-readable supplement, or a separately audited
manual/OCR extraction.

### 28. Audit Combined Bcc Biocide Phenotype-Genome Evidence

Open:

```text
notebooks/24_bcc_combined_phenotype_genome_audit.ipynb
```

This notebook extracts row-level Bcc biocide phenotypes from public PMC tables
in Moore et al. (2009) and Kim/Wang et al. (2015), then queries NCBI Assembly
for exact or near-exact strain-to-genome matches. The current result is 378
row-level phenotype rows across 75 Bcc strains, but only 9 strains with
straightforward public assembly matches. That is strong biological support and
pharmaceutical relevance, but still too small for a supervised Bcc model.

This notebook downloads the remaining compact public FASTA assemblies with a
disk-space reserve, creates the full assembly-backed manifest, and runs the
first serious grouped-test character-k-mer baseline without writing a large
window-level CSV.

### 13. Run A Staged DNABERT2 Smoke Test

Open:

```text
notebooks/07_run_cpu_dnabert2_smoke_embeddings.ipynb
```

This notebook prepares a local PyTorch-attention fallback snapshot and validates
frozen DNABERT2 embedding extraction on four isolates and 16 windows per
isolate before attempting a larger transformer inference run. Use a GPU
environment for scaling if local inference is too slow.

### 14. Run Capped Full-Cohort DNABERT2 Embeddings

Open:

```text
notebooks/08_run_full_capped_dnabert2_embeddings.ipynb
```

This notebook builds a compact 16-window-per-isolate table, extracts frozen
DNABERT2 embeddings with Apple MPS, and evaluates the lineage-held-out
embedding baseline.

### 15. Evaluate Repeated Lineage-Aware Splits

Open:

```text
notebooks/09_evaluate_dnabert2_repeated_grouped_splits.ipynb
```

This notebook reuses the frozen embeddings, evaluates ten grouped split seeds,
selects thresholds using development isolates only, and reports test metric
distributions.

### 16. Compare Repeated K-mer Evaluation and Audit QAC Metadata

Open:

```text
notebooks/10_compare_repeated_kmer_and_audit_qac_metadata.ipynb
```

This notebook evaluates the capped assembly k-mer baseline across the same ten
lineage-aware splits, compares it with DNABERT2, and audits the supplementary
`no_known_QAC_gene` field. The field marks tolerant isolates without an observed
known QAC determinant; it is not a positive known-gene feature.

### 17. Complete the Bundled Feasibility Experiments

Run notebooks `11` through `13` in order. Notebook `11` screens four known QAC
determinants with BLAST+, evaluates the interpretable rule baseline, and reviews
unexplained tolerant isolates. Notebook `12` evaluates 18 DNABERT2
representation strategies from one reusable 64-window embedding extraction.
Notebook `13` builds the final comparison table and defines the stop/go rule for
an external-validation phase.

### 18. Audit the Harrand Use-Level Sanitizer-Survival Dataset

Open:

```text
notebooks/14_audit_harrand_use_level_sanitizer_survival_dataset.ipynb
```

This notebook downloads only small public audit files from Europe PMC and
Cornell eCommons. It inventories the 501-isolate cohort, records the compact
assembly-deposit size without downloading genomes, and checks whether public
isolate-level BC, PAA, and NaOCl log-reduction labels are available before any
external-validation modeling begins.

### 19. Audit MEGARes and Run a Homology-Aware Fragment Baseline

Open:

```text
notebooks/15_audit_megares_and_run_homology_aware_fragment_baseline.ipynb
```

This notebook begins the exploratory pharmaceutical-water extension. It
downloads the official MEGARes v3 nucleotide FASTA and annotations, constructs a
biocide/metal/multi-compound fragment manifest, and evaluates a k-mer baseline
with entire MEGARes gene groups held out. The output is fragment-level genomic
risk prioritization, not validated isolate-level sanitizer-survival prediction.

### 20. Benchmark MEGARes Similarity Robustness and Calibration

Open:

```text
notebooks/16_benchmark_megares_similarity_robustness_and_calibration.ipynb
```

This notebook builds BLAST-based nucleotide-similarity clusters, generates
centered partial-fragment and substitution-noise conditions, extracts reusable
DNABERT2 embeddings, and compares calibrated k-mer, random-forest, XGBoost, and
DNABERT2 strategies across repeated similarity-cluster-held-out splits.

### 21. Narrow the MEGARes Task and Audit BacMet Evidence

Open:

```text
notebooks/17_narrow_megares_biocide_task_and_audit_bacmet.ipynb
```

This notebook reuses the notebook 16 artifacts to evaluate a narrower binary
task: biocide-relevant resistance fragments versus metal-only resistance
fragments. It also audits the official BacMet experimentally confirmed protein
release and compound mapping table. BacMet is retained as an interpretable
annotation layer, not claimed as an independent DNABERT nucleotide-validation
set.

### 22. Plan a Disk-Capped Pharmaceutical-Water Genome Pilot

Open:

```text
notebooks/18_plan_disk_capped_pharma_water_genome_pilot.ipynb
```

This notebook queries the official NCBI Datasets API and writes a reviewable
metadata-only RefSeq pilot plan across four pharmaceutical-water-relevant taxa.
It estimates the FASTA footprint and keeps assembly download explicitly gated.

### 23. Download and Screen the Pharmaceutical-Water Pilot

Open:

```text
notebooks/19_download_and_screen_pharma_water_pilot.ipynb
```

This notebook retrieves only the reviewed RefSeq assemblies, extracts genomic
FASTA files without retaining duplicate ZIP packages, screens interpretable
MEGARes nucleotide matches, and scores reviewable 500 bp candidate regions with
the calibrated partial-fragment model. Scores are exploratory prioritization
signals, not isolate-level disinfectant-tolerance predictions.

### 24. Harden Pilot Evidence and Selectively Score Intact Regions

Open:

```text
notebooks/20_harden_pharma_water_evidence_and_selective_dnabert2.ipynb
```

This notebook downloads compact NCBI protein and GFF annotations, runs
AMRFinderPlus with the `--plus` catalog, screens experimentally supported
BacMet biocide proteins, merges nearby evidence onto the pilot regions, and
applies DNABERT2 selectively to a small set of intact candidate regions.

### 25. Build the Publication-Readiness Dossier

Open:

```text
notebooks/21_publication_readiness_and_final_dossier.ipynb
```

This notebook defines the final publication claim and stop rule. It creates the
benchmark table, candidate-region table, evidence-enrichment summary, and
claim-boundary text for a computational feasibility manuscript.

### 26. Audit Literature Phenotype Datasets

Open:

```text
notebooks/22_literature_phenotype_dataset_audit.ipynb
```

This notebook surveys public disinfectant-tolerance phenotype datasets that
could extend the supervised benchmark beyond *Listeria*. It records which
studies have row-level labels, public genome accessions, and suitable endpoints.

### 27. Verify Pottier 2023 *Pseudomonas* DDAC Materials

Open:

```text
notebooks/23_verify_pottier_2023_ddac_manifest.ipynb
```

This notebook checks whether the public *P. aeruginosa* DDAC materials can be
converted into a supervised phenotype-genome manifest. The current result is
blocked because the public phenotype files are aggregate figure/image materials,
not row-level strain labels.

### 28. Audit Bcc Phenotype-Genome Linkage

Open:

```text
notebooks/24_bcc_combined_phenotype_genome_audit.ipynb
```

This notebook combines row-level biocide phenotype records for Burkholderia
cepacia complex strains and searches for public assembly matches. The current
result is useful context, but insufficient exact public genome matches for a
second supervised model.

### 29. Audit *Acinetobacter* and *Klebsiella* Candidates

Open:

```text
notebooks/25_acinetobacter_klebsiella_phenotype_genome_audit.ipynb
```

This notebook audits *A. baumannii* and *K. pneumoniae* disinfectant-tolerance
papers. It finds one row-level *A. baumannii* MIC table, but no confirmed public
genome linkage, and a *K. pneumoniae* study with public genomes but aggregate
biocide MICs. The stop-rule decision is not to start a second supervised
organism model without author-provided or recovered row-level phenotype-genome
data.

### 30. Audit Multi-Organism ML-Ready Biocide Datasets

Open:

```text
notebooks/26_audit_multiorganism_ml_ready_biocide_datasets.ipynb
```

This notebook ranks additional candidate organisms and datasets by supervised
ML readiness. The immediate audit batch is *Proteus mirabilis*, *Enterococcus
faecium*, *Enterococcus faecalis*, and a mixed hospital-environment
chlorhexidine preprint. Smaller *Pseudomonas*, food-drain, and *Staphylococcus*
datasets are treated as external-case candidates unless row-level
phenotype-genome linkage is stronger than expected.

### 31. Verify Immediate Multi-Organism Candidates

Open:

```text
notebooks/27_verify_immediate_multiorganism_candidates.ipynb
```

This notebook checks public metadata for the immediate audit batch. It finds
*Proteus mirabilis* as the only candidate ready for row-level manifest extraction
from automated metadata, because a machine-readable Table S1 supplement and
public SRA signal are exposed. *E. faecium* and the hospital-environment
preprint remain near misses, and *E. faecalis* requires a manual retry because
the automated PMC fetch was blocked by a provider challenge page.

### 32. Extract the Proteus 2025 Candidate Manifest

Open:

```text
notebooks/28_extract_proteus_2025_candidate_manifest.ipynb
```

This notebook attempts row-level manifest extraction for the *Proteus mirabilis*
chlorhexidine/cationic-biocide candidate. The manually downloaded Table S1 XLSX
provides 78 isolate/accession rows, but no row-level MIC phenotype columns were
found. Treat Proteus as an accession/genomics near miss unless row-level MIC
labels are recovered from authors or another source.

## Research Progression

Proceed only when the previous stage is reproducible:

1. Process a 100-200 isolate subset.
2. Reproduce the k-mer baseline.
3. Compare frozen DNABERT2 embeddings.
4. Scale to the full labeled dataset.
5. Add phylogeny-aware nested cross-validation.
6. Compare with pan-genome features and known QAC-gene rules.
7. Evaluate independent datasets.
8. Test use-level sanitizer-survival outcomes as a separate endpoint.
9. Add the exploratory pharmaceutical-water extension.

## Scientific Guardrails

- Do not synthesize MIC or sanitizer-survival labels.
- Use synthetic sequence noise only for robustness experiments.
- Keep low-level MIC tolerance separate from use-level sanitizer survival.
- Report grouped and external validation results, not random-split results
  alone.
- Treat influential model regions as hypotheses, not causal mechanisms.
- Frame the pharmaceutical-water extension as research-grade risk
  prioritization.

## Source Workflows

- [Gmeiner et al. 2025](https://www.nature.com/articles/s41598-025-94321-6):
  public *Listeria* WGS and QAC phenotype benchmark with phylogeny-aware ML.
- [LmonoDisinfectML](https://github.com/agmei/LmonoDisinfectML):
  published conventional-model scripts.
- [Harrand et al. 2026](https://journals.asm.org/doi/10.1128/aem.01060-25):
  separate public use-level sanitizer-survival dataset.
- [DNABERT2](https://github.com/MAGICS-LAB/DNABERT_2):
  official PyTorch and Hugging Face implementation.
- [ENA Browser Tools](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-tools.html):
  public sequence retrieval guidance.
- [Panaroo](https://gthlab.au/panaroo/) and
  [PopPUNK](https://poppunk.readthedocs.io/):
  pan-genome analysis and lineage clustering.
- [MEGARes](https://www.meglab.org/megares/) and
  [AMR++](https://github.com/Microbial-Ecology-Group/AMRplusplus):
  interpretable resistome profiling for the pharmaceutical-water extension.
- [BacMet](http://bacmet.biomedicine.gu.se/):
  experimentally confirmed biocide- and metal-resistance annotations.
