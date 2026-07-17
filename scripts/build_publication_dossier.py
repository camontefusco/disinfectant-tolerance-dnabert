#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def load_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def evidence_flag(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(False, index=frame.index)
    return frame[column].notna() & frame[column].astype(str).ne("")


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else 0.0


def markdown_table(frame: pd.DataFrame) -> str:
    display = frame.copy().astype(object).where(pd.notna(frame), "")
    for column in display.columns:
        if pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].round(3)
    rows = ["| " + " | ".join(map(str, display.columns)) + " |"]
    rows.append("| " + " | ".join(["---"] * len(display.columns)) + " |")
    for row in display.astype(str).itertuples(index=False):
        rows.append("| " + " | ".join(row) + " |")
    return "\n".join(rows)


def source_enrichment(frame: pd.DataFrame, flag: pd.Series) -> dict[str, float | str]:
    hit = frame["source"].eq("megares_hit")
    background = frame["source"].eq("background_sample")
    hit_rate = safe_ratio(float((flag & hit).sum()), float(hit.sum()))
    background_rate = safe_ratio(float((flag & background).sum()), float(background.sum()))
    if background_rate == 0 and hit_rate > 0:
        fold_enrichment: float | str = "not estimable; background=0"
    else:
        fold_enrichment = safe_ratio(hit_rate, background_rate)
    return {
        "hit_rate": hit_rate,
        "background_rate": background_rate,
        "rate_difference": hit_rate - background_rate,
        "fold_enrichment": fold_enrichment,
    }


def benchmark_rows(binary_report: dict) -> pd.DataFrame:
    summary = pd.DataFrame(binary_report["summary"])
    keep = summary[
        summary["model"].isin(["dnabert2_xgboost", "kmer_random_forest"])
        & summary["condition"].isin(["full", "noise_1%", "partial_500bp"])
    ].copy()
    return keep[
        [
            "model", "condition", "calibrated_balanced_accuracy_mean",
            "calibrated_balanced_accuracy_std", "calibrated_roc_auc_mean",
            "calibrated_roc_auc_std", "calibrated_average_precision_mean",
        ]
    ].sort_values(["condition", "calibrated_balanced_accuracy_mean"], ascending=[True, False])


def candidate_table(evidence: pd.DataFrame, limit: int = 40) -> pd.DataFrame:
    table = evidence.copy()
    table["has_amrfinderplus"] = evidence_flag(table, "amrfinder_element_symbol")
    table["has_bacmet_biocide"] = evidence_flag(table, "bacmet_bacmet_id")
    table["has_bacmet_qac"] = table.get("bacmet_has_qac_annotation", pd.Series("", index=table.index)).astype(str).str.lower().eq("true")
    table["evidence_channels"] = (
        table["has_amrfinderplus"].astype(int)
        + table["has_bacmet_biocide"].astype(int)
        + table["has_bacmet_qac"].astype(int)
        + table["best_megares_class"].astype(str).str.len().gt(0).astype(int)
        + table["selective_dnabert2_biocide_relevance_score"].notna().astype(int)
    )
    columns = [
        "region_id", "taxon_query", "assembly_accession", "best_megares_class",
        "exploratory_biocide_relevance_score", "selective_dnabert2_biocide_relevance_score",
        "evidence_channels", "amrfinder_element_symbol", "amrfinder_class",
        "bacmet_gene_name", "bacmet_has_qac_annotation",
    ]
    return table.sort_values(
        ["evidence_channels", "selective_dnabert2_biocide_relevance_score", "exploratory_biocide_relevance_score"],
        ascending=[False, False, False],
    )[columns].head(limit)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", default="data/manifests/pharma_water_pilot_evidence_matrix_with_dnabert2.csv")
    parser.add_argument("--binary-report", default="reports/megares_binary_biocide_relevance_benchmark.json")
    parser.add_argument("--listeria-summary", default="reports/final_feasibility_summary.csv")
    parser.add_argument("--amrfinder-report", default="reports/pharma_water_pilot_amrfinderplus.json")
    parser.add_argument("--evidence-report", default="reports/pharma_water_pilot_evidence_matrix.json")
    parser.add_argument("--top-candidates", default="reports/publication_top_candidate_regions.csv")
    parser.add_argument("--benchmark-table", default="reports/publication_benchmark_table.csv")
    parser.add_argument("--output", default="reports/publication_readiness_dossier.md")
    args = parser.parse_args()

    evidence = pd.read_csv(args.evidence)
    binary = load_json(args.binary_report)
    listeria = pd.read_csv(args.listeria_summary)
    amrfinder = load_json(args.amrfinder_report)
    evidence_report = load_json(args.evidence_report)
    top = candidate_table(evidence)
    benchmark = benchmark_rows(binary)
    Path(args.top_candidates).parent.mkdir(parents=True, exist_ok=True)
    top.to_csv(args.top_candidates, index=False)
    benchmark.to_csv(args.benchmark_table, index=False)

    hit = evidence["source"].eq("megares_hit")
    background = evidence["source"].eq("background_sample")
    partial_summary = evidence.groupby("source")["exploratory_biocide_relevance_score"].agg(["count", "median", "mean", "max"]).round(3)
    enrichments = {
        "AMRFinderPlus nearby": source_enrichment(evidence, evidence_flag(evidence, "amrfinder_element_symbol")),
        "BacMet biocide nearby": source_enrichment(evidence, evidence_flag(evidence, "bacmet_bacmet_id")),
        "BacMet QAC nearby": source_enrichment(evidence, evidence.get("bacmet_has_qac_annotation", pd.Series("", index=evidence.index)).astype(str).str.lower().eq("true")),
    }
    selective = evidence[evidence["selective_dnabert2_biocide_relevance_score"].notna()]
    correlation = float(selective["exploratory_biocide_relevance_score"].corr(selective["selective_dnabert2_biocide_relevance_score"]))

    lines = [
        "# Publication-readiness dossier",
        "",
        "## Recommended final goal",
        "",
        "Submit this as a bounded computational feasibility/methods manuscript:",
        "",
        "> A leakage-aware DNABERT2 and annotation-anchored workflow for prioritizing disinfectant-resistance-associated genomic regions, validated first on published *Listeria monocytogenes* benzalkonium-chloride MIC labels and then applied cautiously to public pharmaceutical-water-relevant opportunists.",
        "",
        "Do **not** claim validated pharmaceutical-water disinfectant tolerance prediction. The pharmaceutical-water section is exploratory genomic risk prioritization.",
        "",
        "## Stop rule",
        "",
        "Stop adding new organisms, new model families, and new public datasets. The remaining work is manuscript assembly: figures, tables, limitations, and reproducibility cleanup.",
        "",
        "## Core result blocks",
        "",
        "1. Published-label benchmark: grouped *Listeria* BC-MIC feasibility with DNABERT2, k-mer baselines, and known-QAC determinant context.",
        "2. Curated-reference benchmark: MEGARes biocide-relevant versus metal-only fragments under nucleotide-similarity-held splits, partial-fragment tests, noise tests, and calibration.",
        "3. Pharmaceutical-water pilot: twenty RefSeq assemblies from Bcc, *Ralstonia pickettii*, *Pseudomonas aeruginosa*, and *Stenotrophomonas maltophilia* scored with MEGARes, AMRFinderPlus, BacMet, k-mer, and selective DNABERT2 evidence.",
        "",
        "## MEGARes benchmark snapshot",
        "",
        markdown_table(benchmark),
        "",
        "## Published-label Listeria benchmark snapshot",
        "",
        markdown_table(listeria),
        "",
        "## Pharmaceutical-water evidence snapshot",
        "",
        f"- Scored regions: {len(evidence):,}",
        f"- MEGARes-hit regions: {int(hit.sum()):,}",
        f"- Background regions: {int(background.sum()):,}",
        f"- Regions with nearby AMRFinderPlus evidence: {evidence_report['regions_with_nearby_amrfinderplus']:,}",
        f"- Regions with nearby BacMet biocide-protein evidence: {evidence_report['regions_with_nearby_bacmet_biocide_protein']:,}",
        f"- Regions with nearby BacMet QAC-protein evidence: {evidence_report['regions_with_nearby_bacmet_qac_protein']:,}",
        f"- AMRFinderPlus database version: {amrfinder['database_version']}",
        f"- Selective partial-versus-intact score correlation: {correlation:.3f}",
        "",
        "### Region score distribution",
        "",
        markdown_table(partial_summary.reset_index()),
        "",
        "### Evidence enrichment over deterministic background",
        "",
        markdown_table(pd.DataFrame(enrichments).T.reset_index().rename(columns={"index": "evidence"})),
        "",
        "## Strongest candidate-region classes",
        "",
        markdown_table(top.head(20)),
        "",
        "## Publication claim boundary",
        "",
        "Defensible: the workflow prioritizes genomic regions enriched for biocide/stress-response evidence and shows concordance between simple sequence models, DNABERT2 embeddings, and curated annotations.",
        "",
        "Not defensible without wet-lab data: predicting sanitizer survival, pharmaceutical-water persistence, MIC values for the pilot organisms, or causal resistance mechanisms.",
        "",
        "## Remaining manuscript-only tasks",
        "",
        "1. Convert this dossier into four main figures and two supplemental tables.",
        "2. Write methods around leakage-aware splits, calibration, and database-version pinning.",
        "3. Keep the pharmaceutical-water result as exploratory prioritization with explicit no-phenotype limitation.",
        "4. Add a short future-validation paragraph: test top-ranked Bcc, *P. aeruginosa*, and *S. maltophilia* regions against measured QAC/PAA/NaOCl tolerance if collaborators become available.",
    ]
    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.top_candidates}")
    print(f"Wrote {args.benchmark_table}")


if __name__ == "__main__":
    main()
