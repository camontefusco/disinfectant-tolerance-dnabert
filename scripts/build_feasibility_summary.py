#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reports", default="reports")
    parser.add_argument("--output", default="reports/final_feasibility_summary.csv")
    parser.add_argument("--markdown", default="reports/final_feasibility_summary.md")
    args = parser.parse_args()

    reports = Path(args.reports)
    rows = []
    for filename, model in [
        ("kmer_repeated_grouped_evaluation.json", "Capped k-mer TF-IDF"),
        ("dnabert2_repeated_grouped_evaluation.json", "DNABERT2 16-window mean + logistic regression"),
    ]:
        report = json.loads((reports / filename).read_text())
        rows.append(
            {
                "model": model,
                "balanced_accuracy_mean": report["summary"]["balanced_accuracy"]["mean"],
                "balanced_accuracy_std": report["summary"]["balanced_accuracy"]["std"],
                "roc_auc_mean": report["summary"]["roc_auc"]["mean"],
                "roc_auc_std": report["summary"]["roc_auc"]["std"],
            }
        )
    gene_path = reports / "qac_gene_rule_baseline.json"
    if gene_path.exists():
        report = json.loads(gene_path.read_text())
        rows.append(
            {
                "model": "Known-QAC determinant BLAST rule",
                "balanced_accuracy_mean": report["balanced_accuracy"],
                "balanced_accuracy_std": None,
                "roc_auc_mean": report["roc_auc"],
                "roc_auc_std": None,
            }
        )
    strategy_path = reports / "dnabert2_strategy_evaluation.json"
    if strategy_path.exists():
        report = json.loads(strategy_path.read_text())
        best = max(report["summary"], key=lambda row: row["roc_auc_mean"])
        rows.append(
            {
                "model": f"Best sensitivity strategy: {best['cap']} windows, {best['aggregation']}, {best['classifier']}",
                "balanced_accuracy_mean": best["balanced_accuracy_mean"],
                "balanced_accuracy_std": best["balanced_accuracy_std"],
                "roc_auc_mean": best["roc_auc_mean"],
                "roc_auc_std": best["roc_auc_std"],
            }
        )

    frame = pd.DataFrame(rows).sort_values("roc_auc_mean", ascending=False)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    columns = list(frame.columns)
    markdown_table = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.round(3).fillna("").iterrows():
        markdown_table.append("| " + " | ".join(str(row[column]) for column in columns) + " |")
    markdown = [
        "# Disinfectant-tolerance DNABERT2 feasibility summary",
        "",
        *markdown_table,
        "",
        "## Interpretation",
        "",
        "- Treat these results as a feasibility package, not a deployable model.",
        "- The endpoint is low-level BC MIC tolerance, not use-level sanitizer survival.",
        "- External validation remains the next requirement for a stronger manuscript claim.",
    ]
    Path(args.markdown).write_text("\n".join(markdown) + "\n")
    print(frame.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
