#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, classification_report, roc_auc_score
from sklearn.pipeline import Pipeline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.io import read_fasta


def assembly_text(
    path: str | Path,
    window_size: int,
    maximum_windows: int,
) -> str:
    sequence = "N".join(sequence for _, sequence in read_fasta(path))
    if len(sequence) <= window_size:
        return sequence
    maximum_start = len(sequence) - window_size
    window_count = min(maximum_windows, maximum_start + 1)
    if window_count == 1:
        starts = [0]
    else:
        starts = [
            round(index * maximum_start / (window_count - 1))
            for index in range(window_count)
        ]
    return "N".join(sequence[start : start + window_size] for start in starts)


def attach_splits(manifest: pd.DataFrame, splits: pd.DataFrame) -> pd.DataFrame:
    if "split" in manifest.columns:
        data = manifest.copy()
    else:
        data = manifest.merge(
            splits[["isolate_id", "split"]],
            on="isolate_id",
            how="left",
            validate="one_to_one",
        )
    if data["split"].isna().any():
        raise ValueError("Every manifest isolate must have a split assignment")
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--splits", default="data/processed/isolate_splits.csv")
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--output", default="reports/kmer_assembly_baseline.json")
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    manifest = pd.read_csv(args.manifest, dtype={"label": int})
    splits = pd.read_csv(args.splits)
    data = attach_splits(manifest, splits)
    data["sequence"] = data["assembly_path"].map(
        lambda path: assembly_text(
            path,
            window_size=config["window_size"],
            maximum_windows=config["kmer_max_windows_per_isolate"],
        )
    )

    train = data[data["split"] == "train"]
    test = data[data["split"] == "test"]
    if train.empty or test.empty:
        raise ValueError("Train and test splits must both contain isolates")

    low, high = config["kmer_range"]
    model = Pipeline(
        [
            (
                "tfidf",
                TfidfVectorizer(
                    analyzer="char",
                    ngram_range=(low, high),
                    lowercase=False,
                    max_features=100_000,
                ),
            ),
            (
                "classifier",
                LogisticRegression(
                    class_weight="balanced",
                    max_iter=2_000,
                    random_state=config["random_seed"],
                ),
            ),
        ]
    )
    model.fit(train["sequence"], train["label"])
    probabilities = model.predict_proba(test["sequence"])[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    result = {
        "model": "assembly_character_kmer_tfidf_logistic_regression",
        "balanced_accuracy": balanced_accuracy_score(test["label"], predictions),
        "roc_auc": roc_auc_score(test["label"], probabilities),
        "classification_report": classification_report(
            test["label"], predictions, output_dict=True, zero_division=0
        ),
        "train_isolates": len(train),
        "test_isolates": len(test),
        "window_size": config["window_size"],
        "maximum_windows_per_isolate": config["kmer_max_windows_per_isolate"],
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
