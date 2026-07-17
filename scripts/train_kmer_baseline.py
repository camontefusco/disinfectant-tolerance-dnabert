#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, classification_report, roc_auc_score
from sklearn.pipeline import Pipeline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows", default="data/processed/windows.csv")
    parser.add_argument("--splits", default="data/processed/isolate_splits.csv")
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--output", default="reports/kmer_baseline.json")
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    windows = pd.read_csv(args.windows)
    splits = pd.read_csv(args.splits)[["isolate_id", "split"]]
    data = windows.merge(splits, on="isolate_id", validate="many_to_one")
    data = (
        data.groupby(["isolate_id", "label", "split"], as_index=False)["sequence"]
        .agg("N".join)
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
        "model": "character_kmer_tfidf_logistic_regression",
        "balanced_accuracy": balanced_accuracy_score(test["label"], predictions),
        "roc_auc": roc_auc_score(test["label"], probabilities),
        "classification_report": classification_report(
            test["label"], predictions, output_dict=True, zero_division=0
        ),
        "train_isolates": len(train),
        "test_isolates": len(test),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

