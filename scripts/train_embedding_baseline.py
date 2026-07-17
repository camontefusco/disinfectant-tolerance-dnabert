#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--embeddings", default="data/processed/dnabert2_mean_embeddings.npz")
    parser.add_argument("--splits", default="data/processed/isolate_splits.csv")
    parser.add_argument("--output", default="reports/dnabert2_embedding_baseline.json")
    args = parser.parse_args()

    archive = np.load(args.embeddings)
    vectors = pd.DataFrame(archive["embeddings"])
    vectors.insert(0, "isolate_id", archive["isolate_ids"].astype(str))
    splits = pd.read_csv(args.splits)[["isolate_id", "label", "split"]]
    data = splits.merge(vectors, on="isolate_id", validate="one_to_one")
    train = data[data["split"] == "train"]
    test = data[data["split"] == "test"]
    feature_columns = vectors.columns[1:]
    if train.empty or test.empty:
        raise ValueError("Train and test splits must both contain isolates")

    model = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("classifier", LogisticRegression(class_weight="balanced", max_iter=2_000)),
        ]
    )
    model.fit(train[feature_columns], train["label"])
    probabilities = model.predict_proba(test[feature_columns])[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    result = {
        "model": "frozen_dnabert2_mean_embeddings_logistic_regression",
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

