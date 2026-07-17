#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.splits import SplitFractions, assign_grouped_splits
try:
    from scripts.evaluate_embedding_repeated_grouped_splits import best_threshold
except ModuleNotFoundError:
    from evaluate_embedding_repeated_grouped_splits import best_threshold


def aggregate_embeddings(
    isolate_ids: np.ndarray,
    embeddings: np.ndarray,
    cap: int,
    aggregation: str,
) -> pd.DataFrame:
    rows = []
    for isolate_id in sorted(set(isolate_ids.astype(str))):
        vectors = embeddings[isolate_ids.astype(str) == isolate_id]
        if len(vectors) > cap:
            indices = np.linspace(0, len(vectors) - 1, cap).round().astype(int)
            vectors = vectors[indices]
        if aggregation == "mean":
            vector = vectors.mean(axis=0)
        elif aggregation == "median":
            vector = np.median(vectors, axis=0)
        else:
            raise ValueError(f"Unsupported aggregation: {aggregation}")
        rows.append([isolate_id, *vector])
    return pd.DataFrame(rows, columns=["isolate_id", *range(embeddings.shape[1])])


def make_model(classifier: str, seed: int):
    if classifier == "logistic_regression":
        return Pipeline(
            [
                ("scaler", StandardScaler()),
                ("classifier", LogisticRegression(class_weight="balanced", max_iter=2_000, random_state=seed)),
            ]
        )
    if classifier == "linear_svm":
        return Pipeline(
            [
                ("scaler", StandardScaler()),
                ("classifier", LinearSVC(class_weight="balanced", random_state=seed, dual="auto", max_iter=10_000)),
            ]
        )
    if classifier == "random_forest":
        return RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=seed, n_jobs=-1)
    raise ValueError(f"Unsupported classifier: {classifier}")


def scores(model, frame: pd.DataFrame, feature_columns: list) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(frame[feature_columns])[:, 1]
    return model.decision_function(frame[feature_columns])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--window-embeddings", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--caps", default="16,32,64")
    parser.add_argument("--aggregations", default="mean,median")
    parser.add_argument("--classifiers", default="logistic_regression,linear_svm,random_forest")
    parser.add_argument("--seeds", default="42,43,44,45,46,47,48,49,50,51")
    parser.add_argument("--randomization-strength", type=float, default=0.2)
    parser.add_argument("--output", default="reports/dnabert2_strategy_evaluation.json")
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    manifest = pd.read_csv(args.manifest, dtype={"label": int}).fillna("")
    archive = np.load(args.window_embeddings)
    isolate_ids = archive["isolate_ids"].astype(str)
    embeddings = archive["embeddings"]
    fractions = SplitFractions(**config["split_fractions"])
    row_records = manifest[["isolate_id", "label", "group"]].astype(str).to_dict(orient="records")
    seeds = [int(value) for value in args.seeds.split(",")]
    results = []

    for cap in [int(value) for value in args.caps.split(",")]:
        for aggregation in args.aggregations.split(","):
            vectors = aggregate_embeddings(isolate_ids, embeddings, cap, aggregation)
            feature_columns = list(vectors.columns[1:])
            data = manifest[["isolate_id", "label", "group"]].merge(vectors, on="isolate_id", validate="one_to_one")
            for classifier in args.classifiers.split(","):
                for seed in seeds:
                    assignments = assign_grouped_splits(
                        row_records,
                        fractions=fractions,
                        seed=seed,
                        randomization_strength=args.randomization_strength,
                    )
                    split_data = data.copy()
                    split_data["split"] = split_data["group"].map(assignments)
                    train = split_data[split_data["split"] == "train"]
                    dev = split_data[split_data["split"] == "dev"]
                    test = split_data[split_data["split"] == "test"]
                    model = make_model(classifier, seed)
                    model.fit(train[feature_columns], train["label"])
                    threshold = best_threshold(dev["label"], scores(model, dev, feature_columns))
                    probabilities = scores(model, test, feature_columns)
                    predictions = (probabilities >= threshold).astype(int)
                    results.append(
                        {
                            "cap": cap,
                            "aggregation": aggregation,
                            "classifier": classifier,
                            "seed": seed,
                            "threshold_selected_on_dev": threshold,
                            "balanced_accuracy": balanced_accuracy_score(test["label"], predictions),
                            "roc_auc": roc_auc_score(test["label"], probabilities),
                        }
                    )
                print(f"Finished cap={cap}, aggregation={aggregation}, classifier={classifier}")

    frame = pd.DataFrame(results)
    summary = (
        frame.groupby(["cap", "aggregation", "classifier"])[["balanced_accuracy", "roc_auc"]]
        .agg(["mean", "std", "min", "max"])
        .reset_index()
    )
    summary.columns = ["_".join(str(part) for part in column if part) for column in summary.columns]
    report = {"results": results, "summary": summary.to_dict(orient="records")}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(summary.sort_values("roc_auc_mean", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
