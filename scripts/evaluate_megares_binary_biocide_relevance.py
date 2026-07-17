#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, f1_score, log_loss, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

try:
    from scripts.evaluate_megares_similarity_robustness import (
        CONDITIONS,
        expected_calibration_error,
        load_embeddings,
        make_model,
    )
except ModuleNotFoundError:
    from evaluate_megares_similarity_robustness import (
        CONDITIONS,
        expected_calibration_error,
        load_embeddings,
        make_model,
    )


def select_binary_task(data: pd.DataFrame) -> pd.DataFrame:
    positive = data["type"].eq("Biocides") | (
        data["type"].eq("Multi-compound") & data["class"].str.contains("biocide", case=False, na=False)
    )
    negative = data["type"].eq("Metals")
    selected = data[positive | negative].copy()
    selected["target"] = np.where(positive.loc[selected.index], "biocide_relevant", "metal_only")
    selected["label"] = selected["target"].eq("biocide_relevant").astype(int)
    return selected


def split_data(data: pd.DataFrame, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    for offset in range(100):
        outer = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed + offset)
        train_dev, test = next(outer.split(data, data["label"], groups=data["similarity_cluster"]))
        inner_data = data.iloc[train_dev]
        inner = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed + 10_000 + offset)
        train_relative, dev_relative = next(inner.split(inner_data, inner_data["label"], groups=inner_data["similarity_cluster"]))
        train, dev = train_dev[train_relative], train_dev[dev_relative]
        if all(data.iloc[indices]["label"].nunique() == 2 for indices in (train, dev, test)):
            return train, dev, test
    raise ValueError(f"Could not produce a two-class similarity-aware split for seed {seed}")


def metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    positive = probabilities[:, 1]
    predictions = (positive >= 0.5).astype(int)
    return {
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "f1": float(f1_score(labels, predictions)),
        "roc_auc": float(roc_auc_score(labels, positive)),
        "average_precision": float(average_precision_score(labels, positive)),
        "log_loss": float(log_loss(labels, probabilities, labels=[0, 1])),
        "brier_score": float(brier_score_loss(labels, positive)),
        "expected_calibration_error": expected_calibration_error(labels, probabilities),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/manifests/megares_similarity_clustered_manifest.csv")
    parser.add_argument("--conditions-dir", default="data/processed/megares_conditions")
    parser.add_argument("--embeddings-dir", default="data/processed/megares_embeddings")
    parser.add_argument("--seeds", default="42,43,44,45,46")
    parser.add_argument("--output", default="reports/megares_binary_biocide_relevance_benchmark.json")
    args = parser.parse_args()

    data = select_binary_task(pd.read_csv(args.manifest))
    sequence_conditions = {
        condition: pd.read_csv(Path(args.conditions_dir) / f"{condition}.csv").set_index("isolate_id")["sequence"]
        for condition in CONDITIONS
    }
    embedding_conditions = {
        condition: load_embeddings(Path(args.embeddings_dir) / f"{condition}.npz").set_index("fragment_id")
        for condition in CONDITIONS
    }
    model_names = (
        "kmer_logistic_regression", "kmer_random_forest", "kmer_xgboost",
        "dnabert2_logistic_regression", "dnabert2_random_forest", "dnabert2_xgboost",
    )
    evaluations = []
    for seed in [int(value) for value in args.seeds.split(",")]:
        train_indices, dev_indices, test_indices = split_data(data, seed)
        train, dev, test = data.iloc[train_indices], data.iloc[dev_indices], data.iloc[test_indices]
        for model_name in model_names:
            is_embedding = model_name.startswith("dnabert2")
            source = embedding_conditions["full"] if is_embedding else sequence_conditions["full"]
            model = make_model(model_name, seed)
            model.fit(source.loc[train["fragment_id"]], train["label"])
            calibrated = CalibratedClassifierCV(FrozenEstimator(model), method="sigmoid")
            calibrated.fit(source.loc[dev["fragment_id"]], dev["label"])
            for condition in CONDITIONS:
                condition_data = embedding_conditions[condition] if is_embedding else sequence_conditions[condition]
                test_x = condition_data.loc[test["fragment_id"]]
                evaluations.append(
                    {
                        "seed": seed, "model": model_name, "condition": condition,
                        "train_fragments": len(train), "dev_fragments": len(dev), "test_fragments": len(test),
                        "raw": metrics(test["label"].to_numpy(), model.predict_proba(test_x)),
                        "calibrated": metrics(test["label"].to_numpy(), calibrated.predict_proba(test_x)),
                    }
                )
            print(f"Finished seed={seed}, model={model_name}")
    flat = pd.DataFrame(
        [
            {
                "seed": row["seed"], "model": row["model"], "condition": row["condition"],
                **{f"raw_{key}": value for key, value in row["raw"].items()},
                **{f"calibrated_{key}": value for key, value in row["calibrated"].items()},
            }
            for row in evaluations
        ]
    )
    summary = flat.groupby(["model", "condition"]).agg(["mean", "std"]).reset_index()
    summary.columns = ["_".join(str(value) for value in column if value) for column in summary.columns]
    result = {
        "claim_boundary": "Exploratory biocide-relevant fragment prioritization, not isolate-level tolerance prediction.",
        "task": {
            "positive": "MEGARes Biocides plus Multi-compound classes containing biocide",
            "negative": "MEGARes Metals only",
            "excluded": "Drug and metal resistance Multi-compound records without biocide annotation",
            "class_counts": data["target"].value_counts().to_dict(),
        },
        "similarity_split": "Connected components at >=90% nucleotide identity and >=80% bidirectional coverage.",
        "evaluations": evaluations,
        "summary": summary.to_dict(orient="records"),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(summary.sort_values(["condition", "calibrated_balanced_accuracy_mean"], ascending=[True, False]).to_string(index=False))


if __name__ == "__main__":
    main()
