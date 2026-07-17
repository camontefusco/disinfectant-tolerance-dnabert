#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, log_loss
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
CONDITIONS = ("full", "partial_500bp", "noise_1%")
CLASSES = ("Biocides", "Metals", "Multi-compound")


def brier_score(labels: np.ndarray, probabilities: np.ndarray, class_count: int) -> float:
    expected = np.eye(class_count)[labels]
    return float(np.mean(np.sum((probabilities - expected) ** 2, axis=1)))


def expected_calibration_error(labels: np.ndarray, probabilities: np.ndarray, bins: int = 10) -> float:
    predictions = probabilities.argmax(axis=1)
    confidence = probabilities.max(axis=1)
    correct = predictions == labels
    total = len(labels)
    error = 0.0
    for lower in np.linspace(0.0, 1.0, bins, endpoint=False):
        upper = lower + 1.0 / bins
        selected = (confidence >= lower) & (confidence < upper if upper < 1.0 else confidence <= upper)
        if selected.any():
            error += selected.sum() / total * abs(correct[selected].mean() - confidence[selected].mean())
    return float(error)


def metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = probabilities.argmax(axis=1)
    return {
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "macro_f1": float(f1_score(labels, predictions, average="macro")),
        "log_loss": float(log_loss(labels, probabilities, labels=np.arange(len(CLASSES)))),
        "brier_score": brier_score(labels, probabilities, len(CLASSES)),
        "expected_calibration_error": expected_calibration_error(labels, probabilities),
    }


def make_model(name: str, seed: int):
    if "xgboost" in name:
        from xgboost import XGBClassifier

    if name == "kmer_logistic_regression":
        return Pipeline(
            [
                ("tfidf", TfidfVectorizer(analyzer="char", ngram_range=(4, 6), lowercase=False, max_features=100_000)),
                ("classifier", LogisticRegression(class_weight="balanced", max_iter=2_000, random_state=seed)),
            ]
        )
    if name in {"kmer_random_forest", "kmer_xgboost"}:
        classifier = (
            RandomForestClassifier(n_estimators=300, class_weight="balanced", n_jobs=2, random_state=seed)
            if name == "kmer_random_forest"
            else XGBClassifier(n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.9, colsample_bytree=0.9, n_jobs=2, random_state=seed)
        )
        return Pipeline(
            [
                ("tfidf", TfidfVectorizer(analyzer="char", ngram_range=(4, 6), lowercase=False, max_features=50_000)),
                ("svd", TruncatedSVD(n_components=128, random_state=seed)),
                ("classifier", classifier),
            ]
        )
    if name == "dnabert2_logistic_regression":
        return Pipeline(
            [
                ("scaler", StandardScaler()),
                ("classifier", LogisticRegression(class_weight="balanced", max_iter=2_000, random_state=seed)),
            ]
        )
    if name == "dnabert2_random_forest":
        return RandomForestClassifier(n_estimators=300, class_weight="balanced", n_jobs=2, random_state=seed)
    if name == "dnabert2_xgboost":
        return XGBClassifier(n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.9, colsample_bytree=0.9, n_jobs=2, random_state=seed)
    raise ValueError(f"Unknown model: {name}")


def split_data(data: pd.DataFrame, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    for offset in range(100):
        outer = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed + offset)
        train_dev, test = next(outer.split(data, data["target"], groups=data["similarity_cluster"]))
        inner_data = data.iloc[train_dev]
        inner = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed + 10_000 + offset)
        train_relative, dev_relative = next(inner.split(inner_data, inner_data["target"], groups=inner_data["similarity_cluster"]))
        train = train_dev[train_relative]
        dev = train_dev[dev_relative]
        if all(data.iloc[indices]["target"].nunique() == len(CLASSES) for indices in (train, dev, test)):
            return train, dev, test
    raise ValueError(f"Could not produce a three-class similarity-aware split for seed {seed}")


def load_embeddings(path: Path) -> pd.DataFrame:
    archive = np.load(path)
    frame = pd.DataFrame(archive["embeddings"])
    frame.insert(0, "fragment_id", archive["isolate_ids"].astype(str))
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/manifests/megares_similarity_clustered_manifest.csv")
    parser.add_argument("--conditions-dir", default="data/processed/megares_conditions")
    parser.add_argument("--embeddings-dir", default="data/processed/megares_embeddings")
    parser.add_argument("--seeds", default="42,43,44,45,46")
    parser.add_argument("--output", default="reports/megares_similarity_robustness_benchmark.json")
    args = parser.parse_args()

    data = pd.read_csv(args.manifest)
    data["target"] = data["type"].where(data["type"].isin(CLASSES), "Multi-compound")
    encoder = LabelEncoder().fit(CLASSES)
    data["label"] = encoder.transform(data["target"])
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
        train = data.iloc[train_indices]
        dev = data.iloc[dev_indices]
        test = data.iloc[test_indices]
        for model_name in model_names:
            is_embedding = model_name.startswith("dnabert2")
            if is_embedding:
                full = embedding_conditions["full"]
                train_x = full.loc[train["fragment_id"]]
                dev_x = full.loc[dev["fragment_id"]]
            else:
                full = sequence_conditions["full"]
                train_x = full.loc[train["fragment_id"]]
                dev_x = full.loc[dev["fragment_id"]]
            model = make_model(model_name, seed)
            model.fit(train_x, train["label"])
            calibrated = CalibratedClassifierCV(FrozenEstimator(model), method="sigmoid")
            calibrated.fit(dev_x, dev["label"])
            for condition in CONDITIONS:
                condition_data = embedding_conditions[condition] if is_embedding else sequence_conditions[condition]
                test_x = condition_data.loc[test["fragment_id"]]
                raw = metrics(test["label"].to_numpy(), model.predict_proba(test_x))
                adjusted = metrics(test["label"].to_numpy(), calibrated.predict_proba(test_x))
                evaluations.append(
                    {
                        "seed": seed,
                        "model": model_name,
                        "condition": condition,
                        "train_fragments": len(train),
                        "dev_fragments": len(dev),
                        "test_fragments": len(test),
                        "train_similarity_clusters": int(train["similarity_cluster"].nunique()),
                        "dev_similarity_clusters": int(dev["similarity_cluster"].nunique()),
                        "test_similarity_clusters": int(test["similarity_cluster"].nunique()),
                        "raw": raw,
                        "calibrated": adjusted,
                    }
                )
            print(f"Finished seed={seed}, model={model_name}")
    flat_rows = []
    for row in evaluations:
        flat_rows.append(
            {
                "seed": row["seed"], "model": row["model"], "condition": row["condition"],
                **{f"raw_{key}": value for key, value in row["raw"].items()},
                **{f"calibrated_{key}": value for key, value in row["calibrated"].items()},
            }
        )
    frame = pd.DataFrame(flat_rows)
    summary = frame.groupby(["model", "condition"]).agg(["mean", "std"]).reset_index()
    summary.columns = ["_".join(str(value) for value in column if value) for column in summary.columns]
    result = {
        "claim_boundary": "Exploratory resistance-associated fragment classification, not isolate-level sanitizer-survival prediction.",
        "similarity_split": "Connected components at >=90% nucleotide identity and >=80% bidirectional coverage.",
        "conditions": list(CONDITIONS),
        "evaluations": evaluations,
        "summary": summary.to_dict(orient="records"),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(summary.sort_values(["condition", "calibrated_balanced_accuracy_mean"], ascending=[True, False]).to_string(index=False))


if __name__ == "__main__":
    main()
