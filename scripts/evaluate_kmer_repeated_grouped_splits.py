#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from disinfectant_tolerance.splits import SplitFractions, assign_grouped_splits
try:
    from scripts.evaluate_embedding_repeated_grouped_splits import best_threshold
    from scripts.train_kmer_baseline_from_assemblies import assembly_text
except ModuleNotFoundError:
    from evaluate_embedding_repeated_grouped_splits import best_threshold
    from train_kmer_baseline_from_assemblies import assembly_text


def make_model(low: int, high: int, seed: int) -> Pipeline:
    return Pipeline(
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
                    random_state=seed,
                ),
            ),
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--seeds", default="42,43,44,45,46,47,48,49,50,51")
    parser.add_argument("--randomization-strength", type=float, default=0.2)
    parser.add_argument("--output", default="reports/kmer_repeated_grouped_evaluation.json")
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    manifest = pd.read_csv(args.manifest, dtype={"label": int}).fillna("")
    data = manifest[["isolate_id", "label", "group", "assembly_path"]].copy()
    data["sequence"] = data["assembly_path"].map(
        lambda path: assembly_text(
            path,
            window_size=config["window_size"],
            maximum_windows=config["kmer_max_windows_per_isolate"],
        )
    )

    fractions = SplitFractions(**config["split_fractions"])
    rows = manifest[["isolate_id", "label", "group"]].astype(str).to_dict(orient="records")
    low, high = config["kmer_range"]

    evaluations = []
    for seed in (int(value) for value in args.seeds.split(",")):
        assignments = assign_grouped_splits(
            rows,
            fractions=fractions,
            seed=seed,
            randomization_strength=args.randomization_strength,
        )
        split_data = data.copy()
        split_data["split"] = split_data["group"].map(assignments)
        train = split_data[split_data["split"] == "train"]
        dev = split_data[split_data["split"] == "dev"]
        test = split_data[split_data["split"] == "test"]
        if min(train["label"].nunique(), dev["label"].nunique(), test["label"].nunique()) < 2:
            raise ValueError(f"Seed {seed} produced a split without both labels")

        model = make_model(low, high, seed)
        model.fit(train["sequence"], train["label"])
        dev_probabilities = model.predict_proba(dev["sequence"])[:, 1]
        threshold = best_threshold(dev["label"], dev_probabilities)
        test_probabilities = model.predict_proba(test["sequence"])[:, 1]
        test_predictions = (test_probabilities >= threshold).astype(int)
        evaluations.append(
            {
                "seed": seed,
                "threshold_selected_on_dev": threshold,
                "balanced_accuracy": balanced_accuracy_score(test["label"], test_predictions),
                "roc_auc": roc_auc_score(test["label"], test_probabilities),
                "train_isolates": len(train),
                "dev_isolates": len(dev),
                "test_isolates": len(test),
                "test_sensitive": int((test["label"] == 0).sum()),
                "test_tolerant": int((test["label"] == 1).sum()),
            }
        )
        print(f"Finished seed {seed}")

    result_frame = pd.DataFrame(evaluations)
    summary = {
        metric: {
            "mean": float(result_frame[metric].mean()),
            "std": float(result_frame[metric].std(ddof=1)),
            "min": float(result_frame[metric].min()),
            "max": float(result_frame[metric].max()),
        }
        for metric in ("balanced_accuracy", "roc_auc")
    }
    result = {
        "model": "capped_assembly_character_kmer_tfidf_logistic_regression",
        "split_strategy": "repeated_lineage_grouped_train_dev_test",
        "randomization_strength": args.randomization_strength,
        "window_size": config["window_size"],
        "maximum_windows_per_isolate": config["kmer_max_windows_per_isolate"],
        "summary": summary,
        "evaluations": evaluations,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
