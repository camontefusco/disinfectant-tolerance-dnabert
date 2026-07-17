#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, classification_report
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline


def normalized_target(value: str) -> str:
    return value if value in {"Biocides", "Metals"} else "Multi-compound"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/manifests/megares_biocide_metal_manifest.csv")
    parser.add_argument("--output", default="reports/megares_homology_aware_kmer_baseline.json")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    data = pd.read_csv(args.manifest)
    data["target"] = data["type"].map(normalized_target)
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=args.seed)
    train_indices, test_indices = next(splitter.split(data, data["target"], groups=data["homology_group"]))
    train = data.iloc[train_indices]
    test = data.iloc[test_indices]
    model = Pipeline(
        [
            ("tfidf", TfidfVectorizer(analyzer="char", ngram_range=(4, 6), lowercase=False, max_features=100_000)),
            ("classifier", LogisticRegression(class_weight="balanced", max_iter=2_000, random_state=args.seed)),
        ]
    )
    model.fit(train["sequence"], train["target"])
    predictions = model.predict(test["sequence"])
    result = {
        "model": "megares_fragment_character_kmer_tfidf_logistic_regression",
        "split_strategy": "group_shuffle_split_by_megares_gene_group",
        "claim_boundary": "Exploratory fragment-family classification, not isolate-level sanitizer-survival prediction.",
        "train_fragments": len(train),
        "test_fragments": len(test),
        "train_groups": int(train["homology_group"].nunique()),
        "test_groups": int(test["homology_group"].nunique()),
        "balanced_accuracy": balanced_accuracy_score(test["target"], predictions),
        "classification_report": classification_report(test["target"], predictions, output_dict=True, zero_division=0),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
