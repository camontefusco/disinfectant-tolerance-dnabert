#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator

try:
    from scripts.evaluate_megares_binary_biocide_relevance import select_binary_task, split_data
    from scripts.evaluate_megares_similarity_robustness import load_embeddings, make_model
except ModuleNotFoundError:
    from evaluate_megares_binary_biocide_relevance import select_binary_task, split_data
    from evaluate_megares_similarity_robustness import load_embeddings, make_model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--megares", default="data/manifests/megares_similarity_clustered_manifest.csv")
    parser.add_argument("--megares-embeddings", default="data/processed/megares_embeddings/full.npz")
    parser.add_argument("--candidate-manifest", default="data/manifests/pharma_water_selective_dnabert2_candidates.csv")
    parser.add_argument("--candidate-embeddings", default="data/processed/pharma_water_selective_dnabert2_embeddings.npz")
    parser.add_argument("--evidence", default="data/manifests/pharma_water_pilot_evidence_matrix.csv")
    parser.add_argument("--output", default="data/manifests/pharma_water_pilot_evidence_matrix_with_dnabert2.csv")
    args = parser.parse_args()

    data = select_binary_task(pd.read_csv(args.megares))
    training_embeddings = load_embeddings(Path(args.megares_embeddings)).set_index("fragment_id")
    candidates = pd.read_csv(args.candidate_manifest)
    candidate_embeddings = load_embeddings(Path(args.candidate_embeddings)).set_index("fragment_id")
    probabilities = []
    for seed in [42, 43, 44, 45, 46]:
        train_indices, dev_indices, _ = split_data(data, seed)
        train, dev = data.iloc[train_indices], data.iloc[dev_indices]
        model = make_model("dnabert2_xgboost", seed)
        model.fit(training_embeddings.loc[train["fragment_id"]], train["label"])
        calibrated = CalibratedClassifierCV(FrozenEstimator(model), method="sigmoid")
        calibrated.fit(training_embeddings.loc[dev["fragment_id"]], dev["label"])
        probabilities.append(calibrated.predict_proba(candidate_embeddings.loc[candidates["region_id"]])[:, 1])
    matrix = np.column_stack(probabilities)
    candidates["selective_dnabert2_biocide_relevance_score"] = matrix.mean(axis=1)
    candidates["selective_dnabert2_score_std_across_splits"] = matrix.std(axis=1)
    evidence = pd.read_csv(args.evidence)
    output = evidence.merge(
        candidates[
            [
                "region_id", "aligned_length", "selective_dnabert2_biocide_relevance_score",
                "selective_dnabert2_score_std_across_splits",
            ]
        ],
        on="region_id",
        how="left",
    )
    output.to_csv(args.output, index=False)
    print(f"Wrote {output['selective_dnabert2_biocide_relevance_score'].notna().sum()} selective DNABERT2 scores to {args.output}")


if __name__ == "__main__":
    main()
