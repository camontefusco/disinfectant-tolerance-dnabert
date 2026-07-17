#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


BASES = np.array(list("ACGT"))


def centered_fragment(sequence: str, length: int) -> str:
    if len(sequence) <= length:
        return sequence
    start = (len(sequence) - length) // 2
    return sequence[start : start + length]


def noisy_sequence(sequence: str, rate: float, key: str) -> str:
    if not sequence or rate <= 0:
        return sequence
    seed = int(hashlib.sha256(key.encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    output = np.array(list(sequence))
    count = max(1, round(len(sequence) * rate))
    indices = rng.choice(len(output), size=min(count, len(output)), replace=False)
    for index in indices:
        choices = BASES[BASES != output[index]]
        output[index] = rng.choice(choices) if len(choices) else rng.choice(BASES)
    return "".join(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/manifests/megares_similarity_clustered_manifest.csv")
    parser.add_argument("--output-dir", default="data/processed/megares_conditions")
    parser.add_argument("--partial-length", type=int, default=500)
    parser.add_argument("--noise-rate", type=float, default=0.01)
    args = parser.parse_args()

    data = pd.read_csv(args.manifest)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    conditions = {
        "full": data["sequence"],
        f"partial_{args.partial_length}bp": data["sequence"].map(lambda value: centered_fragment(value, args.partial_length)),
        f"noise_{args.noise_rate:.0%}": [
            noisy_sequence(sequence, args.noise_rate, fragment_id)
            for sequence, fragment_id in zip(data["sequence"], data["fragment_id"])
        ],
    }
    for condition, sequences in conditions.items():
        output = pd.DataFrame({"isolate_id": data["fragment_id"], "sequence": sequences})
        path = output_dir / f"{condition}.csv"
        output.to_csv(path, index=False)
        print(f"Wrote {len(output)} {condition} fragments to {path}")


if __name__ == "__main__":
    main()
