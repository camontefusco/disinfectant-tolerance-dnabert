#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer


def hidden_state(model_output: object) -> torch.Tensor:
    if hasattr(model_output, "last_hidden_state"):
        return model_output.last_hidden_state
    if isinstance(model_output, tuple):
        return model_output[0]
    raise TypeError(f"Unsupported model output type: {type(model_output).__name__}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows", default="data/processed/windows.csv")
    parser.add_argument("--config", default="configs/feasibility.json")
    parser.add_argument("--output", default="data/processed/dnabert2_mean_embeddings.npz")
    parser.add_argument("--window-output")
    parser.add_argument("--model")
    parser.add_argument("--limit-isolates", type=int)
    parser.add_argument("--maximum-windows-per-isolate", type=int)
    parser.add_argument("--batch-size", type=int)
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    data = pd.read_csv(args.windows)
    if args.limit_isolates is not None:
        isolate_ids = data["isolate_id"].drop_duplicates().head(args.limit_isolates)
        data = data[data["isolate_id"].isin(isolate_ids)]
    maximum_windows = (
        args.maximum_windows_per_isolate
        if args.maximum_windows_per_isolate is not None
        else config["maximum_windows_per_isolate"]
    )
    data = data.groupby("isolate_id", group_keys=False).head(
        maximum_windows
    )
    model_name = args.model or config["dnabert_model"]
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
    if torch.cuda.is_available():
        device = "cuda"
    elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        device = "mps"
    else:
        device = "cpu"
    model.to(device).eval()
    print(
        f"Embedding {data['isolate_id'].nunique()} isolates and {len(data)} windows "
        f"on {device}"
    )

    embeddings: dict[str, list[np.ndarray]] = defaultdict(list)
    batch_size = args.batch_size or config["dnabert_batch_size"]
    completed_windows = 0
    with torch.no_grad():
        for start in range(0, len(data), batch_size):
            batch = data.iloc[start : start + batch_size]
            tokens = tokenizer(
                batch["sequence"].tolist(),
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max(256, config["window_size"]),
            ).to(device)
            hidden = hidden_state(model(**tokens))
            mask = tokens["attention_mask"].unsqueeze(-1)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
            for isolate_id, embedding in zip(batch["isolate_id"], pooled.cpu().numpy()):
                embeddings[isolate_id].append(embedding)
            completed_windows += len(batch)
            if completed_windows % 100 == 0 or completed_windows == len(data):
                print(f"Embedded {completed_windows}/{len(data)} windows")

    isolate_ids = sorted(embeddings)
    matrix = np.stack([np.mean(embeddings[key], axis=0) for key in isolate_ids])
    np.savez_compressed(args.output, isolate_ids=np.array(isolate_ids), embeddings=matrix)
    print(f"Wrote {len(isolate_ids)} isolate embeddings to {args.output}")
    if args.window_output:
        window_isolate_ids = []
        window_embeddings = []
        for isolate_id in isolate_ids:
            for embedding in embeddings[isolate_id]:
                window_isolate_ids.append(isolate_id)
                window_embeddings.append(embedding)
        np.savez_compressed(
            args.window_output,
            isolate_ids=np.array(window_isolate_ids),
            embeddings=np.stack(window_embeddings),
        )
        print(f"Wrote {len(window_embeddings)} window embeddings to {args.window_output}")


if __name__ == "__main__":
    main()
