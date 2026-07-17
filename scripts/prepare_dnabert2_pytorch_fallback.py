#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import shutil

from huggingface_hub import snapshot_download


OPTIONAL_TRITON_IMPORT = """try:
    from .flash_attn_triton import flash_attn_qkvpacked_func
except ImportError as e:
    flash_attn_qkvpacked_func = None
"""
PYTORCH_FALLBACK = """# Triton flash attention is optional. Use the official PyTorch fallback.
flash_attn_qkvpacked_func = None
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-id", default="zhihan1996/DNABERT-2-117M")
    parser.add_argument("--output-dir", default="data/external/dnabert2_pytorch_fallback")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    snapshot = Path(snapshot_download(repo_id=args.repo_id))
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(snapshot, output_dir, dirs_exist_ok=True, symlinks=False)

    layers_path = output_dir / "bert_layers.py"
    layers = layers_path.read_text()
    if OPTIONAL_TRITON_IMPORT in layers:
        layers_path.write_text(layers.replace(OPTIONAL_TRITON_IMPORT, PYTORCH_FALLBACK))
    elif PYTORCH_FALLBACK not in layers:
        raise RuntimeError("Could not identify the optional Triton import in bert_layers.py")

    flash_path = output_dir / "flash_attn_triton.py"
    if flash_path.exists():
        flash_path.unlink()
    print(f"Prepared DNABERT2 PyTorch fallback snapshot at {output_dir}")


if __name__ == "__main__":
    main()

