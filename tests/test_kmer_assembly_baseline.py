from pathlib import Path

import pandas as pd

from scripts.train_kmer_baseline_from_assemblies import assembly_text, attach_splits


def test_assembly_manifest_can_keep_existing_split(tmp_path: Path) -> None:
    manifest = pd.DataFrame(
        [
            {"isolate_id": "a", "assembly_path": "a.fasta", "label": 0, "split": "train"},
            {"isolate_id": "b", "assembly_path": "b.fasta", "label": 1, "split": "test"},
        ]
    )
    conflicting_splits = pd.DataFrame(
        [
            {"isolate_id": "a", "split": "test"},
            {"isolate_id": "b", "split": "train"},
        ]
    )
    result = attach_splits(manifest, conflicting_splits)
    assert list(result["split"]) == ["train", "test"]


def test_assembly_text_caps_evenly_spaced_windows(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fasta"
    fasta.write_text(">contig\n" + "A" * 20 + "C" * 20 + "\n")
    text = assembly_text(fasta, window_size=10, maximum_windows=3)
    assert text == "A" * 10 + "N" + "A" * 5 + "C" * 5 + "N" + "C" * 10
