from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from scripts.extract_proteus_2025_manifest import (
    build_manifest,
    infer_header_row,
    is_valid_xlsx,
    read_normalized_sheet,
)


ROOT = Path(__file__).resolve().parents[1]


def test_infer_header_row_handles_preface_rows() -> None:
    raw = pd.DataFrame(
        [
            ["Table S1", None, None],
            ["notes", None, None],
            ["Isolate ID", "SRA accession", "CHD MIC"],
            ["PM001", "SRR1", 32],
        ]
    )

    assert infer_header_row(raw) == 2


def test_build_manifest_extracts_accessions_and_phenotypes() -> None:
    frame = pd.DataFrame(
        {
            "Isolate ID": ["PM001", "PM002"],
            "SRA accession": ["SRR000001", "SRR000002"],
            "BioSample accession": ["SAMN1", "SAMN2"],
            "CHD MIC": [32, 128],
            "Octenidine MIC": [16, ""],
        }
    )

    isolates, phenotypes = build_manifest(frame, "Table S1")

    assert isolates["isolate_id"].nunique() == 2
    assert "sra_accession" in isolates.columns
    assert len(phenotypes) == 3
    assert set(phenotypes["compound_or_endpoint"]) == {"CHD MIC", "Octenidine MIC"}


def test_script_reports_manual_download_when_no_valid_xlsx(tmp_path: Path) -> None:
    fake = tmp_path / "fake.xlsx"
    fake.write_text("<html>Preparing to download ...</html>")
    report = tmp_path / "report.json"

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "extract_proteus_2025_manifest.py"),
            "--supplement",
            str(fake),
            "--sheet-inventory",
            str(tmp_path / "sheets.csv"),
            "--isolate-output",
            str(tmp_path / "isolates.csv"),
            "--phenotype-output",
            str(tmp_path / "phenotypes.csv"),
            "--report-output",
            str(report),
        ],
        cwd=ROOT,
        check=True,
    )

    payload = json.loads(report.read_text())
    assert payload["decision"] == "blocked_manual_supplement_download_required"
    assert payload["valid_supplement_found"] is False


def test_script_extracts_from_valid_synthetic_workbook(tmp_path: Path) -> None:
    workbook = tmp_path / "table_s1.xlsx"
    pd.DataFrame(
        [
            ["Supplementary Table S1", None, None, None],
            ["Isolate ID", "SRA accession", "BioSample accession", "CHD MIC"],
            ["PM001", "SRR000001", "SAMN1", 64],
            ["PM002", "SRR000002", "SAMN2", 128],
        ]
    ).to_excel(workbook, index=False, header=False)

    assert is_valid_xlsx(workbook)
    frame = read_normalized_sheet(workbook, "Sheet1")
    assert "Isolate ID" in frame.columns

    report = tmp_path / "report.json"
    phenotypes = tmp_path / "phenotypes.csv"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "extract_proteus_2025_manifest.py"),
            "--supplement",
            str(workbook),
            "--sheet-inventory",
            str(tmp_path / "sheets.csv"),
            "--isolate-output",
            str(tmp_path / "isolates.csv"),
            "--phenotype-output",
            str(phenotypes),
            "--report-output",
            str(report),
        ],
        cwd=ROOT,
        check=True,
    )

    payload = json.loads(report.read_text())
    extracted = pd.read_csv(phenotypes)
    assert payload["decision"] == "candidate_manifest_extracted"
    assert payload["candidate_isolates"] == 2
    assert len(extracted) == 2
