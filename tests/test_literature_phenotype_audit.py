from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_literature_audit_identifies_only_one_highest_priority_candidate(tmp_path: Path) -> None:
    output_csv = tmp_path / "audit.csv"
    output_json = tmp_path / "audit.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "audit_literature_phenotype_datasets.py"),
            "--output-csv",
            str(output_csv),
            "--output-json",
            str(output_json),
        ],
        cwd=ROOT,
        check=True,
    )

    audit = pd.read_csv(output_csv)
    report = json.loads(output_json.read_text())
    highest = audit[audit["supervised_ml_readiness"].eq("highest_priority_candidate")]

    assert len(highest) == 1
    assert highest.iloc[0]["study_key"] == "pottier_2023_paeruginosa_ddac"
    assert highest.iloc[0]["accession_or_project"].find("PRJNA884650") >= 0
    assert report["summary"]["highest_priority_candidates"] == 1


def test_literature_audit_does_not_mark_bcc_without_genomes_as_ready(tmp_path: Path) -> None:
    output_csv = tmp_path / "audit.csv"
    output_json = tmp_path / "audit.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "audit_literature_phenotype_datasets.py"),
            "--output-csv",
            str(output_csv),
            "--output-json",
            str(output_json),
        ],
        cwd=ROOT,
        check=True,
    )

    audit = pd.read_csv(output_csv)
    bcc = audit[audit["organism_group"].str.contains("Burkholderia cepacia complex", regex=False)]

    assert not bcc.empty
    assert not bcc["supervised_ml_readiness"].eq("highest_priority_candidate").any()
    assert bcc["genome_extract_status"].isin(["not_confirmed"]).any()
