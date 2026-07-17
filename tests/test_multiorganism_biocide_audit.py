from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def run_audit(tmp_path: Path) -> tuple[pd.DataFrame, dict]:
    output_csv = tmp_path / "multiorganism.csv"
    output_json = tmp_path / "multiorganism.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "audit_multiorganism_ml_ready_biocide_datasets.py"),
            "--output-csv",
            str(output_csv),
            "--output-json",
            str(output_json),
        ],
        cwd=ROOT,
        check=True,
    )
    return pd.read_csv(output_csv), json.loads(output_json.read_text())


def test_multiorganism_audit_prioritizes_more_than_proteus(tmp_path: Path) -> None:
    audit, report = run_audit(tmp_path)
    immediate = audit[audit["audit_batch"].eq("immediate")]

    assert report["summary"]["immediate_audit_candidates"] >= 4
    assert "proteus_2025_chx_cationic" in immediate["study_key"].tolist()
    assert "enterococcus_faecium_2019_chx" in immediate["study_key"].tolist()
    assert "hospital_environment_chx_2024" in immediate["study_key"].tolist()


def test_multiorganism_audit_keeps_small_or_mechanistic_sets_out_of_training(tmp_path: Path) -> None:
    audit, _ = run_audit(tmp_path)
    not_training_ready = audit[
        audit["study_key"].isin(["serratia_hri_2020_qac", "ecoli_2021_bzk_evolution", "pseudomonas_2021_triclosan"])
    ]

    assert not not_training_ready.empty
    assert not not_training_ready["audit_batch"].eq("immediate").any()
    assert not not_training_ready["ml_readiness"].eq("top_audit_candidate").any()


def test_multiorganism_audit_records_accession_sources_without_inventing_labels(tmp_path: Path) -> None:
    audit, _ = run_audit(tmp_path)
    accession_leads = audit[audit["has_reported_bioproject"]]

    assert len(accession_leads) >= 4
    assert audit.loc[audit["study_key"].eq("food_drain_biofilm_2025"), "accession_or_project"].iloc[0] == "PRJNA1194757"
    assert audit.loc[audit["study_key"].eq("ecoli_2021_bzk_evolution"), "row_level_labels_status"].iloc[0] == "lab_evolution_design"
