from __future__ import annotations

import pandas as pd

from scripts.audit_bcc_combined_phenotype_genome import (
    build_candidate_manifest,
    build_report,
    extract_kim_2015,
    extract_moore_2009,
)


def test_extract_moore_2009_row_level_values() -> None:
    table = pd.DataFrame(
        [
            ["B. cepacia (7)", None, None, None, None, None, None],
            ["LMG 17997", "100", "100", "70-100", "2000", "400", "200"],
        ]
    )
    rows = extract_moore_2009([pd.DataFrame(), table])

    bzk = rows[(rows["strain_id"].eq("LMG 17997")) & rows["compound"].eq("benzalkonium chloride")]

    assert len(rows) == 6
    assert bzk.iloc[0]["value"] == "200"
    assert bzk.iloc[0]["species_reported"] == "B. cepacia"


def test_extract_kim_2015_row_level_values() -> None:
    columns = pd.MultiIndex.from_tuples(
        [
            ("Species name", "Species name", "Species name"),
            ("Strain number", "Strain number", "Strain number"),
            ("Isolation source", "Isolation source", "Isolation source"),
            ("CHX (ug/ml)", "x", "Initial"),
            ("CHX (ug/ml)", "x", "40"),
            ("BZK (ug/ml)", "x", "Initial"),
            ("BZK (ug/ml)", "x", "40"),
        ]
    )
    table = pd.DataFrame([["B. cepacia", "PC783", "Onion", "100", "100", "50", "50"]], columns=columns)
    rows = extract_kim_2015([table])

    assert len(rows) == 4
    assert set(rows["compound"]) == {"chlorhexidine", "benzalkonium chloride"}
    assert set(rows["timepoint"]) == {"baseline", "day_40"}


def test_candidate_manifest_requires_match_flags() -> None:
    phenotypes = pd.DataFrame(
        [
            {
                "study_key": "x",
                "strain_id": "J2315",
                "compound": "benzalkonium chloride",
                "endpoint": "MIC",
                "value": "350",
            }
        ]
    )
    matches = pd.DataFrame(
        [
            {
                "query_strain": "J2315",
                "assembly_accession": "GCF_000009485.1",
                "exact_norm_match": True,
                "contains_norm_match": True,
            }
        ]
    )

    manifest = build_candidate_manifest(phenotypes, matches)
    report = build_report(phenotypes, matches, manifest)

    assert len(manifest) == 1
    assert report["candidate_manifest_unique_strains"] == 1
