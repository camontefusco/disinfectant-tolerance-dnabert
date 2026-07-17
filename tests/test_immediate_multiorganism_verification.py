from __future__ import annotations

import pandas as pd

from scripts.verify_immediate_multiorganism_candidates import (
    build_report,
    extract_bioprojects,
    find_supplement_links,
    term_hits,
)


def test_extract_bioprojects_deduplicates_projects() -> None:
    value = "Data in PRJNA1154625 and PRJNA475751; repeated PRJNA1154625."

    assert extract_bioprojects(value) == ["PRJNA1154625", "PRJNA475751"]


def test_find_supplement_links_flags_machine_readable_files() -> None:
    html = """
    <html><body>
      <a href="/articles/supplementary_table_s1.xlsx">Table S1</a>
      <a href="/articles/figure_s2.tif">Supplementary Figure S2</a>
      <a href="https://github.com/example/project">Code</a>
    </body></html>
    """

    links = find_supplement_links(html, "https://example.org/paper/")

    assert len(links) == 3
    assert any(row["spreadsheet_candidate"] for row in links)
    assert any("github.com" in str(row["href"]) for row in links)


def test_term_hits_is_case_insensitive() -> None:
    hits = term_hits("The CHLORHEXIDINE MIC is listed in Table S1.", ["chlorhexidine", "Table S1", "PRJNA"])

    assert hits["chlorhexidine"] is True
    assert hits["Table S1"] is True
    assert hits["PRJNA"] is False


def test_build_report_blocks_when_no_candidate_is_ready() -> None:
    candidates = pd.DataFrame(
        [
            {
                "study_key": "proteus_2025_chx_cationic",
                "can_build_without_manual_extraction": False,
                "verification_status": "sequence_public_but_row_level_table_not_found",
            },
            {
                "study_key": "enterococcus_faecium_2019_chx",
                "can_build_without_manual_extraction": False,
                "verification_status": "not_ready_from_automated_public_metadata",
            },
        ]
    )

    report = build_report(candidates)

    assert report["candidate_manifest_extraction_ready"] == 0
    assert report["decision"] == "no_immediate_candidate_manifest_ready_from_automated_metadata"
    assert report["sequence_public_but_row_level_table_not_found"] == 1
