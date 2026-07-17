from __future__ import annotations

import pandas as pd

from scripts.verify_pottier_2023_ddac_manifest import build_decision, parse_assembly_summaries


def test_parse_assembly_summaries_extracts_strain_and_accession() -> None:
    payload = {
        "result": {
            "uids": ["1"],
            "1": {
                "assemblyaccession": "GCF_000000001.1",
                "synonym": {"genbank": "GCA_000000001.1", "refseq": "GCF_000000001.1"},
                "assemblyname": "ASM_TEST",
                "organism": "Pseudomonas aeruginosa",
                "speciesname": "Pseudomonas aeruginosa",
                "biosampleaccn": "SAMN00000001",
                "biosource": {"infraspecieslist": [{"sub_type": "strain", "sub_value": "H-16-00133-1-1"}]},
                "coverage": "85.47",
                "assemblystatus": "Scaffold",
                "contign50": 199916,
                "scaffoldn50": 244776,
                "ftppath_refseq": "ftp://example/refseq",
                "ftppath_genbank": "ftp://example/genbank",
            },
        },
    }

    parsed = parse_assembly_summaries([payload])

    assert parsed.iloc[0]["assembly_accession"] == "GCF_000000001.1"
    assert parsed.iloc[0]["strain"] == "H-16-00133-1-1"
    assert parsed.iloc[0]["biosample_accession"] == "SAMN00000001"


def test_build_decision_blocks_when_ddac_files_are_not_machine_readable() -> None:
    figshare = pd.DataFrame(
        [
            {
                "article_id": 21534780,
                "title": "SUPPLEMENTARY DATA 6. Antimicrobial resistance for the study panel (n=180).",
                "doi": "10.6084/m9.figshare.21534780.v2",
                "file_name": "Supplementary_data6.tiff",
                "file_mimetype": "image/tiff",
                "machine_readable_table_file": False,
                "mentions_ddac": True,
            }
        ]
    )
    assemblies = pd.DataFrame(
        [
            {
                "assembly_accession": "GCF_000000001.1",
                "strain": "H-16-00133-1-1",
            }
        ]
    )

    decision = build_decision(figshare, assemblies)

    assert decision["ncbi_assembly_records"] == 1
    assert decision["machine_readable_ddac_files"] == 0
    assert decision["can_build_supervised_manifest_without_manual_extraction"] is False
    assert decision["decision"] == "blocked_missing_machine_readable_row_level_ddac_labels"
