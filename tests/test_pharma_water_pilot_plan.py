import pandas as pd

from scripts.plan_pharma_water_ncbi_pilot import dataset_report_url, select_pilot


def test_dataset_report_url_encodes_taxon() -> None:
    assert "Pseudomonas%20aeruginosa" in dataset_report_url("Pseudomonas aeruginosa", 5)


def test_select_pilot_prefers_complete_genome() -> None:
    rows = pd.DataFrame(
        [
            {"taxon_query": "x", "assembly_accession": "GCF_2", "assembly_level": "Contig", "checkm_completeness": 100, "checkm_contamination": 0},
            {"taxon_query": "x", "assembly_accession": "GCF_1", "assembly_level": "Complete Genome", "checkm_completeness": 90, "checkm_contamination": 1},
        ]
    )
    assert select_pilot(rows, 1)["assembly_accession"].tolist() == ["GCF_1"]
