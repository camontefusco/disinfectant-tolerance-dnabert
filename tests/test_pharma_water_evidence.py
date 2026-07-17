from scripts.build_pharma_water_evidence_matrix import interval_distance
from scripts.download_pharma_water_ncbi_annotations import package_url


def test_annotation_package_requests_proteins_and_gff() -> None:
    url = package_url("GCF_000006765.1")
    assert "PROT_FASTA" in url
    assert "GENOME_GFF" in url


def test_interval_distance_handles_overlap_and_gap() -> None:
    assert interval_distance(10, 20, 15, 30) == 0
    assert interval_distance(10, 20, 30, 40) == 10
