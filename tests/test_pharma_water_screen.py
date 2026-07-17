from scripts.download_pharma_water_ncbi_pilot import package_url
from scripts.screen_pharma_water_pilot import centered_window


def test_ncbi_package_url_requests_genome_fasta() -> None:
    url = package_url("GCF_000006765.1")
    assert "GCF_000006765.1" in url
    assert "GENOME_FASTA" in url


def test_centered_window_respects_sequence_bounds() -> None:
    sequence = "A" * 1000
    start, end, fragment = centered_window(sequence, center=20, length=500)
    assert (start, end, len(fragment)) == (0, 500, 500)
