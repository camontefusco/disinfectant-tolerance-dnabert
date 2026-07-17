from disinfectant_tolerance.ena import (
    extract_accessions,
    https_urls,
    preferred_accessions,
    sum_semicolon_ints,
)


def test_extract_accessions_handles_parenthetical_alias() -> None:
    assert extract_accessions("ERS4774026 (SAMEA7006197)") == [
        "ERS4774026",
        "SAMEA7006197",
    ]


def test_preferred_accessions_puts_direct_run_first() -> None:
    assert preferred_accessions("SRS123 SRR456") == ["SRR456", "SRS123"]


def test_https_urls_normalizes_ena_paths() -> None:
    assert https_urls("ftp.sra.ebi.ac.uk/a.fastq.gz;ftp://example.org/b.fastq.gz") == [
        "https://ftp.sra.ebi.ac.uk/a.fastq.gz",
        "https://example.org/b.fastq.gz",
    ]


def test_sum_semicolon_ints() -> None:
    assert sum_semicolon_ints("10;20;") == 30


def test_https_urls_handles_single_assembly_path() -> None:
    assert https_urls("ftp.sra.ebi.ac.uk/vol1/analysis/example.fa.gz") == [
        "https://ftp.sra.ebi.ac.uk/vol1/analysis/example.fa.gz"
    ]
