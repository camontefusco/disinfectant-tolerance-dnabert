from pathlib import Path

from scripts.screen_qac_determinants import best_hits


def test_qac_screen_returns_expected_columns_when_no_hits(tmp_path: Path) -> None:
    references = tmp_path / "references.fasta"
    references.write_text(">qacH\nAAAAAA\n")
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig\nCCCCCC\n")
    frame = best_hits(references, assembly)
    assert frame.empty
    assert {"gene", "coverage"}.issubset(frame.columns)
