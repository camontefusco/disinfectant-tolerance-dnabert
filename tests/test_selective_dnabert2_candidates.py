from scripts.prepare_selective_dnabert2_candidates import reverse_complement


def test_reverse_complement_preserves_unknown_base() -> None:
    assert reverse_complement("AAGTN") == "NACTT"
