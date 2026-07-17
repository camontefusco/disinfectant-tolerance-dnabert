from scripts.audit_harrand_2026_feasibility import normalize


def test_harrand_overlap_normalizer_ignores_punctuation() -> None:
    assert normalize("FSL C9-0001") == normalize("FSL_C9_0001")
