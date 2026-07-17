from scripts.audit_bacmet_experimental import BIOCIDE_PATTERN, QAC_PATTERN


def test_bacmet_compound_patterns_find_qac_and_biocide_examples() -> None:
    compound = "Benzylkonium Chloride (BAC) [class: Quaternary Ammonium Compounds (QACs)]"
    assert QAC_PATTERN.search(compound)
    assert BIOCIDE_PATTERN.search(compound)
