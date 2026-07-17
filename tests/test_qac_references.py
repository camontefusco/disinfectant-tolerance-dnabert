from scripts.fetch_qac_references import REFERENCES


def test_qac_reference_panel_contains_expected_determinants() -> None:
    assert {reference[0] for reference in REFERENCES} == {"bcrABC", "qacH", "emrE", "emrC"}
