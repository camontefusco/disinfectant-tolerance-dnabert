import pandas as pd

from scripts.build_publication_dossier import safe_ratio, source_enrichment


def test_safe_ratio_handles_zero_denominator() -> None:
    assert safe_ratio(1, 0) == 0.0


def test_source_enrichment_compares_hits_and_background() -> None:
    frame = pd.DataFrame({"source": ["megares_hit", "megares_hit", "background_sample"]})
    flag = pd.Series([True, False, False])
    result = source_enrichment(frame, flag)
    assert result["hit_rate"] == 0.5
    assert result["background_rate"] == 0.0
