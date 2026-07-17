from scripts.audit_megares_v3_extension import TARGET_TYPES
from scripts.train_megares_homology_aware_baseline import normalized_target


def test_megares_subset_includes_biocides_metals_and_multi_compound() -> None:
    assert TARGET_TYPES == {"Biocides", "Metals", "Multi-compound"}


def test_multi_compound_target_is_preserved() -> None:
    assert normalized_target("Multi-compound") == "Multi-compound"
    assert normalized_target("Biocides") == "Biocides"
