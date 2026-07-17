from scripts.evaluate_kmer_repeated_grouped_splits import make_model


def test_repeated_kmer_model_uses_requested_ngram_range() -> None:
    model = make_model(4, 6, seed=42)
    assert model["tfidf"].ngram_range == (4, 6)
    assert model["classifier"].random_state == 42
