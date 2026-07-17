from scripts.prepare_megares_fragment_conditions import centered_fragment, noisy_sequence


def test_centered_fragment_uses_middle_of_long_sequence() -> None:
    assert centered_fragment("AAAACCCC", 4) == "AACC"


def test_noisy_sequence_is_deterministic_and_changes_sequence() -> None:
    first = noisy_sequence("A" * 100, 0.01, "fragment")
    second = noisy_sequence("A" * 100, 0.01, "fragment")
    assert first == second
    assert first != "A" * 100
