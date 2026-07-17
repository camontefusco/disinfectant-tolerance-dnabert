from scripts.build_sampled_windows import evenly_spaced


def test_evenly_spaced_keeps_boundaries() -> None:
    assert evenly_spaced(list(range(10)), 4) == [0, 3, 6, 9]


def test_evenly_spaced_does_not_pad_short_input() -> None:
    assert evenly_spaced([1, 2], 4) == [1, 2]

