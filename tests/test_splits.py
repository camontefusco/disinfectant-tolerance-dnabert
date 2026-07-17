from disinfectant_tolerance.splits import assign_grouped_splits


def test_groups_never_cross_splits() -> None:
    rows = [
        {"isolate_id": f"iso-{index}", "group": f"group-{index // 2}", "label": str(index % 2)}
        for index in range(30)
    ]
    assignments = assign_grouped_splits(rows)
    assert set(assignments) == {f"group-{index}" for index in range(15)}
    assert set(assignments.values()) == {"train", "dev", "test"}


def test_randomized_grouped_splits_vary_across_seeds() -> None:
    rows = [
        {"isolate_id": f"iso-{index}", "group": f"group-{index // 2}", "label": str(index % 2)}
        for index in range(60)
    ]
    first = assign_grouped_splits(rows, seed=42, randomization_strength=0.2)
    second = assign_grouped_splits(rows, seed=43, randomization_strength=0.2)
    assert first != second
