from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import random


@dataclass(frozen=True)
class SplitFractions:
    train: float = 0.7
    dev: float = 0.15
    test: float = 0.15

    def validate(self) -> None:
        if abs(self.train + self.dev + self.test - 1.0) > 1e-9:
            raise ValueError("Split fractions must sum to 1")
        if min(self.train, self.dev, self.test) <= 0:
            raise ValueError("Split fractions must be positive")


def assign_grouped_splits(
    rows: list[dict[str, str]],
    fractions: SplitFractions = SplitFractions(),
    seed: int = 42,
    randomization_strength: float = 0.0,
) -> dict[str, str]:
    """Assign entire lineage groups to one split with approximate class balance."""
    fractions.validate()
    if randomization_strength < 0:
        raise ValueError("randomization_strength must be non-negative")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["group"]].append(row)

    rng = random.Random(seed)
    groups = list(grouped.items())
    rng.shuffle(groups)
    groups.sort(key=lambda item: len(item[1]), reverse=True)

    split_names = ("train", "dev", "test")
    targets = {
        "train": fractions.train,
        "dev": fractions.dev,
        "test": fractions.test,
    }
    totals = Counter(row["label"] for row in rows)
    assigned_counts = {name: Counter() for name in split_names}
    assigned_sizes = Counter()
    assignments: dict[str, str] = {}

    def cost(split: str, group_rows: list[dict[str, str]]) -> float:
        proposed_sizes = assigned_sizes.copy()
        proposed_sizes[split] += len(group_rows)
        proposed_counts = {
            name: assigned_counts[name].copy()
            for name in split_names
        }
        proposed_counts[split].update(row["label"] for row in group_rows)

        size_cost = sum(
            abs(proposed_sizes[name] - len(rows) * targets[name])
            / max(1, len(rows) * targets[name])
            for name in split_names
        )
        class_cost = sum(
            abs(proposed_counts[name][label] - totals[label] * targets[name])
            / max(1, totals[label] * targets[name])
            for name in split_names
            for label in totals
        )
        return size_cost + class_cost

    for group, group_rows in groups:
        split = min(
            split_names,
            key=lambda name: cost(name, group_rows)
            + rng.uniform(0.0, randomization_strength),
        )
        assignments[group] = split
        assigned_sizes[split] += len(group_rows)
        assigned_counts[split].update(row["label"] for row in group_rows)
    return assignments
