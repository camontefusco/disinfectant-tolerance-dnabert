from dataclasses import dataclass

import torch

from scripts.embed_dnabert2 import hidden_state


@dataclass
class ModelOutput:
    last_hidden_state: torch.Tensor


def test_hidden_state_accepts_transformers_output() -> None:
    expected = torch.ones((1, 2, 3))
    assert hidden_state(ModelOutput(expected)) is expected


def test_hidden_state_accepts_dnabert2_tuple() -> None:
    expected = torch.ones((1, 2, 3))
    assert hidden_state((expected,)) is expected
