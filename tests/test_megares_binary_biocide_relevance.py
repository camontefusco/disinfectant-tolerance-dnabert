import pandas as pd

from scripts.evaluate_megares_binary_biocide_relevance import select_binary_task


def test_binary_task_selects_biocide_relevant_and_metal_only_fragments() -> None:
    data = pd.DataFrame(
        [
            {"type": "Biocides", "class": "Biocide resistance"},
            {"type": "Multi-compound", "class": "Drug and biocide resistance"},
            {"type": "Metals", "class": "Metal resistance"},
            {"type": "Multi-compound", "class": "Drug and metal resistance"},
        ]
    )
    selected = select_binary_task(data)
    assert selected["target"].tolist() == ["biocide_relevant", "biocide_relevant", "metal_only"]
