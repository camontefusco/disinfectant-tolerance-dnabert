from pathlib import Path

from disinfectant_tolerance.windows import iter_windows


def test_iter_windows_keeps_large_tail(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fasta"
    fasta.write_text(">contig_a\n" + "ACGT" * 700 + "\n")
    windows = list(
        iter_windows(
            "iso-1",
            fasta,
            label=1,
            group="cluster-a",
            window_size=999,
            step=999,
            minimum_tail_size=500,
        )
    )
    assert [len(window.sequence) for window in windows] == [999, 999, 802]
    assert all(window.isolate_id == "iso-1" for window in windows)

