from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .io import read_fasta


@dataclass(frozen=True)
class Window:
    isolate_id: str
    contig_id: str
    start: int
    end: int
    sequence: str
    label: int
    group: str


def iter_windows(
    isolate_id: str,
    assembly_path: str | Path,
    label: int,
    group: str,
    window_size: int = 999,
    step: int = 999,
    minimum_tail_size: int = 500,
) -> Iterator[Window]:
    if window_size <= 0 or step <= 0:
        raise ValueError("window_size and step must be positive")
    for contig_id, sequence in read_fasta(assembly_path):
        if len(sequence) < minimum_tail_size:
            continue
        starts = list(range(0, max(1, len(sequence) - window_size + 1), step))
        if not starts:
            starts = [0]
        last_end = starts[-1] + window_size
        if len(sequence) - last_end >= minimum_tail_size:
            starts.append(last_end)
        for start in starts:
            fragment = sequence[start : start + window_size]
            if len(fragment) < minimum_tail_size:
                continue
            yield Window(
                isolate_id=isolate_id,
                contig_id=contig_id,
                start=start,
                end=start + len(fragment),
                sequence=fragment,
                label=label,
                group=group,
            )

