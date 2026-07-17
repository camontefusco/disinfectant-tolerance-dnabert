from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator


REQUIRED_MANIFEST_COLUMNS = {"isolate_id", "assembly_path", "label", "group"}


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Manifest is empty: {path}")
    missing = REQUIRED_MANIFEST_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
    for row in rows:
        if row["label"] not in {"0", "1"}:
            raise ValueError(f"Expected binary label for {row['isolate_id']}")
        if not row["group"]:
            raise ValueError(f"Missing group for {row['isolate_id']}")
    return rows


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    parts: list[str] = []
    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()

