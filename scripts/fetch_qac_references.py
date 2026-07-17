#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


REFERENCES = [
    ("bcrABC", "JX023284.1", None, None, 1),
    ("qacH", "MK944277.1", None, None, 1),
    ("emrE", "CP001602.2", 1850347, 1850670, 2),
    ("emrC", "LT732640.1", 1575, 1961, 1),
]


def fetch_reference(
    gene: str,
    accession: str,
    start: int | None,
    stop: int | None,
    strand: int,
) -> str:
    params: dict[str, str | int] = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    if start is not None and stop is not None:
        params.update({"seq_start": start, "seq_stop": stop, "strand": strand})
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(params)
    text = urlopen(url, timeout=60).read().decode()
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"NCBI did not return FASTA for {gene}: {accession}")
    sequence = "".join(lines[1:]).upper()
    if not sequence:
        raise RuntimeError(f"NCBI returned an empty sequence for {gene}: {accession}")
    return f">{gene}|{accession}\n{sequence}\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/references/qac_determinants.fasta")
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = [fetch_reference(*reference) for reference in REFERENCES]
    output.write_text("".join(records))
    print(f"Wrote {len(records)} QAC determinant references to {output}")


if __name__ == "__main__":
    main()
