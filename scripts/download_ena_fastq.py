#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shutil
from urllib.parse import urlencode
from urllib.request import urlopen


ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def fastq_urls(accession: str) -> list[str]:
    query = urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp",
            "format": "tsv",
        }
    )
    with urlopen(f"{ENA_API}?{query}") as response:
        rows = list(csv.DictReader(line.decode("utf-8") for line in response))
    urls: list[str] = []
    for row in rows:
        for url in row.get("fastq_ftp", "").split(";"):
            if url:
                urls.append("https://" + url.removeprefix("ftp://"))
    return urls


def download(url: str, output_dir: Path) -> None:
    output = output_dir / Path(url).name
    if output.exists():
        print(f"Exists: {output}")
        return
    print(f"Downloading: {url}")
    with urlopen(url) as response, output.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", default="data/raw/fastq")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with Path(args.manifest).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if args.limit is not None:
        rows = rows[: args.limit]
    for row in rows:
        accession = row.get("ena_accession", "")
        if not accession:
            print(f"Skipping {row.get('isolate_id', '?')}: no ENA accession")
            continue
        urls = fastq_urls(accession)
        if not urls:
            print(f"No FASTQ files found for {accession}")
        for url in urls:
            download(url, output_dir)


if __name__ == "__main__":
    main()

