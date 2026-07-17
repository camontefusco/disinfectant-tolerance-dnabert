from __future__ import annotations

import csv
import json
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import urlopen


ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_SEARCH_API = "https://www.ebi.ac.uk/ena/portal/api/search"
ACCESSION_PATTERN = re.compile(r"\b[A-Z]{2,6}\d+(?:\.\d+)?\b")
RUN_PREFIXES = ("ERR", "SRR", "DRR")


def extract_accessions(value: str) -> list[str]:
    """Extract ordered, unique archive identifiers from a supplementary-table cell."""
    accessions: list[str] = []
    for accession in ACCESSION_PATTERN.findall(value.upper()):
        if accession not in accessions:
            accessions.append(accession)
    return accessions


def preferred_accessions(value: str) -> list[str]:
    """Try direct run accessions before sample or experiment identifiers."""
    accessions = extract_accessions(value)
    return sorted(accessions, key=lambda accession: accession[:3] not in RUN_PREFIXES)


def https_urls(value: str) -> list[str]:
    urls: list[str] = []
    for url in value.split(";"):
        url = url.strip()
        if not url:
            continue
        if url.startswith("ftp://"):
            url = "https://" + url.removeprefix("ftp://")
        elif not url.startswith(("http://", "https://")):
            url = "https://" + url
        urls.append(url)
    return urls


def sum_semicolon_ints(value: str) -> int:
    return sum(int(part) for part in value.split(";") if part.strip().isdigit())


def query_read_runs(accession: str, timeout: int = 60) -> list[dict[str, str]]:
    fields = (
        "run_accession,sample_accession,secondary_sample_accession,"
        "study_accession,scientific_name,fastq_ftp,fastq_bytes"
    )
    query = urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": fields,
            "format": "tsv",
        }
    )
    with urlopen(f"{ENA_API}?{query}", timeout=timeout) as response:
        return list(
            csv.DictReader(
                (line.decode("utf-8") for line in response),
                delimiter="\t",
            )
        )


def query_assembly_analyses(
    sample_accession: str,
    timeout: int = 60,
) -> list[dict[str, str]]:
    fields = (
        "analysis_accession,analysis_type,sample_accession,submitted_ftp,"
        "submitted_bytes,generated_ftp,generated_bytes,assembly_type,"
        "assembly_quality,assembly_software"
    )
    query = urlencode(
        {
            "result": "analysis",
            "query": f'sample_accession="{sample_accession}"',
            "fields": fields,
            "format": "tsv",
        }
    )
    with urlopen(f"{ENA_SEARCH_API}?{query}", timeout=timeout) as response:
        rows = list(
            csv.DictReader(
                (line.decode("utf-8") for line in response),
                delimiter="\t",
            )
        )
    return [row for row in rows if row.get("analysis_type") == "SEQUENCE_ASSEMBLY"]


def resolve_accession_cell(
    value: str,
    cache: dict[str, list[dict[str, str]]] | None = None,
) -> list[dict[str, str]]:
    cache = cache if cache is not None else {}
    resolved: dict[str, dict[str, str]] = {}
    for accession in preferred_accessions(value):
        if accession not in cache:
            cache[accession] = query_read_runs(accession)
        for row in cache[accession]:
            run_accession = row.get("run_accession", "")
            if run_accession:
                resolved[run_accession] = row
        if resolved:
            break
    return list(resolved.values())


def load_cache(path: str | Path) -> dict[str, list[dict[str, str]]]:
    path = Path(path)
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def save_cache(path: str | Path, cache: dict[str, list[dict[str, str]]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=2, sort_keys=True))
