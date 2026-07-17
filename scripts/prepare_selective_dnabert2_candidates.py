#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                records[name] = "".join(parts).upper()
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name:
        records[name] = "".join(parts).upper()
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regions", default="data/manifests/pharma_water_pilot_scored_fragments.csv")
    parser.add_argument("--hits", default="data/manifests/pharma_water_pilot_megares_hits.csv")
    parser.add_argument("--assemblies", default="data/manifests/pharma_water_downloaded_assemblies.csv")
    parser.add_argument("--per-taxon", type=int, default=25)
    parser.add_argument("--output", default="data/processed/pharma_water_selective_dnabert2_candidates.csv")
    parser.add_argument("--manifest", default="data/manifests/pharma_water_selective_dnabert2_candidates.csv")
    args = parser.parse_args()

    regions = pd.read_csv(args.regions).fillna("")
    selected = (
        regions[regions["source"].eq("megares_hit")]
        .sort_values("exploratory_biocide_relevance_score", ascending=False)
        .groupby("taxon_query", sort=False)
        .head(args.per_taxon)
        .copy()
    )
    hits = pd.read_csv(args.hits)
    assemblies = pd.read_csv(args.assemblies).set_index("assembly_accession")
    rows = []
    fasta_cache: dict[str, dict[str, str]] = {}
    for region in selected.itertuples():
        matches = hits[
            hits["assembly_accession"].eq(region.assembly_accession)
            & hits["contig"].eq(region.contig)
            & hits["fragment_id"].eq(region.best_megares_fragment_id)
        ].sort_values("bitscore", ascending=False)
        if matches.empty:
            raise ValueError(f"Could not resolve intact MEGARes hit for {region.region_id}")
        hit = matches.iloc[0]
        if region.assembly_accession not in fasta_cache:
            fasta_cache[region.assembly_accession] = read_fasta(Path(assemblies.loc[region.assembly_accession, "assembly_path"]))
        contig = fasta_cache[region.assembly_accession][region.contig]
        start, end = min(int(hit["sstart"]), int(hit["send"])) - 1, max(int(hit["sstart"]), int(hit["send"]))
        sequence = contig[start:end]
        strand = "+" if int(hit["sstart"]) <= int(hit["send"]) else "-"
        if strand == "-":
            sequence = reverse_complement(sequence)
        rows.append(
            {
                "isolate_id": region.region_id, "sequence": sequence, "region_id": region.region_id,
                "taxon_query": region.taxon_query, "assembly_accession": region.assembly_accession,
                "contig": region.contig, "intact_start": start, "intact_end": end, "strand": strand,
                "aligned_length": len(sequence), "best_megares_fragment_id": region.best_megares_fragment_id,
                "partial_kmer_score": region.exploratory_biocide_relevance_score,
            }
        )
    output = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output[["isolate_id", "sequence"]].to_csv(args.output, index=False)
    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    output.drop(columns=["sequence"]).to_csv(args.manifest, index=False)
    print(f"Wrote {len(output)} selective intact candidate regions")


if __name__ == "__main__":
    main()
