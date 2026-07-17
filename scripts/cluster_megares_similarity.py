#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess

import pandas as pd


class UnionFind:
    def __init__(self, values: list[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/manifests/megares_biocide_metal_manifest.csv")
    parser.add_argument("--output", default="data/manifests/megares_similarity_clustered_manifest.csv")
    parser.add_argument("--work-dir", default="data/processed/megares_similarity")
    parser.add_argument("--minimum-identity", type=float, default=90.0)
    parser.add_argument("--minimum-coverage", type=float, default=0.8)
    args = parser.parse_args()

    if shutil.which("makeblastdb") is None or shutil.which("blastn") is None:
        raise RuntimeError("NCBI BLAST+ is required")
    data = pd.read_csv(args.manifest)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta = work_dir / "megares_subset.fasta"
    database = work_dir / "megares_subset"
    hits_path = work_dir / "all_vs_all.tsv"
    with fasta.open("w") as handle:
        for row in data.itertuples():
            handle.write(f">{row.fragment_id}\n{row.sequence}\n")
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(database)], check=True)
    subprocess.run(
        [
            "blastn", "-query", str(fasta), "-db", str(database),
            "-outfmt", "6 qseqid sseqid pident length qlen slen",
            "-perc_identity", str(args.minimum_identity),
            "-qcov_hsp_perc", str(args.minimum_coverage * 100),
            "-max_hsps", "1", "-out", str(hits_path),
        ],
        check=True,
    )
    hits = pd.read_csv(
        hits_path,
        sep="\t",
        names=["query", "subject", "identity", "alignment_length", "query_length", "subject_length"],
    )
    hits["query_coverage"] = hits["alignment_length"] / hits["query_length"]
    hits["subject_coverage"] = hits["alignment_length"] / hits["subject_length"]
    accepted = hits[
        (hits["identity"] >= args.minimum_identity)
        & (hits["query_coverage"] >= args.minimum_coverage)
        & (hits["subject_coverage"] >= args.minimum_coverage)
    ]
    union_find = UnionFind(data["fragment_id"].tolist())
    for row in accepted.itertuples():
        union_find.union(row.query, row.subject)
    roots = {fragment_id: union_find.find(fragment_id) for fragment_id in data["fragment_id"]}
    root_names = {root: f"simcluster_{index:04d}" for index, root in enumerate(sorted(set(roots.values())), start=1)}
    data["similarity_cluster"] = data["fragment_id"].map(lambda value: root_names[roots[value]])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(output, index=False)
    print(f"Wrote {len(data)} fragments across {data['similarity_cluster'].nunique()} similarity clusters to {output}")


if __name__ == "__main__":
    main()
