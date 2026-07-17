#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator

try:
    from scripts.evaluate_megares_binary_biocide_relevance import select_binary_task, split_data
    from scripts.evaluate_megares_similarity_robustness import make_model
except ModuleNotFoundError:
    from evaluate_megares_binary_biocide_relevance import select_binary_task, split_data
    from evaluate_megares_similarity_robustness import make_model


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


def centered_window(sequence: str, center: int, length: int = 500) -> tuple[int, int, str]:
    start = max(0, min(len(sequence) - length, center - length // 2))
    end = min(len(sequence), start + length)
    return start, end, sequence[start:end]


def deterministic_background(records: dict[str, str], count: int, key: str) -> list[tuple[str, int, int, str]]:
    choices = [(name, sequence) for name, sequence in records.items() if len(sequence) >= 500]
    if not choices:
        return []
    rng = np.random.default_rng(int(hashlib.sha256(key.encode()).hexdigest()[:16], 16))
    output = []
    for _ in range(count):
        name, sequence = choices[int(rng.integers(0, len(choices)))]
        start = int(rng.integers(0, max(1, len(sequence) - 499)))
        output.append((name, start, start + 500, sequence[start : start + 500]))
    return output


def build_megares_fasta(data: pd.DataFrame, path: Path) -> None:
    with path.open("w") as handle:
        for row in data.itertuples():
            handle.write(f">{row.fragment_id}\n{row.sequence}\n")


def blast_hits(assembly_path: Path, accession: str, query: Path, work_dir: Path) -> pd.DataFrame:
    database = work_dir / accession
    hits_path = work_dir / f"{accession}.tsv"
    subprocess.run(["makeblastdb", "-in", str(assembly_path), "-dbtype", "nucl", "-out", str(database)], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        [
            "blastn", "-query", str(query), "-db", str(database), "-task", "megablast",
            "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
            "-perc_identity", "80", "-qcov_hsp_perc", "70", "-max_hsps", "1", "-out", str(hits_path),
        ],
        check=True,
    )
    names = ["fragment_id", "contig", "identity", "alignment_length", "query_length", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    if hits_path.stat().st_size == 0:
        return pd.DataFrame(columns=names)
    return pd.read_csv(hits_path, sep="\t", names=names)


def train_partial_kmer_ensemble(clustered: pd.DataFrame, conditions: Path, seeds: list[int]):
    training = select_binary_task(clustered)
    partial = pd.read_csv(conditions / "partial_500bp.csv").set_index("isolate_id")["sequence"]
    models = []
    for seed in seeds:
        train_indices, dev_indices, _ = split_data(training, seed)
        train, dev = training.iloc[train_indices], training.iloc[dev_indices]
        model = make_model("kmer_random_forest", seed)
        model.fit(partial.loc[train["fragment_id"]], train["label"])
        calibrated = CalibratedClassifierCV(FrozenEstimator(model), method="sigmoid")
        calibrated.fit(partial.loc[dev["fragment_id"]], dev["label"])
        models.append(calibrated)
    return models


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", default="data/manifests/pharma_water_downloaded_assemblies.csv")
    parser.add_argument("--megares", default="data/manifests/megares_similarity_clustered_manifest.csv")
    parser.add_argument("--conditions-dir", default="data/processed/megares_conditions")
    parser.add_argument("--work-dir", default="data/processed/pharma_water_pilot_screen")
    parser.add_argument("--fragments", default="data/manifests/pharma_water_pilot_scored_fragments.csv")
    parser.add_argument("--hits", default="data/manifests/pharma_water_pilot_megares_hits.csv")
    parser.add_argument("--report", default="reports/pharma_water_pilot_screen.json")
    parser.add_argument("--background-per-assembly", type=int, default=20)
    args = parser.parse_args()

    if shutil.which("makeblastdb") is None or shutil.which("blastn") is None:
        raise RuntimeError("NCBI BLAST+ is required")
    assemblies = pd.read_csv(args.assemblies)
    megares = pd.read_csv(args.megares)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    query = work_dir / "megares_subset.fasta"
    build_megares_fasta(megares, query)
    metadata = megares.set_index("fragment_id")[["type", "class", "mechanism", "group"]]
    all_hits, rows = [], []
    for assembly in assemblies.itertuples():
        records = read_fasta(Path(assembly.assembly_path))
        hits = blast_hits(Path(assembly.assembly_path), assembly.assembly_accession, query, work_dir)
        hits["assembly_accession"] = assembly.assembly_accession
        hits["taxon_query"] = assembly.taxon_query
        all_hits.append(hits)
        for hit in hits.itertuples():
            sequence = records[hit.contig]
            center = (min(hit.sstart, hit.send) + max(hit.sstart, hit.send)) // 2
            start, end, fragment = centered_window(sequence, center - 1)
            rows.append(
                {
                    "assembly_accession": assembly.assembly_accession, "taxon_query": assembly.taxon_query,
                    "contig": hit.contig, "start": start, "end": end, "source": "megares_hit",
                    "best_megares_fragment_id": hit.fragment_id, "best_megares_identity": hit.identity,
                    "best_megares_class": metadata.loc[hit.fragment_id, "class"], "sequence": fragment,
                }
            )
        for contig, start, end, fragment in deterministic_background(records, args.background_per_assembly, assembly.assembly_accession):
            rows.append(
                {
                    "assembly_accession": assembly.assembly_accession, "taxon_query": assembly.taxon_query,
                    "contig": contig, "start": start, "end": end, "source": "background_sample",
                    "best_megares_fragment_id": "", "best_megares_identity": "", "best_megares_class": "", "sequence": fragment,
                }
            )
        print(f"{assembly.assembly_accession}: {len(hits)} MEGARes hits and {args.background_per_assembly} background windows")
    hit_frame = pd.concat(all_hits, ignore_index=True).merge(metadata, on="fragment_id", how="left")
    hit_output = Path(args.hits)
    hit_output.parent.mkdir(parents=True, exist_ok=True)
    hit_frame.to_csv(hit_output, index=False)
    fragments = pd.DataFrame(rows).drop_duplicates(["assembly_accession", "contig", "start", "end", "source"]).reset_index(drop=True)
    models = train_partial_kmer_ensemble(megares, Path(args.conditions_dir), [42, 43, 44, 45, 46])
    probabilities = np.column_stack([model.predict_proba(fragments["sequence"])[:, 1] for model in models])
    fragments.insert(0, "region_id", [f"PWR{i:06d}" for i in range(1, len(fragments) + 1)])
    fragments["exploratory_biocide_relevance_score"] = probabilities.mean(axis=1)
    fragments["score_std_across_splits"] = probabilities.std(axis=1)
    fragments["bacmet_layer_status"] = "protein_reference_annotation_only"
    fragments["amrfinderplus_status"] = "not_run"
    fragments.drop(columns=["sequence"]).to_csv(args.fragments, index=False)
    report = {
        "claim_boundary": "Exploratory fragment prioritization, not isolate-level sanitizer-survival or tolerance prediction.",
        "assemblies": len(assemblies),
        "megares_hits": len(hit_frame),
        "scored_regions": len(fragments),
        "scored_regions_by_source": fragments["source"].value_counts().to_dict(),
        "scores": {
            "model": "Five-seed calibrated 500 bp k-mer random-forest ensemble",
            "median": float(fragments["exploratory_biocide_relevance_score"].median()),
            "maximum": float(fragments["exploratory_biocide_relevance_score"].max()),
        },
        "interpretability": {
            "megares": "BLAST nucleotide matches included per candidate region.",
            "bacmet": "Protein-based BacMet evidence remains an annotation reference; direct nucleotide screening was not claimed.",
            "amrfinderplus": "Not run yet; add in the next annotation-hardening step.",
        },
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
