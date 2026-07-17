#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def interval_distance(left_start: int, left_end: int, right_start: int, right_end: int) -> int:
    if left_end < right_start:
        return right_start - left_end
    if right_end < left_start:
        return left_start - right_end
    return 0


def nearest_evidence(regions: pd.DataFrame, evidence: pd.DataFrame, prefix: str, max_distance: int) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame({"region_id": regions["region_id"]})
    rows = []
    for region in regions.itertuples():
        matches = evidence[
            evidence["assembly_accession"].eq(region.assembly_accession)
            & evidence["contig"].eq(region.contig)
        ].copy()
        if matches.empty:
            rows.append({"region_id": region.region_id})
            continue
        matches["distance"] = [
            interval_distance(region.start, region.end, int(start), int(end))
            for start, end in zip(matches["evidence_start"], matches["evidence_end"])
        ]
        match = matches.sort_values(["distance", "evidence_score"], ascending=[True, False]).iloc[0]
        row = {"region_id": region.region_id}
        if match["distance"] <= max_distance:
            row.update({f"{prefix}_{column}": match[column] for column in matches.columns if column not in {"assembly_accession", "contig"}})
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regions", default="data/manifests/pharma_water_pilot_scored_fragments.csv")
    parser.add_argument("--amrfinder", default="data/manifests/pharma_water_pilot_amrfinderplus_hits.csv")
    parser.add_argument("--bacmet", default="data/manifests/pharma_water_pilot_bacmet_biocide_hits.csv")
    parser.add_argument("--output", default="data/manifests/pharma_water_pilot_evidence_matrix.csv")
    parser.add_argument("--report", default="reports/pharma_water_pilot_evidence_matrix.json")
    parser.add_argument("--max-distance", type=int, default=2000)
    args = parser.parse_args()

    regions = pd.read_csv(args.regions).fillna("")
    amr = pd.read_csv(args.amrfinder).fillna("")
    amr_evidence = pd.DataFrame(
        {
            "assembly_accession": amr["assembly_accession"], "contig": amr["Contig id"],
            "evidence_start": amr["Start"], "evidence_end": amr["Stop"], "evidence_score": amr["% Identity to reference"],
            "element_symbol": amr["Element symbol"], "element_name": amr["Element name"],
            "scope": amr["Scope"], "type": amr["Type"], "class": amr["Class"], "subclass": amr["Subclass"],
        }
    )
    bacmet = pd.read_csv(args.bacmet).fillna("")
    bacmet_evidence = pd.DataFrame(
        {
            "assembly_accession": bacmet["assembly_accession"], "contig": bacmet["contig"],
            "evidence_start": bacmet[["sstart", "send"]].min(axis=1), "evidence_end": bacmet[["sstart", "send"]].max(axis=1),
            "evidence_score": bacmet["identity"], "bacmet_id": bacmet["BacMet_ID"], "gene_name": bacmet["Gene_name"],
            "compound": bacmet["Compound"], "has_qac_annotation": bacmet["has_qac_annotation"],
        }
    ) if not bacmet.empty else pd.DataFrame()
    output = regions.merge(nearest_evidence(regions, amr_evidence, "amrfinder", args.max_distance), on="region_id", how="left")
    output = output.merge(nearest_evidence(regions, bacmet_evidence, "bacmet", args.max_distance), on="region_id", how="left")
    output.to_csv(args.output, index=False)
    report = {
        "claim_boundary": "Exploratory evidence matrix for region prioritization, not tolerance prediction.",
        "regions": len(output),
        "nearby_distance_bp": args.max_distance,
        "regions_with_nearby_amrfinderplus": int(output["amrfinder_element_symbol"].notna().sum()),
        "regions_with_nearby_bacmet_biocide_protein": int(output["bacmet_bacmet_id"].notna().sum()),
        "regions_with_nearby_bacmet_qac_protein": int(
            output.get("bacmet_has_qac_annotation", pd.Series(dtype="string")).astype("string").str.lower().eq("true").sum()
        ),
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
