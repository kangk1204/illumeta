#!/usr/bin/env python3
"""Build non-overwriting GSE66351 AD-vs-CTRL stratified configs.

The source GEO metadata uses generic characteristics columns.  This helper emits
analysis-ready IlluMeta configs using the same explicit schema as the existing
bulk AD-vs-CTRL config so downstream runs do not depend on positional metadata
interpretation.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

OUTPUT_COLUMNS = [
    "primary_group",
    "geo_accession",
    "cell_type",
    "diagnosis",
    "braak_stage",
    "brain_region",
    "age",
    "sex",
    "donor_id",
    "Sentrix_ID",
    "Sentrix_Position",
    "source_name_ch1",
    "description",
]

# configure.tsv column mapping for the GSE66351 source metadata.
FIELD_MAP = {
    "cell_type": "characteristics_ch1",
    "diagnosis": "characteristics_ch1.1",
    "braak_stage": "characteristics_ch1.2",
    "brain_region": "characteristics_ch1.3",
    "age": "characteristics_ch1.4",
    "sex": "characteristics_ch1.5",
    "donor_id": "characteristics_ch1.6",
    "Sentrix_ID": "characteristics_ch1.7",
    "Sentrix_Position": "characteristics_ch1.8",
}


def normalise_row(row: dict[str, str]) -> dict[str, str]:
    out = {"primary_group": row[FIELD_MAP["diagnosis"]], "geo_accession": row["geo_accession"]}
    for target, source in FIELD_MAP.items():
        out[target] = row[source]
    out["source_name_ch1"] = row["source_name_ch1"]
    out["description"] = row["description"]
    return out


def write_config(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="projects/GSE66351/configure.tsv", help="GSE66351 configure.tsv source metadata")
    parser.add_argument("--out-dir", default="projects/GSE66351", help="Directory for generated configs")
    args = parser.parse_args()

    source = Path(args.source)
    out_dir = Path(args.out_dir)
    with source.open(newline="") as handle:
        rows = [normalise_row(r) for r in csv.DictReader(handle, delimiter="\t")]

    neuron = [r for r in rows if r["cell_type"] == "Neuron" and r["brain_region"] == "Occipital cortex"]
    glia = [r for r in rows if r["cell_type"] == "Glia" and r["brain_region"] == "Occipital cortex"]
    frontal_bulk = [r for r in rows if r["cell_type"] == "bulk" and r["brain_region"] == "Frontal cortex"]
    temporal_bulk = [r for r in rows if r["cell_type"] == "bulk" and r["brain_region"] == "Temporal cortex"]

    configs = {
        "configure_AD_vs_CTRL_neuron_occipital.tsv": neuron,
        "configure_AD_vs_CTRL_glia_occipital.tsv": glia,
    }
    for name, config_rows in configs.items():
        write_config(config_rows, out_dir / name)

    summary_rows = []
    for label, config_rows, action in [
        ("Occipital cortex / Neuron", neuron, "config_created_highest_value_run_candidate"),
        ("Occipital cortex / Glia", glia, "config_created_verified_not_run_by_worker2"),
        ("Frontal cortex / bulk", frontal_bulk, "documented_bulk_region_count_not_run"),
        ("Temporal cortex / bulk", temporal_bulk, "documented_bulk_region_count_not_run"),
    ]:
        counts = Counter(r["primary_group"] for r in config_rows)
        summary_rows.append({
            "stratum": label,
            "AD": str(counts.get("AD", 0)),
            "CTRL": str(counts.get("CTRL", 0)),
            "n_total": str(len(config_rows)),
            "action": action,
        })

    summary_path = out_dir / "gse66351_stratified_config_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stratum", "AD", "CTRL", "n_total", "action"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    for row in summary_rows:
        print("\t".join(row[k] for k in ["stratum", "AD", "CTRL", "n_total", "action"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
