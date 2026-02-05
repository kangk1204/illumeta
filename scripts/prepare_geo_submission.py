#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import shutil


def describe_file(path: str) -> str:
    name = os.path.basename(path)
    if name.endswith("_BetaMatrix_PreFilter.tsv.gz"):
        return "Pre-filter beta matrix (after normalization, before probe QC)"
    if name.endswith("_MvalueMatrix_PreFilter.tsv.gz"):
        return "Pre-filter M-value matrix (after normalization, before probe QC)"
    if name.endswith("_MethylatedSignal_PreFilter.tsv.gz"):
        return "Pre-filter methylated signal matrix (after normalization, before probe QC)"
    if name.endswith("_UnmethylatedSignal_PreFilter.tsv.gz"):
        return "Pre-filter unmethylated signal matrix (after normalization, before probe QC)"
    if name.endswith("_DetectionP_PreFilter.tsv.gz"):
        return "Pre-filter detection P-value matrix (all samples)"
    if name.endswith("_BetaMatrix.tsv.gz"):
        return "Processed beta matrix used for modeling"
    if name.endswith("_MvalueMatrix.tsv.gz"):
        return "Processed M-value matrix (logit of beta) used for statistics"
    if name.endswith("_DMPs_full.csv"):
        return "Full DMP results table"
    if name.endswith("_DMRs.csv"):
        return "DMR results table"
    if name == "QC_Summary.csv":
        return "Sample/probe QC summary"
    if name == "summary.json":
        return "Run summary (primary branch/mode/DMR status)"
    if name == "methods.md":
        return "Auto-generated methods summary"
    if name == "analysis_parameters.json":
        return "Run parameters and thresholds"
    return "Supplementary output"


def collect_files(result_dir: str):
    patterns = [
        "*_BetaMatrix_PreFilter.tsv.gz",
        "*_MvalueMatrix_PreFilter.tsv.gz",
        "*_MethylatedSignal_PreFilter.tsv.gz",
        "*_UnmethylatedSignal_PreFilter.tsv.gz",
        "*_DetectionP_PreFilter.tsv.gz",
        "*_BetaMatrix.tsv.gz",
        "*_MvalueMatrix.tsv.gz",
        "*_DMPs_full.csv",
        "*_DMRs.csv",
        "QC_Summary.csv",
        "summary.json",
        "methods.md",
        "analysis_parameters.json",
    ]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(result_dir, pat)))
    return sorted(set(files))


def write_manifest(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "size_bytes", "description"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Prepare a GEO submission bundle from an IlluMeta result directory.")
    parser.add_argument("--result-dir", required=True, help="IlluMeta analysis output directory")
    parser.add_argument("--out-dir", default=None, help="Output directory for submission bundle (default: <result-dir>/geo_submission)")
    parser.add_argument("--copy", action="store_true", help="Copy files into out-dir (otherwise manifest only)")
    args = parser.parse_args()

    result_dir = os.path.abspath(args.result_dir)
    out_dir = os.path.abspath(args.out_dir) if args.out_dir else os.path.join(result_dir, "geo_submission")
    os.makedirs(out_dir, exist_ok=True)

    files = collect_files(result_dir)
    if not files:
        print("No GEO submission files found. Ensure analysis has completed.")
        return 1

    rows = []
    for path in files:
        size = os.path.getsize(path)
        desc = describe_file(path)
        rel = os.path.basename(path)
        rows.append({"file": rel, "size_bytes": str(size), "description": desc})
        if args.copy:
            shutil.copy2(path, os.path.join(out_dir, rel))

    manifest_path = os.path.join(out_dir, "geo_submission_manifest.tsv")
    write_manifest(rows, manifest_path)
    print(f"Wrote manifest: {manifest_path}")
    if args.copy:
        print(f"Copied {len(files)} files to: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
