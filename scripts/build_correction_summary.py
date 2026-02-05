#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
from datetime import datetime


PIPELINES = ("Minfi", "Sesame", "Sesame_Native")
BATCH_FALLBACK_PATTERNS = ("sentrix", "slide", "array", "plate", "chip", "batch", "position")


def safe_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def load_json(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}


def parse_metrics_csv(path):
    metrics = {}
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                key = (row.get("metric") or "").strip()
                if not key:
                    continue
                metrics[key] = (row.get("value") or "").strip()
    except OSError:
        return {}
    return metrics


def parse_pvca(path):
    entries = []
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                term = (row.get("term") or "").strip()
                prop = safe_float(row.get("proportion"))
                if term:
                    entries.append((term, prop))
    except OSError:
        return []
    return entries


def sum_term_proportions(entries, matcher):
    total = 0.0
    matched = False
    for term, prop in entries:
        if prop is None:
            continue
        if matcher(term):
            total += prop
            matched = True
    return total if matched else ""


def batch_matcher(batch_candidate):
    candidate = (batch_candidate or "").strip().lower()
    if candidate:
        return lambda term: candidate in term.lower()
    return lambda term: any(pat in term.lower() for pat in BATCH_FALLBACK_PATTERNS)


def group_matcher():
    return lambda term: "primary_group" in term.lower()


def find_analysis_parameters(analysis_dir):
    primary = os.path.join(analysis_dir, "analysis_parameters.json")
    if os.path.exists(primary):
        return load_json(primary)
    fallback = os.path.join(analysis_dir, "results", "logs", "analysis_parameters.json")
    if os.path.exists(fallback):
        return load_json(fallback)
    return {}


def extract_dataset_id(path):
    matches = re.findall(r"GSE\\d+", path)
    return matches[-1] if matches else ""


def walk_metrics(root_dir, pipelines):
    metrics_files = []
    skip_dirs = {
        "idat", "cache", "results", "logs", "__pycache__", ".git",
        "minfi", "sesame", "sesame_native", "consensus", "reports",
    }
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d not in skip_dirs]
        for pipeline in pipelines:
            name = f"{pipeline}_Metrics.csv"
            if name in filenames:
                metrics_files.append(os.path.join(dirpath, name))
    return metrics_files


def build_rows(root_dir, pipelines):
    rows = []
    for path in walk_metrics(root_dir, pipelines):
        analysis_dir = os.path.dirname(path)
        pipeline = os.path.basename(path).replace("_Metrics.csv", "")
        dataset_id = extract_dataset_id(analysis_dir)
        analysis_name = os.path.basename(analysis_dir)

        metrics = parse_metrics_csv(path)
        params = find_analysis_parameters(analysis_dir)

        pvca_before = parse_pvca(os.path.join(analysis_dir, f"{pipeline}_PVCA.csv"))
        pvca_after = parse_pvca(os.path.join(analysis_dir, f"{pipeline}_AfterCorrection_PVCA.csv"))

        batch_candidate = metrics.get("batch_candidate", "")
        pvca_batch_before = sum_term_proportions(pvca_before, batch_matcher(batch_candidate))
        pvca_batch_after = sum_term_proportions(pvca_after, batch_matcher(batch_candidate))
        pvca_group_before = sum_term_proportions(pvca_before, group_matcher())
        pvca_group_after = sum_term_proportions(pvca_after, group_matcher())

        def delta(before, after):
            if before == "" or after == "":
                return ""
            return before - after

        row = {
            "dataset_id": dataset_id,
            "analysis_name": analysis_name,
            "analysis_dir": analysis_dir,
            "pipeline": pipeline,
            "n_samples": metrics.get("n_samples", ""),
            "n_cpgs": metrics.get("n_cpgs", ""),
            "lambda": metrics.get("lambda", ""),
            "batch_method_applied": metrics.get("batch_method_applied", ""),
            "batch_candidate": batch_candidate,
            "batch_tier": metrics.get("batch_tier", ""),
            "batch_sig_p_lt_0.05_before": metrics.get("batch_sig_p_lt_0.05_before", ""),
            "batch_sig_p_lt_0.05_after": metrics.get("batch_sig_p_lt_0.05_after", ""),
            "batch_min_p_before": metrics.get("batch_min_p_before", ""),
            "batch_min_p_after": metrics.get("batch_min_p_after", ""),
            "group_min_p_before": metrics.get("group_min_p_before", ""),
            "group_min_p_after": metrics.get("group_min_p_after", ""),
            "primary_result_mode": metrics.get("primary_result_mode", ""),
            "tier3_batch": metrics.get("tier3_batch", ""),
            "tier3_primary_lambda": metrics.get("tier3_primary_lambda", ""),
            "perm_ks_p_median": metrics.get("perm_ks_p_median", ""),
            "perm_lambda_median": metrics.get("perm_lambda_median", ""),
            "pvca_batch_before": pvca_batch_before,
            "pvca_batch_after": pvca_batch_after,
            "pvca_batch_delta": delta(pvca_batch_before, pvca_batch_after),
            "pvca_group_before": pvca_group_before,
            "pvca_group_after": pvca_group_after,
            "pvca_group_delta": delta(pvca_group_before, pvca_group_after),
            "array_type": params.get("array_type", ""),
            "tissue": params.get("tissue", ""),
            "preset": params.get("preset", ""),
        }
        rows.append(row)
    return rows


def write_tsv(path, rows):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path, rows, root_dir):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    dataset_count = len({row["dataset_id"] for row in rows if row["dataset_id"]})
    lines = [
        "# IlluMeta Correction Summary",
        "",
        f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- Root: {root_dir}",
        f"- Rows: {len(rows)}",
        f"- Datasets: {dataset_count}",
        f"- Pipelines: {', '.join(PIPELINES)}",
        "",
        "See `correction_summary.tsv` for the full table.",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize IlluMeta batch correction metrics across analyses."
    )
    parser.add_argument("--root", default="projects", help="Root directory to scan (default: projects)")
    parser.add_argument("--output", default="benchmarks/correction_summary.tsv", help="Output TSV path")
    parser.add_argument("--markdown", default="benchmarks/correction_summary.md", help="Output markdown summary")
    args = parser.parse_args()

    rows = build_rows(args.root, PIPELINES)
    rows.sort(key=lambda r: (r["dataset_id"], r["analysis_name"], r["pipeline"]))
    write_tsv(args.output, rows)
    write_markdown(args.markdown, rows, args.root)

    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
