#!/usr/bin/env python3
import argparse
import csv
import json
import os
from datetime import datetime


PIPELINES = ("Minfi", "Sesame", "Sesame_Native")
FIELDNAMES = [
    "gse_id",
    "pipeline",
    "result_dir",
    "group_counts",
    "n_con",
    "n_test",
    "n_samples",
    "n_cpgs",
    "qc_samples_total",
    "qc_samples_passed",
    "qc_pass_rate",
    "qc_probes_raw",
    "qc_probes_final",
    "qc_probe_retention",
    "lambda",
    "lambda_guard_status",
    "lambda_guard_action",
    "lambda_guard_threshold",
    "lambda_guard_lambda",
    "batch_method_applied",
    "batch_candidate",
    "batch_tier",
    "batch_sig_p_lt_0.05_before",
    "batch_sig_p_lt_0.05_after",
    "batch_sig_reduction",
    "batch_sig_reduction_frac",
    "batch_min_p_before",
    "batch_min_p_after",
    "group_min_p_before",
    "group_min_p_after",
    "perm_ks_p_median",
    "perm_lambda_median",
    "perm_mean_sig",
    "perm_max_sig",
    "vp_primary_group_mean",
    "dmp_sig",
    "dmr_count",
    "pvca_batch_before",
    "pvca_batch_after",
    "pvca_batch_delta",
    "pvca_group_before",
    "pvca_group_after",
    "pvca_group_delta",
    "n_covariates_used",
    "n_sv_used",
    "tier3_batch",
    "tier3_primary_lambda",
    "primary_tier3_meta_method",
    "primary_tier3_meta_i2_median",
    "primary_branch",
    "primary_result_mode",
    "primary_dmr_status",
    "primary_dmr_reason",
    "primary_branch_override",
    "primary_branch_reason",
    "sesame_dyebias_mode",
    "sesame_dyebias_note",
    "sesame_typeinorm_enabled",
    "intersection_logFC_correlation",
    "intersection_jaccard_overlap",
    "intersection_n_both",
    "intersection_native_logFC_correlation",
    "intersection_native_jaccard_overlap",
    "intersection_native_n_both",
    "array_type",
    "preset",
    "tissue",
    "tissue_source",
]


def read_tsv(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    return rows


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def write_tsv(path, rows, fieldnames):
    _ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None


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


def parse_kv_csv(path):
    rows = {}
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                key = (row.get("metric") or "").strip()
                value = (row.get("value") or "").strip()
                if key:
                    rows[key] = value
    except OSError:
        return {}
    return rows


def _to_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def _ratio(numer, denom):
    n = _to_float(numer)
    d = _to_float(denom)
    if n is None or d in (None, 0):
        return ""
    return str(n / d)


def parse_group_counts(path):
    if not os.path.exists(path):
        return ""
    out = []
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                group = (row.get("group") or "").strip()
                count = (row.get("n") or "").strip()
                if group and count:
                    out.append(f"{group}:{count}")
    except OSError:
        return ""
    return ";".join(out)


def parse_pvca(path):
    entries = []
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                term = (row.get("term") or "").strip()
                prop = _to_float(row.get("proportion"))
                if term:
                    entries.append((term, prop))
    except OSError:
        return []
    return entries


def _sum_terms(entries, matcher):
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
    fallback = ("sentrix", "slide", "array", "plate", "chip", "batch", "position")
    return lambda term: any(pat in term.lower() for pat in fallback)


def group_matcher():
    return lambda term: "primary_group" in term.lower()


def count_rows(path):
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            return max(0, sum(1 for _ in handle) - 1)
    except OSError:
        return ""


def parse_decision_ledger(path):
    primary_branch = ""
    tier3_batch = ""
    if not os.path.exists(path):
        return primary_branch, tier3_batch
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                stage = (row.get("stage") or "").strip()
                decision = (row.get("decision") or "").strip()
                value = (row.get("value") or "").strip()
                if stage == "consensus" and decision == "primary_branch":
                    primary_branch = value
                if stage == "confounding" and decision == "tier3_batch":
                    tier3_batch = value
    except OSError:
        return primary_branch, tier3_batch
    return primary_branch, tier3_batch


def resolve_analysis_dir(row, projects_root):
    analysis_dir = (row.get("analysis_dir") or "").strip()
    if analysis_dir:
        candidate = analysis_dir
        if not os.path.isabs(candidate):
            candidate = os.path.abspath(candidate)
            if not os.path.exists(candidate) and projects_root:
                alt = os.path.join(projects_root, analysis_dir)
                if os.path.exists(alt):
                    candidate = alt
        if os.path.exists(os.path.join(candidate, "summary.json")):
            return candidate
    return None


def find_latest_result_dir(projects_root, gse_id):
    base_dir = os.path.join(projects_root, gse_id)
    if not os.path.isdir(base_dir):
        return None
    candidates = []
    summary_path = os.path.join(base_dir, "summary.json")
    if os.path.exists(summary_path):
        candidates.append((summary_path, base_dir))
    for name in os.listdir(base_dir):
        path = os.path.join(base_dir, name)
        if not os.path.isdir(path):
            continue
        summary = os.path.join(path, "summary.json")
        if os.path.exists(summary):
            candidates.append((summary, path))
    if not candidates:
        return None
    candidates.sort(key=lambda item: os.path.getmtime(item[0]), reverse=True)
    return candidates[0][1]


def build_rows(rows, projects_root):
    out_rows = []
    missing = []
    for row in rows:
        gse_id = (row.get("gse_id") or "").strip()
        if not gse_id:
            continue
        result_dir = resolve_analysis_dir(row, projects_root)
        if not result_dir:
            result_dir = find_latest_result_dir(projects_root, gse_id)
        if not result_dir:
            missing.append(gse_id)
            continue

        summary = load_json(os.path.join(result_dir, "summary.json")) or {}
        params = load_json(os.path.join(result_dir, "analysis_parameters.json")) or {}
        qc_summary = parse_kv_csv(os.path.join(result_dir, "QC_Summary.csv"))
        group_counts = parse_group_counts(os.path.join(result_dir, "Input_Group_Distribution.csv"))
        primary_branch, tier3_batch = parse_decision_ledger(os.path.join(result_dir, "decision_ledger.tsv"))
        intersection_metrics = parse_metrics_csv(os.path.join(result_dir, "Intersection_Comparison_Metrics.csv"))
        intersection_native_metrics = parse_metrics_csv(os.path.join(result_dir, "Intersection_Native_Comparison_Metrics.csv"))
        qc_samples_total = qc_summary.get("Total_samples_input", "")
        qc_samples_passed = qc_summary.get("Samples_passed_QC", "")
        qc_pass_rate = _ratio(qc_samples_passed, qc_samples_total)
        qc_probes_raw = qc_summary.get("Total_probes_raw", "")
        qc_probes_final = qc_summary.get("Probes_final", "")
        qc_probe_retention = _ratio(qc_probes_final, qc_probes_raw)

        for pipeline in PIPELINES:
            metrics_path = os.path.join(result_dir, f"{pipeline}_Metrics.csv")
            if not os.path.exists(metrics_path):
                continue
            metrics = parse_metrics_csv(metrics_path)
            batch_sig_before = metrics.get("batch_sig_p_lt_0.05_before", "")
            batch_sig_after = metrics.get("batch_sig_p_lt_0.05_after", "")
            batch_sig_reduction = ""
            b_before = _to_float(batch_sig_before)
            b_after = _to_float(batch_sig_after)
            if b_before is not None and b_after is not None:
                batch_sig_reduction = str(b_before - b_after)
            batch_sig_reduction_frac = ""
            if b_before not in (None, 0) and b_after is not None:
                batch_sig_reduction_frac = str((b_before - b_after) / b_before)

            pvca_before = parse_pvca(os.path.join(result_dir, f"{pipeline}_PVCA.csv"))
            pvca_after = parse_pvca(os.path.join(result_dir, f"{pipeline}_AfterCorrection_PVCA.csv"))
            batch_candidate = metrics.get("batch_candidate", "")
            pvca_batch_before = _sum_terms(pvca_before, batch_matcher(batch_candidate))
            pvca_batch_after = _sum_terms(pvca_after, batch_matcher(batch_candidate))
            pvca_group_before = _sum_terms(pvca_before, group_matcher())
            pvca_group_after = _sum_terms(pvca_after, group_matcher())
            pvca_batch_delta = "" if pvca_batch_before == "" or pvca_batch_after == "" else pvca_batch_before - pvca_batch_after
            pvca_group_delta = "" if pvca_group_before == "" or pvca_group_after == "" else pvca_group_before - pvca_group_after

            dmp_sig = ""
            if pipeline == "Minfi":
                dmp_sig = _to_float(summary.get("minfi_up")) or 0
                dmp_sig += _to_float(summary.get("minfi_down")) or 0
            elif pipeline == "Sesame":
                dmp_sig = _to_float(summary.get("sesame_up")) or 0
                dmp_sig += _to_float(summary.get("sesame_down")) or 0
            elif pipeline == "Sesame_Native":
                dmp_sig = _to_float(summary.get("sesame_native_up")) or 0
                dmp_sig += _to_float(summary.get("sesame_native_down")) or 0

            dmr_count = count_rows(os.path.join(result_dir, f"{pipeline}_DMRs.csv"))

            out_rows.append(
                {
                    "gse_id": gse_id,
                    "pipeline": pipeline,
                    "result_dir": result_dir,
                    "group_counts": group_counts,
                    "n_con": summary.get("n_con", ""),
                    "n_test": summary.get("n_test", ""),
                    "n_samples": metrics.get("n_samples", ""),
                    "n_cpgs": metrics.get("n_cpgs", ""),
                    "qc_samples_total": qc_samples_total,
                    "qc_samples_passed": qc_samples_passed,
                    "qc_pass_rate": qc_pass_rate,
                    "qc_probes_raw": qc_probes_raw,
                    "qc_probes_final": qc_probes_final,
                    "qc_probe_retention": qc_probe_retention,
                    "lambda": metrics.get("lambda", ""),
                    "lambda_guard_status": metrics.get("lambda_guard_status", ""),
                    "lambda_guard_action": metrics.get("lambda_guard_action", ""),
                    "lambda_guard_threshold": metrics.get("lambda_guard_threshold", ""),
                    "lambda_guard_lambda": metrics.get("lambda_guard_lambda", ""),
                    "batch_method_applied": metrics.get("batch_method_applied", ""),
                    "batch_candidate": batch_candidate,
                    "batch_tier": metrics.get("batch_tier", ""),
                    "batch_sig_p_lt_0.05_before": batch_sig_before,
                    "batch_sig_p_lt_0.05_after": batch_sig_after,
                    "batch_sig_reduction": batch_sig_reduction,
                    "batch_sig_reduction_frac": batch_sig_reduction_frac,
                    "batch_min_p_before": metrics.get("batch_min_p_before", ""),
                    "batch_min_p_after": metrics.get("batch_min_p_after", ""),
                    "group_min_p_before": metrics.get("group_min_p_before", ""),
                    "group_min_p_after": metrics.get("group_min_p_after", ""),
                    "perm_ks_p_median": metrics.get("perm_ks_p_median", ""),
                    "perm_lambda_median": metrics.get("perm_lambda_median", ""),
                    "perm_mean_sig": metrics.get("perm_mean_sig", ""),
                    "perm_max_sig": metrics.get("perm_max_sig", ""),
                    "vp_primary_group_mean": metrics.get("vp_primary_group_mean", ""),
                    "dmp_sig": dmp_sig,
                    "dmr_count": dmr_count,
                    "pvca_batch_before": pvca_batch_before,
                    "pvca_batch_after": pvca_batch_after,
                    "pvca_batch_delta": pvca_batch_delta,
                    "pvca_group_before": pvca_group_before,
                    "pvca_group_after": pvca_group_after,
                    "pvca_group_delta": pvca_group_delta,
                    "n_covariates_used": metrics.get("n_covariates_used", ""),
                    "n_sv_used": metrics.get("n_sv_used", ""),
                    "primary_branch": primary_branch,
                    "tier3_batch": tier3_batch,
                    "tier3_primary_lambda": metrics.get("tier3_primary_lambda", ""),
                    "primary_tier3_meta_method": summary.get("primary_tier3_meta_method", ""),
                    "primary_tier3_meta_i2_median": summary.get("primary_tier3_meta_i2_median", ""),
                    "primary_result_mode": summary.get("primary_result_mode", ""),
                    "primary_dmr_status": summary.get("primary_dmr_status", ""),
                    "primary_dmr_reason": summary.get("primary_dmr_reason", ""),
                    "primary_branch_override": summary.get("primary_branch_override", ""),
                    "primary_branch_reason": summary.get("primary_branch_reason", ""),
                    "sesame_dyebias_mode": summary.get("sesame_dyebias_mode", ""),
                    "sesame_dyebias_note": summary.get("sesame_dyebias_note", ""),
                    "sesame_typeinorm_enabled": params.get("sesame_typeinorm_enabled", ""),
                    "intersection_logFC_correlation": intersection_metrics.get("logFC_correlation", ""),
                    "intersection_jaccard_overlap": intersection_metrics.get("jaccard_overlap", ""),
                    "intersection_n_both": intersection_metrics.get("n_both", ""),
                    "intersection_native_logFC_correlation": intersection_native_metrics.get("logFC_correlation", ""),
                    "intersection_native_jaccard_overlap": intersection_native_metrics.get("jaccard_overlap", ""),
                    "intersection_native_n_both": intersection_native_metrics.get("n_both", ""),
                    "array_type": params.get("array_type", ""),
                    "preset": params.get("preset", ""),
                    "tissue": params.get("tissue", ""),
                    "tissue_source": params.get("tissue_source", ""),
                }
            )
    return out_rows, missing


def write_markdown(path, input_tsv, projects_root, rows, missing):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    gse_ids = sorted({row["gse_id"] for row in rows})
    pipeline_counts = {p: 0 for p in PIPELINES}
    for row in rows:
        pipeline_counts[row["pipeline"]] = pipeline_counts.get(row["pipeline"], 0) + 1

    lines = [
        "# IlluMeta Benchmark Summary",
        "",
        f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- Input TSV: {input_tsv}",
        f"- Projects root: {projects_root}",
        f"- Datasets with results: {len(gse_ids)}",
        f"- Rows written: {len(rows)}",
        "- Pipelines: " + ", ".join(f"{p}={pipeline_counts[p]}" for p in PIPELINES),
    ]
    if missing:
        lines.append(f"- Missing results: {len(missing)} (see console output)")
    lines.append("")
    lines.append("See `benchmark_summary.tsv` for full metrics.")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build a long-form benchmark table from IlluMeta result folders."
    )
    parser.add_argument(
        "--input-tsv",
        default="geo_idat_methylation.tsv",
        help="Input TSV with a gse_id column",
    )
    parser.add_argument(
        "--projects-root",
        default="projects",
        help="Projects root containing GSE output folders",
    )
    parser.add_argument(
        "--output-tsv",
        default=os.path.join("benchmarks", "benchmark_summary.tsv"),
        help="Output TSV path",
    )
    parser.add_argument(
        "--output-md",
        default=os.path.join("benchmarks", "benchmark_summary.md"),
        help="Output Markdown summary path",
    )
    args = parser.parse_args()

    rows = read_tsv(args.input_tsv)
    out_rows, missing = build_rows(rows, args.projects_root)
    if not out_rows:
        print("No benchmark rows found. Check input TSV and projects root.", flush=True)
        return 1

    fieldnames = FIELDNAMES
    write_tsv(args.output_tsv, out_rows, fieldnames)
    write_markdown(args.output_md, args.input_tsv, args.projects_root, out_rows, missing)

    print(f"Wrote {len(out_rows)} rows to {args.output_tsv}")
    if missing:
        print(f"Missing results for {len(missing)} datasets: {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
