#!/usr/bin/env python3
import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime


DEFAULT_VARIANTS = {
    "baseline": [],
    "no_sva": ["--disable-sva"],
    "no_auto_covariates": ["--disable-auto-covariates"],
    "no_sva_no_auto_covariates": ["--disable-sva", "--disable-auto-covariates"],
    "batch_none": ["--batch-method", "none"],
    "batch_combat": ["--batch-method", "combat"],
    "batch_limma": ["--batch-method", "limma"],
    "batch_sva": ["--batch-method", "sva"],
    "no_cross_reactive": ["--unsafe-skip-cross-reactive"],
    "aggressive_preset": ["--preset", "aggressive"],
    "minfi_only": ["--skip-sesame"],
}
EXTRA_VARIANTS = {
    # Optional variants (not included in the default run list).
    "sesame_typeinorm": ["--sesame-typeinorm"],
}

PIPELINES = ("Minfi", "Sesame", "Sesame_Native")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ablation variants of IlluMeta analysis and collect summary metrics."
    )
    parser.add_argument("-i", "--input-dir", type=str, help="Project directory (containing configure.tsv)")
    parser.add_argument("-c", "--config", type=str, help="Path to configure.tsv (overrides --input-dir)")
    parser.add_argument("--group-con", required=True, help="Control group label (case-insensitive)")
    parser.add_argument("--group-test", required=True, help="Test group label (case-insensitive)")
    parser.add_argument(
        "--out-root",
        type=str,
        default=None,
        help="Root directory for ablation outputs (default: ablation_runs/<timestamp>)",
    )
    parser.add_argument(
        "--variants",
        type=str,
        default=",".join(DEFAULT_VARIANTS.keys()),
        help="Comma-separated variant names to run (default: all built-ins). "
             f"Optional: {', '.join(EXTRA_VARIANTS.keys())}",
    )
    parser.add_argument(
        "--extra-args",
        type=str,
        default="",
        help="Extra analysis args appended to every variant (e.g., '--tissue Blood')",
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Skip a variant if its output already contains summary.json",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running analysis",
    )
    return parser.parse_args()


def ensure_out_root(root):
    if root:
        return root
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join("ablation_runs", ts)


def parse_metrics_csv(path):
    metrics = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (row.get("metric") or "").strip()
            if not key:
                continue
            metrics[key] = (row.get("value") or "").strip()
    return metrics


def load_json(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None


def to_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def write_tsv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_variant(cmd, dry_run):
    if dry_run:
        print("DRY RUN:", " ".join(shlex.quote(c) for c in cmd))
        return 0
    return subprocess.call(cmd)


def main():
    args = parse_args()
    if not args.input_dir and not args.config:
        print("Error: provide --input-dir or --config", file=sys.stderr)
        return 2

    out_root = ensure_out_root(args.out_root)
    os.makedirs(out_root, exist_ok=True)
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]

    all_variants = {}
    all_variants.update(DEFAULT_VARIANTS)
    all_variants.update(EXTRA_VARIANTS)
    invalid = [v for v in variants if v not in all_variants]
    if invalid:
        print(f"Error: unknown variants: {', '.join(invalid)}", file=sys.stderr)
        return 2

    base_cmd = [sys.executable, "illumeta.py", "analysis"]
    if args.config:
        base_cmd += ["--config", args.config]
    else:
        base_cmd += ["--input-dir", args.input_dir]
    base_cmd += ["--group_con", args.group_con, "--group_test", args.group_test]
    extra_args = shlex.split(args.extra_args) if args.extra_args else []

    manifest_rows = []
    metrics_rows = []
    counts_rows = []
    params_by_variant = {}

    for variant in variants:
        variant_dir = os.path.join(out_root, variant)
        os.makedirs(variant_dir, exist_ok=True)
        summary_path = os.path.join(variant_dir, "summary.json")
        if args.reuse_existing and os.path.exists(summary_path):
            status = "reused"
        else:
            cmd = base_cmd + ["--output", variant_dir] + all_variants[variant] + extra_args
            status = "ran"
            exit_code = run_variant(cmd, args.dry_run)
            if exit_code != 0 and not args.dry_run:
                status = f"failed({exit_code})"

        manifest_rows.append(
            {
                "variant": variant,
                "output_dir": variant_dir,
                "status": status,
            }
        )

        summary = load_json(summary_path)
        if summary:
            counts_rows.append({"variant": variant, **summary})

        params = load_json(os.path.join(variant_dir, "analysis_parameters.json"))
        if params:
            params_by_variant[variant] = params

        for pipeline in PIPELINES:
            metrics_path = os.path.join(variant_dir, f"{pipeline}_Metrics.csv")
            if not os.path.exists(metrics_path):
                continue
            metrics = parse_metrics_csv(metrics_path)
            for metric, value in metrics.items():
                metrics_rows.append(
                    {
                        "variant": variant,
                        "pipeline": pipeline,
                        "metric": metric,
                        "value": value,
                    }
                )

    write_tsv(os.path.join(out_root, "ablation_manifest.tsv"), manifest_rows, ["variant", "output_dir", "status"])
    if counts_rows:
        count_fields = sorted({k for row in counts_rows for k in row.keys()})
        if "variant" in count_fields:
            count_fields.remove("variant")
        write_tsv(
            os.path.join(out_root, "ablation_counts.tsv"),
            counts_rows,
            ["variant"] + count_fields,
        )
    if metrics_rows:
        write_tsv(
            os.path.join(out_root, "ablation_metrics_long.tsv"),
            metrics_rows,
            ["variant", "pipeline", "metric", "value"],
        )
        # Delta vs baseline
        baseline_metrics = {}
        for row in metrics_rows:
            if row["variant"] != "baseline":
                continue
            key = (row["pipeline"], row["metric"])
            baseline_metrics[key] = to_float(row.get("value"))
        delta_rows = []
        for row in metrics_rows:
            if row["variant"] == "baseline":
                continue
            key = (row["pipeline"], row["metric"])
            base_val = baseline_metrics.get(key)
            curr_val = to_float(row.get("value"))
            delta = ""
            if base_val is not None and curr_val is not None:
                delta = curr_val - base_val
            delta_rows.append(
                {
                    "variant": row["variant"],
                    "pipeline": row["pipeline"],
                    "metric": row["metric"],
                    "value": row.get("value", ""),
                    "baseline_value": "" if base_val is None else base_val,
                    "delta": delta,
                }
            )
        if delta_rows:
            write_tsv(
                os.path.join(out_root, "ablation_metrics_delta.tsv"),
                delta_rows,
                ["variant", "pipeline", "metric", "value", "baseline_value", "delta"],
            )
    if params_by_variant:
        with open(os.path.join(out_root, "ablation_parameters.json"), "w", encoding="utf-8") as handle:
            json.dump(params_by_variant, handle, indent=2, ensure_ascii=True)

    if counts_rows:
        baseline = next((row for row in counts_rows if row.get("variant") == "baseline"), None)
        if baseline:
            delta_rows = []
            baseline_vals = {k: to_float(v) for k, v in baseline.items() if k != "variant"}
            for row in counts_rows:
                if row.get("variant") == "baseline":
                    continue
                out = {"variant": row.get("variant", "")}
                for key, base_val in baseline_vals.items():
                    curr_val = to_float(row.get(key))
                    if base_val is None or curr_val is None:
                        continue
                    out[f"{key}_delta"] = curr_val - base_val
                if len(out) > 1:
                    delta_rows.append(out)
            if delta_rows:
                fieldnames = sorted({k for r in delta_rows for k in r.keys() if k != "variant"})
                write_tsv(
                    os.path.join(out_root, "ablation_counts_delta.tsv"),
                    delta_rows,
                    ["variant"] + fieldnames,
                )

    print(f"[*] Ablation results saved under: {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
