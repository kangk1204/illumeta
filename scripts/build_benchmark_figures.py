#!/usr/bin/env python3
import argparse
import csv
import math
import os


def read_tsv(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


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


def write_md(path, rows, fieldnames):
    _ensure_parent_dir(path)
    lines = []
    lines.append("| " + " | ".join(fieldnames) + " |")
    lines.append("| " + " | ".join(["---"] * len(fieldnames)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(f, "")) for f in fieldnames) + " |")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def to_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def select_primary(rows, primary_only):
    if not primary_only:
        return rows
    grouped = {}
    for row in rows:
        gse_id = (row.get("gse_id") or "").strip()
        if not gse_id:
            continue
        grouped.setdefault(gse_id, []).append(row)
    selected = []
    for gse_id, items in grouped.items():
        primary = None
        for row in items:
            if row.get("pipeline") == row.get("primary_branch"):
                primary = row
                break
        if primary is None:
            for row in items:
                if row.get("pipeline") == "Minfi":
                    primary = row
                    break
        if primary is None:
            primary = items[0]
        selected.append(primary)
    return selected


def build_primary_table(rows):
    fields = [
        "gse_id",
        "pipeline",
        "consensus_primary_view",
        "primary_result_mode",
        "primary_tier3_meta_method",
        "primary_tier3_meta_i2_median",
        "qc_pass_rate",
        "qc_probe_retention",
        "batch_sig_reduction",
        "perm_lambda_median",
        "lambda",
        "vp_primary_group_mean",
        "intersection_native_n_both",
        "intersection_n_both",
        "intersection_native_jaccard_overlap",
        "intersection_jaccard_overlap",
        "primary_dmr_status",
        "primary_branch_override",
        "primary_branch_reason",
    ]
    out_rows = []
    for row in rows:
        record = {field: row.get(field, "") for field in fields}
        # IlluMeta dashboard/reporting treats native intersection as primary, strict as sensitivity.
        record["consensus_primary_view"] = "Intersection_Native"
        out_rows.append(record)
    return out_rows, fields


def plot_summary(rows, out_path):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure generation.")
        return False

    labels = [row.get("gse_id", "") for row in rows]
    n = len(labels)
    if n == 0:
        print("No rows available for plotting.")
        return False

    def get_series(key, default=0.0):
        vals = []
        for row in rows:
            val = to_float(row.get(key))
            if val is None or math.isnan(val):
                vals.append(default)
            else:
                vals.append(val)
        return vals

    qc_pass = get_series("qc_pass_rate")
    probe_ret = get_series("qc_probe_retention")
    batch_red = get_series("batch_sig_reduction")
    perm_lambda = get_series("perm_lambda_median", default=1.0)
    vp_primary = get_series("vp_primary_group_mean")
    jaccard_strict = get_series("intersection_jaccard_overlap")
    jaccard_native = get_series("intersection_native_jaccard_overlap")
    pvca_batch_delta = get_series("pvca_batch_delta")

    fig_w = max(10, n * 0.6)
    fig, axes = plt.subplots(3, 2, figsize=(fig_w, 12))
    x = list(range(n))

    width = 0.35
    ax = axes[0, 0]
    ax.bar([i - width / 2 for i in x], qc_pass, width, label="QC pass rate")
    ax.bar([i + width / 2 for i in x], probe_ret, width, label="Probe retention")
    ax.set_ylim(0, 1)
    ax.set_title("QC Pass Rate & Probe Retention")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.bar(x, batch_red, color="#4c78a8")
    ax.set_title("Batch Signal Reduction (Before - After)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")

    ax = axes[1, 0]
    ax.bar(x, perm_lambda, color="#f58518")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ax.set_title("Permutation Lambda (Median)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")

    ax = axes[1, 1]
    ax.bar(x, vp_primary, color="#54a24b")
    ax.set_title("VariancePartition: primary_group mean")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")

    ax = axes[2, 0]
    ax.bar([i - width / 2 for i in x], jaccard_native, width, color="#b279a2", label="Native (primary)")
    ax.bar([i + width / 2 for i in x], jaccard_strict, width, color="#9d9da1", label="Strict (sensitivity)")
    ax.set_title("Intersection Jaccard (Native vs Strict)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend(fontsize=8)

    ax = axes[2, 1]
    ax.bar(x, pvca_batch_delta, color="#72b7b2")
    ax.set_title("PVCA Batch Delta (Before - After)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")

    fig.tight_layout()
    _ensure_parent_dir(out_path)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Build paper-ready summary table and figure from benchmark_summary.tsv."
    )
    parser.add_argument("--input-tsv", required=True, help="Input benchmark_summary.tsv")
    parser.add_argument("--out-dir", default="benchmarks", help="Output directory")
    parser.add_argument("--all-rows", action="store_true", default=False, help="Use all rows (do not collapse to primary)")
    args = parser.parse_args()

    rows = read_tsv(args.input_tsv)
    rows = select_primary(rows, primary_only=not args.all_rows)
    rows.sort(key=lambda r: r.get("gse_id", ""))

    out_rows, fields = build_primary_table(rows)
    out_tsv = os.path.join(args.out_dir, "benchmark_primary_summary.tsv")
    out_md = os.path.join(args.out_dir, "benchmark_primary_summary.md")
    write_tsv(out_tsv, out_rows, fields)
    write_md(out_md, out_rows, fields)

    fig_path = os.path.join(args.out_dir, "benchmark_primary_summary.png")
    plot_summary(rows, fig_path)

    print(f"Wrote {len(out_rows)} rows to {out_tsv}")
    print(f"Wrote Markdown table to {out_md}")
    print(f"Wrote figure to {fig_path}")


if __name__ == "__main__":
    main()
