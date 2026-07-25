#!/usr/bin/env python3
"""Build the deterministic, synthetic dashboard used for interface QA."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import illumeta  # noqa: E402


DEMO_DIR = ROOT / "docs" / "dashboard_demo"
RESULTS_DIR = DEMO_DIR / "Case_vs_Control_results"
HTML_PATH = DEMO_DIR / "Case_vs_Control_results_index.html"
MANIFEST_PATH = DEMO_DIR / "manifest.json"
SOURCE_DATE_EPOCH = int(os.environ.get("SOURCE_DATE_EPOCH", "1784764800"))

SUMMARY = {
    "n_con": 6,
    "n_test": 6,
    "minfi_up": 4,
    "minfi_down": 3,
    "sesame_up": 3,
    "sesame_down": 3,
    "sesame_native_up": 3,
    "sesame_native_down": 3,
    "intersect_up": 2,
    "intersect_down": 1,
    "intersect_native_up": 3,
    "intersect_native_down": 2,
    "primary_branch": "Minfi",
    "primary_lambda_guard_status": "ok",
}

QC_ROWS = [
    ("Total_samples_input", 12),
    ("Samples_failed_QC", 0),
    ("Samples_passed_QC", 12),
    ("Samples_failed_sex_mismatch", 0),
    ("Sex_mismatch_samples", 0),
    ("Total_probes_raw", 10000),
    ("Probes_failed_detection", 120),
    ("Probes_cross_reactive", 350),
    ("Probes_with_SNPs", 500),
    ("Probes_sex_chromosomes", 200),
    ("Probes_final", 8830),
]

PARAMETERS = {
    "pval_threshold": 0.05,
    "lfc_threshold": 0.5,
    "delta_beta_threshold": 0.0,
    "tissue": "Synthetic dashboard QA fixture",
    "tissue_source": "fixture",
    "array_type": "450K",
    "crf_sample_tier": "small",
    "cell_reference": "not_applicable",
    "cell_reference_platform": "not_applicable",
    "auto_covariates_enabled": False,
    "auto_covariates_exclude_group_associated": True,
    "auto_covariates_group_assoc_p_threshold": 0.05,
    "auto_covariates_max_cor": 0.7,
    "auto_covariate_alpha": 0.01,
    "auto_covariate_max_pcs": 5,
    "sesame_native_na_max_frac": 0.1,
    "overcorrection_guard_ratio": 0.8,
    "undercorrection_guard_min_sv": 1,
    "sva_enabled": True,
    "sva_inclusion_rule": "best_method_or_no_batch",
    "clock_covariates_enabled": False,
    "sv_group_p_threshold": 0.05,
    "sv_group_eta2_threshold": 0.14,
    "pvca_min_samples": 10,
    "dmr_maxgap": 500,
    "dmr_p_cutoff": 0.05,
    "unsafe_skip_cross_reactive": False,
    "cross_reactive_enabled": True,
    "cross_reactive_active": True,
    "cross_reactive_count": 350,
    "sex_check_action": "flag",
    "sex_mismatch_count": 0,
}

METRICS_ROWS = [
    ("lambda", 1.03),
    ("batch_method_applied", "SVA"),
    ("n_samples", 12),
    ("n_cpgs", 8830),
    ("n_covariates_used", 2),
    ("covariates_used", "age;sex"),
    ("n_sv_used", 2),
    ("batch_sig_p_lt_0.05_before", 3),
    ("batch_sig_p_lt_0.05_after", 0),
    ("batch_min_p_before", 0.0001),
    ("batch_min_p_after", 0.35),
    ("group_min_p_before", 0.01),
    ("group_min_p_after", 0.002),
    ("perm_mean_sig", 1),
    ("perm_max_sig", 2),
    ("lambda_ratio", 1.4),
]

CRF_FIXTURES = {
    "CRF_MMC_Summary.csv": (
        (
            "top_k",
            "core_pct",
            "spearman_mean",
            "direction_mean",
            "mmc_composite",
            "methods",
        ),
        (5, 0.68, 0.82, 0.91, 0.79, 2),
    ),
    "CRF_NCS_Summary.csv": (
        (
            "stage",
            "type",
            "n",
            "lambda",
            "sig_rate",
            "lambda_ci_low",
            "lambda_ci_high",
            "ncs_score",
            "ncs_score_tol",
        ),
        ("corrected", "snp", 100, 1.03, 0.048, 0.98, 1.08, 0.92, 0.1),
    ),
    "CRF_RSS_Summary.csv": (
        (
            "top_k",
            "overlap_mean_corr",
            "rss_mean_corr",
            "jaccard_mean_corr",
            "rbo_mean_corr",
            "rbo_p",
            "sign_mean_corr",
        ),
        (5, 0.74, 0.65, 0.59, 0.71, 0.9, 0.93),
    ),
    "CRF_Sample_Tier.csv": (
        ("tier", "min_per_group"),
        ("small", 3),
    ),
}

DMP_COLUMNS = ("CpG", "Gene", "logFC", "Delta_Beta", "adj.P.Val", "Direction")
DMP_ROWS = [
    ("cg90000001", "SYNTH_A", 0.82, 0.12, 0.0004, "Hyper"),
    ("cg90000002", "SYNTH_B", -0.75, -0.11, 0.0012, "Hypo"),
    ("cg90000003", "SYNTH_C", 0.61, 0.08, 0.0048, "Hyper"),
    ("cg90000004", "SYNTH_D", -0.58, -0.07, 0.0091, "Hypo"),
    ("cg90000005", "SYNTH_E", 0.54, 0.06, 0.018, "Hyper"),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interactive_plot_html() -> str:
    points = [
        {
            "cpg": row[0],
            "gene": row[1],
            "logfc": row[2],
            "score": -math.log10(row[4]),
        }
        for row in DMP_ROWS
    ]
    template = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <link rel="icon" href="data:,">
  <title>Synthetic Minfi volcano interaction fixture</title>
  <style>
    :root { color-scheme: light; font-family: Arial, sans-serif; }
    body { margin: 0; color: #17202a; background: #f7f9fb; }
    main { max-width: 900px; margin: 0 auto; padding: 24px; }
    h1 { margin: 0 0 4px; font-size: 24px; }
    p { color: #52606d; }
    .toolbar { display: flex; gap: 8px; margin: 16px 0; }
    button { border: 1px solid #9aa7b2; background: white; padding: 8px 12px; cursor: pointer; }
    button:focus-visible, circle:focus-visible { outline: 3px solid #0b6e99; outline-offset: 2px; }
    .plot-wrap { position: relative; border: 1px solid #cbd5df; background: white; }
    svg { display: block; width: 100%; height: auto; min-height: 360px; }
    .axis { stroke: #607080; stroke-width: 1.2; }
    .guide { stroke: #b5c0ca; stroke-dasharray: 5 5; }
    .point { cursor: crosshair; stroke: white; stroke-width: 2; }
    .point.hyper { fill: #b42318; }
    .point.hypo { fill: #1261a0; }
    #tooltip { position: absolute; display: none; pointer-events: none; padding: 7px 9px;
      color: white; background: #17202a; font-size: 13px; border-radius: 4px; }
    #status { min-height: 20px; color: #243b53; font-weight: 600; }
  </style>
</head>
<body data-illumeta-fixture="interactive-plot">
<main>
  <h1>Synthetic Minfi volcano plot</h1>
  <p><strong>Synthetic interface fixture.</strong> Hover or focus a point to inspect it. These are not study results.</p>
  <div class="toolbar" aria-label="Plot controls">
    <button type="button" id="zoomIn" title="Zoom in">+ Zoom</button>
    <button type="button" id="reset" title="Reset plot view">Reset</button>
  </div>
  <div class="plot-wrap">
    <svg id="plot" viewBox="0 0 720 420" role="img" aria-labelledby="plotTitle plotDesc">
      <title id="plotTitle">Synthetic volcano plot</title>
      <desc id="plotDesc">Five focusable CpG points plotted by log fold change and negative log ten adjusted P value.</desc>
      <line class="axis" x1="70" y1="350" x2="690" y2="350"></line>
      <line class="axis" x1="70" y1="30" x2="70" y2="350"></line>
      <line class="guide" x1="380" y1="30" x2="380" y2="350"></line>
      <text x="340" y="395">logFC</text>
      <text transform="translate(20 250) rotate(-90)">-log10 adjusted P</text>
    </svg>
    <div id="tooltip" role="tooltip"></div>
  </div>
  <p id="status" aria-live="polite">Focus or hover a point.</p>
</main>
<script>
  const points = __POINTS__;
  const svg = document.getElementById("plot");
  const tooltip = document.getElementById("tooltip");
  const status = document.getElementById("status");
  const ns = "http://www.w3.org/2000/svg";
  const x = value => 70 + ((value + 1) / 2) * 620;
  const y = value => 350 - (value / 4) * 300;

  function describe(point) {
    return `${point.cpg} ${point.gene}: logFC ${point.logfc}, adjusted P score ${point.score.toFixed(2)}`;
  }
  function showPoint(point, circle) {
    const message = describe(point);
    const wrap = svg.parentElement.getBoundingClientRect();
    const box = circle.getBoundingClientRect();
    tooltip.textContent = message;
    tooltip.style.left = `${box.left - wrap.left + 10}px`;
    tooltip.style.top = `${box.top - wrap.top - 34}px`;
    tooltip.style.display = "block";
    status.textContent = message;
  }
  function hidePoint() {
    tooltip.style.display = "none";
  }
  for (const point of points) {
    const circle = document.createElementNS(ns, "circle");
    circle.setAttribute("class", `point ${point.logfc >= 0 ? "hyper" : "hypo"}`);
    circle.setAttribute("cx", x(point.logfc));
    circle.setAttribute("cy", y(point.score));
    circle.setAttribute("r", "9");
    circle.setAttribute("tabindex", "0");
    circle.setAttribute("aria-label", describe(point));
    circle.addEventListener("pointerenter", () => showPoint(point, circle));
    circle.addEventListener("pointerleave", hidePoint);
    circle.addEventListener("focus", () => showPoint(point, circle));
    circle.addEventListener("blur", hidePoint);
    svg.appendChild(circle);
  }
  document.getElementById("zoomIn").addEventListener("click", () => {
    svg.setAttribute("viewBox", "100 20 560 330");
    status.textContent = "Plot zoomed in.";
  });
  document.getElementById("reset").addEventListener("click", () => {
    svg.setAttribute("viewBox", "0 0 720 420");
    status.textContent = "Plot view reset.";
  });
</script>
</body>
</html>
"""
    return template.replace(
        "__POINTS__",
        json.dumps(points, separators=(",", ":"), sort_keys=True),
    )


def interactive_table_html(title: str) -> str:
    payload = {
        "columns": DMP_COLUMNS,
        "rows": DMP_ROWS,
    }
    template = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <link rel="icon" href="data:,">
  <title>__TITLE__</title>
  <style>
    :root { color-scheme: light; font-family: Arial, sans-serif; }
    body { margin: 0; color: #17202a; background: #f7f9fb; }
    main { max-width: 1100px; margin: 0 auto; padding: 24px; }
    h1 { margin: 0 0 4px; font-size: 24px; }
    p { color: #52606d; }
    .toolbar { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin: 16px 0; }
    input { min-width: 240px; padding: 8px; border: 1px solid #9aa7b2; }
    button { border: 1px solid #9aa7b2; background: white; padding: 8px 12px; cursor: pointer; }
    button:focus-visible, input:focus-visible { outline: 3px solid #0b6e99; outline-offset: 2px; }
    .table-wrap { overflow-x: auto; border: 1px solid #cbd5df; background: white; }
    table { border-collapse: collapse; width: 100%; min-width: 720px; }
    th, td { border-bottom: 1px solid #d9e0e7; padding: 9px 10px; text-align: left; }
    th { background: #edf2f7; }
    th button { border: 0; background: transparent; padding: 0; font-weight: 700; }
    tbody tr:hover { background: #eef7fb; }
    #status { min-height: 20px; color: #243b53; font-weight: 600; }
  </style>
</head>
<body data-illumeta-fixture="interactive-table">
<main>
  <h1>__TITLE__</h1>
  <p><strong>Synthetic interface fixture.</strong> Search, sort, and export these rows. They are not study results.</p>
  <div class="toolbar">
    <label for="search">Search</label>
    <input id="search" type="search" autocomplete="off" placeholder="CpG, gene, direction">
    <button type="button" id="downloadCsv">Export CSV</button>
  </div>
  <div class="table-wrap">
    <table>
      <thead><tr id="header"></tr></thead>
      <tbody id="body"></tbody>
    </table>
  </div>
  <p id="status" aria-live="polite"></p>
</main>
<script>
  const source = __PAYLOAD__;
  let rows = source.rows.slice();
  let sortIndex = -1;
  let ascending = true;
  const header = document.getElementById("header");
  const body = document.getElementById("body");
  const search = document.getElementById("search");
  const status = document.getElementById("status");

  source.columns.forEach((column, index) => {
    const th = document.createElement("th");
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = `${column} (sort)`;
    button.dataset.index = index;
    button.addEventListener("click", () => sortRows(index));
    th.appendChild(button);
    header.appendChild(th);
  });

  function filteredRows() {
    const query = search.value.trim().toLowerCase();
    if (!query) return rows;
    return rows.filter(row => row.some(value => String(value).toLowerCase().includes(query)));
  }
  function render() {
    const visible = filteredRows();
    body.replaceChildren();
    for (const row of visible) {
      const tr = document.createElement("tr");
      for (const value of row) {
        const td = document.createElement("td");
        td.textContent = value;
        tr.appendChild(td);
      }
      body.appendChild(tr);
    }
    status.textContent = `${visible.length} of ${rows.length} rows shown.`;
  }
  function sortRows(index) {
    if (sortIndex === index) ascending = !ascending;
    else {
      sortIndex = index;
      ascending = true;
    }
    rows.sort((left, right) => {
      const a = left[index];
      const b = right[index];
      const cmp = typeof a === "number" && typeof b === "number"
        ? a - b
        : String(a).localeCompare(String(b));
      return ascending ? cmp : -cmp;
    });
    render();
  }
  function csvCell(value) {
    return `"${String(value).replaceAll('"', '""')}"`;
  }
  function downloadCsv() {
    const csv = [source.columns, ...filteredRows()]
      .map(row => row.map(csvCell).join(","))
      .join("\\n") + "\\n";
    const url = URL.createObjectURL(new Blob([csv], {type: "text/csv;charset=utf-8"}));
    const link = document.createElement("a");
    link.href = url;
    link.download = "illumeta_synthetic_dmps.csv";
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 0);
  }
  search.addEventListener("input", render);
  document.getElementById("downloadCsv").addEventListener("click", downloadCsv);
  render();
</script>
</body>
</html>
"""
    return (
        template.replace("__TITLE__", title)
        .replace(
            "__PAYLOAD__",
            json.dumps(payload, separators=(",", ":"), sort_keys=True),
        )
    )


def write_fixture() -> list[Path]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    summary = RESULTS_DIR / "summary.json"
    parameters = RESULTS_DIR / "analysis_parameters.json"
    methods = RESULTS_DIR / "methods.md"
    session = RESULTS_DIR / "sessionInfo.txt"
    qc = RESULTS_DIR / "QC_Summary.csv"
    metrics = RESULTS_DIR / "Minfi_Metrics.csv"
    plot = RESULTS_DIR / "Minfi_Volcano.html"
    table = RESULTS_DIR / "Minfi_Top_DMPs.html"
    consensus_table = RESULTS_DIR / "Intersection_Native_Consensus_DMPs.html"
    consensus_csv = RESULTS_DIR / "Intersection_Native_Consensus_DMPs.csv"
    crf_paths = []

    summary.write_text(json.dumps(SUMMARY, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    parameters.write_text(json.dumps(PARAMETERS, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    methods.write_text(
        "# Synthetic dashboard QA fixture\n\n"
        "These values exercise the IlluMeta interface and are not study results.\n",
        encoding="utf-8",
    )
    session.write_text("Synthetic fixture; no R analysis was run.\n", encoding="utf-8")
    with qc.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(QC_ROWS)
    with metrics.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(METRICS_ROWS)
    with consensus_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(DMP_COLUMNS)
        writer.writerows(DMP_ROWS)
    for filename, (columns, row) in sorted(CRF_FIXTURES.items()):
        path = RESULTS_DIR / filename
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(columns)
            writer.writerow(row)
        crf_paths.append(path)
    plot.write_text(interactive_plot_html(), encoding="utf-8")
    table.write_text(
        interactive_table_html("Synthetic Minfi top DMPs"),
        encoding="utf-8",
    )
    consensus_table.write_text(
        interactive_table_html("Synthetic native consensus DMPs"),
        encoding="utf-8",
    )
    return [
        summary,
        parameters,
        methods,
        session,
        qc,
        metrics,
        plot,
        table,
        consensus_table,
        consensus_csv,
        *crf_paths,
    ]


def normalize_html() -> None:
    text = HTML_PATH.read_text(encoding="utf-8")
    text, count = re.subn(
        r'const RESULTS_PATH = "[^"]*";',
        'const RESULTS_PATH = "Case_vs_Control_results";',
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError("dashboard RESULTS_PATH marker was not found exactly once")
    disclosure = (
        '<div class="callout warning" role="note"><strong>Synthetic interface fixture.</strong> '
        "The displayed values test dashboard behavior and are not study results.</div>"
    )
    marker = '\t<div class="container">\n\t    <div class="jump-bar">'
    if marker not in text:
        raise RuntimeError("dashboard main-container marker is missing")
    text = text.replace(
        marker,
        '\t<div class="container">\n\t    ' + disclosure + '\n\t    <div class="jump-bar">',
        1,
    )
    text = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
    if str(ROOT) in text or "/tmp/" in text:
        raise RuntimeError("dashboard contains a local absolute path")
    if re.search(r'https?://', text):
        raise RuntimeError("dashboard contains a remote dependency")
    HTML_PATH.write_text(text, encoding="utf-8")


def build() -> None:
    fixture_files = write_fixture()
    illumeta.generate_dashboard(str(RESULTS_DIR), "Case", "Control")
    normalize_html()

    timestamp = datetime.fromtimestamp(SOURCE_DATE_EPOCH, UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    payload = {
        "artifact": str(HTML_PATH.relative_to(ROOT)),
        "artifact_sha256": sha256(HTML_PATH),
        "fixture_only": True,
        "generated_at": timestamp,
        "generator": str(Path(__file__).resolve().relative_to(ROOT)),
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in sorted(fixture_files)
        },
        "study_result": False,
        "source_date_epoch": SOURCE_DATE_EPOCH,
    }
    MANIFEST_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for path in [*fixture_files, HTML_PATH, MANIFEST_PATH]:
        os.utime(path, (SOURCE_DATE_EPOCH, SOURCE_DATE_EPOCH))


if __name__ == "__main__":
    build()
