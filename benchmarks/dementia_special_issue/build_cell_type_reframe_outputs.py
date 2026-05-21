#!/usr/bin/env python3
"""Build Lane 3 cell-type-aware dementia manuscript/report/figure outputs.

This integration script intentionally does not edit analysis configs. It reads
completed IlluMeta source artifacts and writes derived manuscript-facing files
with guarded language for a Methods/YMETH-compatible neuroepigenomics framing.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd


COHORT_ORDER = ["GSE208623", "GSE125895", "GSE134379", "GSE284764", "GSE105109", "GSE66351"]


def metric_csv(path: Path) -> dict[str, str]:
    df = pd.read_csv(path)
    if {"metric", "value"}.issubset(df.columns):
        return {str(r["metric"]): str(r["value"]) for _, r in df.iterrows()}
    return {}


def load_json(path: Path) -> dict:
    with path.open() as fh:
        return json.load(fh)


def source_paths(source_root: Path) -> dict[str, Path]:
    return {
        "cohort": source_root / "benchmarks/dementia_special_issue/cross_cohort/cohort_summary.tsv",
        "branch": source_root / "benchmarks/dementia_special_issue/cross_cohort/branch_metrics.tsv",
        "overlap": source_root / "benchmarks/dementia_special_issue/cross_cohort/consensus_pairwise_overlap.tsv",
        "draft": source_root / "benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_draft.md",
        "report": source_root / "benchmarks/dementia_special_issue/execution_report.md",
        "gse306226_summary": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/summary.json",
        "gse306226_groups": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Input_Group_Distribution.csv",
        "gse306226_strict_metrics": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Intersection_Comparison_Metrics.csv",
        "gse306226_native_metrics": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Intersection_Native_Comparison_Metrics.csv",
        "gse306226_qc": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/QC_Summary.csv",
        "gse306226_correction": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Correction_Adequacy_Summary.csv",
        "gse306226_minfi_metrics": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Minfi_Metrics.csv",
        "gse306226_sesame_metrics": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Sesame_Metrics.csv",
        "gse306226_native_pipeline_metrics": source_root / "projects/GSE306226/Neurons_vs_Microglia_results/Sesame_Native_Metrics.csv",
        "gse66351_summary": source_root / "projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare/summary.json",
    }


def require_inputs(paths: dict[str, Path]) -> None:
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required Lane 3 source artifacts:\n" + "\n".join(missing))


def build_summary(paths: dict[str, Path]) -> pd.DataFrame:
    cohort = pd.read_csv(paths["cohort"], sep="\t")
    branch = pd.read_csv(paths["branch"], sep="\t")
    overlap = pd.read_csv(paths["overlap"], sep="\t")
    g306 = load_json(paths["gse306226_summary"])
    g663 = load_json(paths["gse66351_summary"])
    groups = pd.read_csv(paths["gse306226_groups"])
    strict = metric_csv(paths["gse306226_strict_metrics"])
    native = metric_csv(paths["gse306226_native_metrics"])
    qc = metric_csv(paths["gse306226_qc"])
    correction = metric_csv(paths["gse306226_correction"])
    minfi_metrics = metric_csv(paths["gse306226_minfi_metrics"])
    sesame_metrics = metric_csv(paths["gse306226_sesame_metrics"])
    sesame_native_metrics = metric_csv(paths["gse306226_native_pipeline_metrics"])

    strict_totals = cohort["consensus_strict_up"] + cohort["consensus_strict_down"]
    native_totals = cohort["consensus_native_up"] + cohort["consensus_native_down"]
    zero_strict = int((strict_totals == 0).sum())
    zero_native = int((native_totals == 0).sum())
    lambda_triggered = int((cohort["lambda_guard_status"] == "triggered").sum())
    ad_total = int(cohort["n_AD"].sum())
    control_total = int(cohort["n_control"].sum())
    nonzero_overlap = overlap[overlap["overlap_n"] > 0].copy()

    rows = [
        {
            "axis": "AD public-IDAT dementia cohorts",
            "source": "benchmarks/dementia_special_issue/cross_cohort/cohort_summary.tsv",
            "n_control": control_total,
            "n_test": ad_total,
            "strict_consensus_dmps": int(strict_totals.sum()),
            "native_consensus_dmps": int(native_totals.sum()),
            "lambda_guard": f"{lambda_triggered}/{len(cohort)} cohorts triggered primary lambda guard",
            "branch_concordance": f"strict median logFC r={branch[branch.consensus_mode == 'strict'].logFC_correlation.median():.3f}; native median logFC r={branch[branch.consensus_mode == 'native'].logFC_correlation.median():.3f}",
            "interpretation": f"Six AD/dementia cohorts are workflow evidence with heterogeneity: {zero_strict} strict-zero and {zero_native} native-zero consensus cohorts, plus only {len(nonzero_overlap)} nonzero pairwise CpG-overlap rows.",
        },
        {
            "axis": "GSE66351 bulk AD sensitivity lane",
            "source": "projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare/summary.json",
            "n_control": int(g663["n_con"]),
            "n_test": int(g663["n_test"]),
            "strict_consensus_dmps": int(g663["intersect_up"] + g663["intersect_down"]),
            "native_consensus_dmps": int(g663["intersect_native_up"] + g663["intersect_native_down"]),
            "lambda_guard": f"primary={g663['primary_lambda_guard_status']}; native branch is treated as guarded sensitivity evidence",
            "branch_concordance": "strict consensus absent; native consensus exists but is not promoted to strict replication",
            "interpretation": "Bulk frontal/temporal cortex AD run supports heterogeneity/guardrail narrative rather than a universal AD CpG panel.",
        },
        {
            "axis": "GSE306226 Neurons_vs_Microglia reference axis",
            "source": "projects/GSE306226/Neurons_vs_Microglia_results",
            "n_control": int(g306["n_con"]),
            "n_test": int(g306["n_test"]),
            "strict_consensus_dmps": int(g306["intersect_up"] + g306["intersect_down"]),
            "native_consensus_dmps": int(g306["intersect_native_up"] + g306["intersect_native_down"]),
            "lambda_guard": f"primary={g306['primary_lambda_guard_status']}; Minfi={minfi_metrics.get('lambda_guard_status')}; SeSAMe={sesame_metrics.get('lambda_guard_status')}; native={sesame_native_metrics.get('lambda_guard_status')}",
            "branch_concordance": f"strict logFC r={float(strict['logFC_correlation']):.3f}, Jaccard={float(strict['jaccard_overlap']):.3f}; native logFC r={float(native['logFC_correlation']):.3f}, Jaccard={float(native['jaccard_overlap']):.3f}",
            "interpretation": f"Neural cell-type positive-control/reference axis ({'; '.join(groups['group'].astype(str)+' n='+groups['n'].astype(str))}); QC passed {qc.get('Samples_passed_QC')}/{qc.get('Total_samples_input')} samples, but lambda guard prevents unqualified biological-discovery claims. Correction preservation score={float(correction.get('preservation_score', 'nan')):.3f}.",
        },
    ]
    return pd.DataFrame(rows)


def guarded_markdown(summary: pd.DataFrame, generated_at: str) -> str:
    row_ad = summary.iloc[0]
    row_g663 = summary.iloc[1]
    row_g306 = summary.iloc[2]
    return textwrap.dedent(f"""
    # Cell-type-aware neuroepigenomics reframe for the dementia IlluMeta package

    Generated: {generated_at}

    ## Source-grounded conclusion

    The current package is strongest as a Methods/YMETH-compatible neuroepigenomics workflow paper, not as a universal Alzheimer disease biomarker paper. The completed AD/dementia public-IDAT analyses cover {row_ad.n_control} controls and {row_ad.n_test} AD/test samples across blood, bulk brain, EPIC brain, BS-only brain, and frontal/temporal cortex contexts. They produced {row_ad.strict_consensus_dmps} strict and {row_ad.native_consensus_dmps} native consensus DMPs in aggregate, but the evidence is heterogeneous: {row_ad.lambda_guard}; {row_ad.interpretation[0].lower() + row_ad.interpretation[1:]}

    ## Guarded AD evidence

    - `GSE125895` remains the cleanest brain AD anchor because it has nonzero strict/native consensus and primary lambda guard `ok`.
    - `GSE284764` provides recent EPIC prefrontal-cortex support, but its primary lambda guard triggered, so it remains guarded support.
    - `GSE134379` and `GSE105109` completed branch-level DMP/DMR outputs but produced no consensus DMPs.
    - `GSE66351` should be described as a bulk-cortex sensitivity lane: {row_g663.strict_consensus_dmps} strict and {row_g663.native_consensus_dmps} native consensus DMPs; {row_g663.interpretation}

    ## Neural cell-type reference axis

    `GSE306226/Neurons_vs_Microglia_results` provides the local positive-control/reference axis for cell-type-aware framing. It compares 20 neurons with 20 microglia, passed sample QC, and produced {row_g306.strict_consensus_dmps} strict plus {row_g306.native_consensus_dmps} native Minfi/SeSAMe consensus DMPs. Branch concordance was {row_g306.branch_concordance}. Because {row_g306.lambda_guard}, it should be used as a reference axis demonstrating that IlluMeta can recover a large neural cell-type contrast, not as an unqualified discovery claim.

    ## Manuscript claim guardrail

    Defensible central claim: IlluMeta exposes context dependence in public neurodegeneration methylation data by making raw-IDAT import, dual-pipeline agreement, consensus DMP/DMR outputs, and lambda/QC guards auditable. The manuscript should explicitly avoid phrases such as "universal AD biomarker" or "validated AD signature" unless tied to a specific cohort/context and guard status.
    """).strip() + "\n"


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    if start not in text or end not in text.split(start, 1)[1]:
        return text
    before, rest = text.split(start, 1)
    _old, after = rest.split(end, 1)
    return before + start + "\n\n" + replacement.strip() + "\n\n" + end + after


def revise_markdown(original: str, summary: pd.DataFrame, generated_at: str) -> str:
    title = "# IlluMeta exposes cohort and neural cell-type dependence in public dementia methylation IDAT reanalysis"
    original = original.replace(original.splitlines()[0], title, 1)
    old_abstract = original.split("## Abstract\n", 1)[1].split("\n\n## Keywords", 1)[0]
    row_g306 = summary.iloc[2]
    new_abstract = (
        "Public Alzheimer disease methylation IDAT datasets can benchmark reproducible EWAS workflows, "
        "but secondary reuse often fails at import, quality-control, covariate, cell-composition, and reporting boundaries. "
        "We reframed the IlluMeta dementia package as a Methods/YMETH-compatible, cell-type-aware neuroepigenomics workflow evaluation. "
        "Six completed AD/dementia public-IDAT cohorts showed strong context dependence: GSE208623, GSE125895, and GSE284764 produced nonzero strict consensus DMPs, whereas GSE134379 and GSE105109 produced no consensus DMPs and GSE66351 produced native-only guarded sensitivity evidence. "
        f"As a neural cell-type reference axis, GSE306226 Neurons_vs_Microglia produced {row_g306.strict_consensus_dmps} strict and {row_g306.native_consensus_dmps} native Minfi/SeSAMe consensus DMPs across 20 neurons and 20 microglia, but its lambda guard triggered. "
        "IlluMeta therefore supports auditable public-IDAT reuse by preserving successful signal recovery and negative/guarded evidence, while preventing overclaiming of universal Alzheimer disease methylation biomarkers."
    )
    original = original.replace(old_abstract, new_abstract)
    insert = "\n\n" + guarded_markdown(summary, generated_at).replace("# Cell-type-aware neuroepigenomics reframe for the dementia IlluMeta package", "## Cell-type-aware neuroepigenomics reframe")
    if "## Cell-type-aware neuroepigenomics reframe" not in original:
        original = original.replace("\n## Ethics Statement", insert + "\n\n## Ethics Statement")
    data_availability = """The raw public datasets are available from GEO under GSE208623 [3], GSE125895 [6], GSE134379 [7], GSE284764 [9], GSE105109 [15], GSE66351 [17], and the metadata-blocked GSE197305 [13]. All derived analysis artifacts are local to the IlluMeta project workspace. The dataset triage table is `benchmarks/dementia_special_issue/dataset_triage.tsv`. Full result directories are `projects/GSE208623/AD_vs_Control_full_long_results`, `projects/GSE125895/AD_vs_Control_all_regions_results_lcC`, `projects/GSE134379/AD_vs_Control_all_regions_results_lcC`, `projects/GSE284764/AD_vs_Ctl_PFC_results_lcC_run2`, `projects/GSE105109/AD_vs_Control_BS_only_results_lcC`, and `projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare`; the neural cell-type reference axis is `projects/GSE306226/Neurons_vs_Microglia_results`. Cross-cohort source tables are stored under `benchmarks/dementia_special_issue/cross_cohort`, and the cell-type reframe source table is `benchmarks/dementia_special_issue/cell_type_reframe_summary.tsv`. Figure 2, source data, legend, and manifest are stored under `benchmarks/dementia_special_issue/figures/FIG2`; the cell-type-aware Figure 3 package is stored under `benchmarks/dementia_special_issue/figures/FIG3`."""
    figure_legends = """**Figure 2. Cross-cohort IlluMeta evidence from completed public dementia methylation IDAT analyses.** A, sample composition for the six completed cohorts. B, strict and native Minfi/SeSAMe consensus DMP counts, split by direction. C, branch concordance measured by logFC correlation and significant-set Jaccard overlap; zero-consensus cohorts are annotated as no DMP. D, DMR row counts from Minfi, strict SeSAMe, and native SeSAMe branches. E, nonzero pairwise consensus-CpG overlaps across cohorts; zero-overlap pairs remain in the source workbook. F, strict-consensus CpGs shared between the 450K brain validation cohort and the EPIC prefrontal-cortex extension, with mean logFC values from each cohort.

**Figure 3. Cell-type-aware interpretation axis for the IlluMeta dementia reframe.** Consensus DMP counts are shown for the aggregate six-cohort AD/dementia public-IDAT package, the guarded GSE66351 bulk-cortex sensitivity lane, and the GSE306226 neuron-versus-microglia neural cell-type reference axis. The panel shows scale and context dependence and should not be interpreted as a universal Alzheimer disease biomarker figure."""
    original = replace_between(original, "## Data and code availability", "## Funding", data_availability)
    original = replace_between(original, "## Figure legends", "## References", figure_legends)
    return original


def revise_report(original: str, summary: pd.DataFrame, generated_at: str) -> str:
    block = "\n\n" + guarded_markdown(summary, generated_at).replace("# Cell-type-aware neuroepigenomics reframe for the dementia IlluMeta package", "## Cell-type-aware neuroepigenomics reframe update")
    if "## Cell-type-aware neuroepigenomics reframe update" not in original:
        original = original.replace("\n## Figure package", block + "\n\n## Figure package")
    return original


def make_figure(summary: pd.DataFrame, fig_dir: Path, generated_at: str) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    labels = ["AD cohorts\naggregate", "GSE66351\nbulk sensitivity", "GSE306226\nneuron/microglia"]
    strict = summary["strict_consensus_dmps"].astype(int).to_numpy()
    native = summary["native_consensus_dmps"].astype(int).to_numpy()
    x = range(len(labels))
    plt.rcParams.update({"font.size": 9})
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar([i - 0.18 for i in x], strict, width=0.36, label="strict consensus", color="#5B7F95")
    ax.bar([i + 0.18 for i in x], native, width=0.36, label="native consensus", color="#D08B35")
    for i, val in enumerate(strict):
        ax.text(i - 0.18, val + max(native) * 0.02, f"{val:,}", ha="center", va="bottom", fontsize=8)
    for i, val in enumerate(native):
        ax.text(i + 0.18, val + max(native) * 0.02, f"{val:,}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Consensus DMPs")
    ax.set_title("Cell-type-aware interpretation axis for IlluMeta dementia reframe", loc="left", fontweight="bold")
    ax.text(0.01, 0.98, "GSE306226 is a neural cell-type reference axis; AD cohorts remain heterogeneous/guarded.", transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#47515E")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D4D7DB", linewidth=0.6, alpha=0.7)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.88))
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(fig_dir / f"figure3_cell_type_reframe.{ext}", dpi=300 if ext == "png" else None)
    plt.close(fig)
    manifest = {
        "figure_id": "figure3_cell_type_reframe",
        "generated_at": generated_at,
        "source_table": "benchmarks/dementia_special_issue/cell_type_reframe_summary.tsv",
        "outputs": [f"figure3_cell_type_reframe.{ext}" for ext in ["png", "pdf", "svg"]],
        "claim_guardrail": "GSE306226 is a neural cell-type reference axis with lambda guard triggered; AD cohorts are heterogeneous and not universal biomarker evidence.",
    }
    (fig_dir / "figure3_cell_type_reframe_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (fig_dir / "figure3_cell_type_reframe_legend.md").write_text(
        "**Figure 3. Cell-type-aware interpretation axis for the IlluMeta dementia reframe.** "
        "Consensus DMP counts are shown for the aggregate six-cohort AD/dementia public-IDAT package, "
        "the guarded GSE66351 bulk-cortex sensitivity lane, and the GSE306226 neuron-versus-microglia neural cell-type reference axis. "
        "The panel is intended to show scale and context dependence, not to promote a universal AD biomarker set.\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=Path.cwd(), help="checkout containing source project artifacts")
    parser.add_argument("--output-root", type=Path, default=Path.cwd(), help="checkout/package root for derived outputs")
    args = parser.parse_args()
    paths = source_paths(args.source_root.resolve())
    require_inputs(paths)
    generated_at = datetime.now(timezone.utc).isoformat()
    out_base = args.output_root.resolve() / "benchmarks/dementia_special_issue"
    (out_base / "manuscript").mkdir(parents=True, exist_ok=True)
    (out_base / "figures/FIG3").mkdir(parents=True, exist_ok=True)
    summary = build_summary(paths)
    summary_path = out_base / "cell_type_reframe_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    (out_base / "cell_type_reframe_report.md").write_text(guarded_markdown(summary, generated_at))
    original_draft = paths["draft"].read_text()
    (out_base / "manuscript/dementia_special_issue_full_draft.md").write_text(revise_markdown(original_draft, summary, generated_at))
    original_report = paths["report"].read_text()
    (out_base / "execution_report.md").write_text(revise_report(original_report, summary, generated_at))
    make_figure(summary, out_base / "figures/FIG3", generated_at)
    print(f"Wrote {summary_path}")
    print(f"Wrote {out_base / 'cell_type_reframe_report.md'}")
    print(f"Wrote {out_base / 'manuscript/dementia_special_issue_full_draft.md'}")
    print(f"Wrote {out_base / 'execution_report.md'}")
    print(f"Wrote {out_base / 'figures/FIG3'}")


if __name__ == "__main__":
    main()
