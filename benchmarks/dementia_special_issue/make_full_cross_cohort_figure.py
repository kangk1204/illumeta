#!/usr/bin/env python3
"""Render a manuscript-grade cross-cohort IlluMeta dementia evidence figure."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmarks" / "dementia_special_issue"
CROSS = BASE / "cross_cohort"
FIG_DIR = BASE / "figures" / "FIG2"
FIG_ID = "figure2_dementia_full_cross_cohort_evidence"

COHORT_ORDER = ["GSE208623", "GSE125895", "GSE134379", "GSE284764", "GSE105109", "GSE66351"]
COHORT_LABELS = {
    "GSE208623": "GSE208623\nblood EPIC",
    "GSE125895": "GSE125895\nbrain 450K",
    "GSE134379": "GSE134379\nbrain 450K",
    "GSE284764": "GSE284764\nPFC EPIC",
    "GSE105109": "GSE105109\nBS brain 450K",
    "GSE66351": "GSE66351\nbulk brain 450K",
}

COLORS = {
    "AD": "#B4433F",
    "Control": "#2F7F7B",
    "up": "#B4433F",
    "down": "#2E5EAA",
    "strict": "#5B7F95",
    "native": "#D08B35",
    "minfi": "#406A9A",
    "sesame_strict": "#4E8B5B",
    "sesame_native": "#9A6A3A",
    "grid": "#D4D7DB",
    "text": "#1F2933",
}


def strip_trailing_whitespace(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def load_inputs() -> dict[str, pd.DataFrame]:
    tables = {
        "cohort": pd.read_csv(CROSS / "cohort_summary.tsv", sep="\t"),
        "adjustment": pd.read_csv(CROSS / "adjustment_summary.tsv", sep="\t"),
        "branch": pd.read_csv(CROSS / "branch_metrics.tsv", sep="\t"),
        "overlap": pd.read_csv(CROSS / "consensus_pairwise_overlap.tsv", sep="\t"),
        "top_consensus": pd.read_csv(CROSS / "top_consensus_cpgs.tsv", sep="\t"),
        "top_dmrs": pd.read_csv(CROSS / "top_dmr_regions.tsv", sep="\t"),
    }
    return tables


def shared_cpg_table() -> pd.DataFrame:
    columns = [
        "CpG",
        "Gene",
        "brain_450k_logFC_mean",
        "pfc_epic_logFC_mean",
        "brain_450k_delta_beta_minfi",
        "pfc_epic_delta_beta_minfi",
        "brain_450k_adj_P_selection",
        "pfc_epic_adj_P_selection",
    ]
    brain_450k = pd.read_csv(
        ROOT
        / "projects"
        / "GSE125895"
        / "AD_vs_Control_all_regions_results_lcC"
        / "Intersection_Consensus_DMPs.csv"
    )
    pfc_epic = pd.read_csv(
        ROOT
        / "projects"
        / "GSE284764"
        / "AD_vs_Ctl_PFC_results_lcC_run2"
        / "Intersection_Consensus_DMPs.csv"
    )
    shared = sorted(set(brain_450k["CpG"].astype(str)) & set(pfc_epic["CpG"].astype(str)))
    rows = []
    for cpg in shared:
        brn = brain_450k.loc[brain_450k["CpG"].astype(str) == cpg].iloc[0]
        pfc = pfc_epic.loc[pfc_epic["CpG"].astype(str) == cpg].iloc[0]
        gene = brn.get("Gene")
        if pd.isna(gene) or gene == "":
            gene = pfc.get("Gene")
        rows.append(
            {
                "CpG": cpg,
                "Gene": "" if pd.isna(gene) else gene,
                "brain_450k_logFC_mean": brn["logFC_mean"],
                "pfc_epic_logFC_mean": pfc["logFC_mean"],
                "brain_450k_delta_beta_minfi": brn["Delta_Beta.Minfi"],
                "pfc_epic_delta_beta_minfi": pfc["Delta_Beta.Minfi"],
                "brain_450k_adj_P_selection": brn["adj.P.Val.selection"],
                "pfc_epic_adj_P_selection": pfc["adj.P.Val.selection"],
            }
        )
    return pd.DataFrame(rows, columns=columns)


def set_panel_title(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", color=COLORS["text"], pad=8)


def style_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A8ADB4")
    ax.spines["bottom"].set_color("#A8ADB4")
    ax.tick_params(colors=COLORS["text"], labelsize=8.5)
    ax.grid(True, axis=grid_axis, color=COLORS["grid"], linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def panel_a(ax: plt.Axes, cohort: pd.DataFrame) -> pd.DataFrame:
    data = cohort.set_index("GSE").loc[COHORT_ORDER].reset_index()
    x = np.arange(len(data))
    ax.bar(x, data["n_control"], color=COLORS["Control"], width=0.62, label="Control")
    ax.bar(x, data["n_AD"], bottom=data["n_control"], color=COLORS["AD"], width=0.62, label="AD")
    totals = data["n_control"] + data["n_AD"]
    for i, total in enumerate(totals):
        ax.text(i, total + max(totals) * 0.03, f"n={int(total)}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([COHORT_LABELS[g] for g in data["GSE"]])
    ax.set_ylabel("Samples")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    set_panel_title(ax, "A. Completed public IDAT cohorts")
    style_axis(ax)
    return data[["GSE", "label", "platform", "tissue", "n_control", "n_AD"]]


def panel_b(ax: plt.Axes, cohort: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in cohort.iterrows():
        rows.extend(
            [
                {
                    "GSE": row["GSE"],
                    "mode": "strict",
                    "up": row["consensus_strict_up"],
                    "down": row["consensus_strict_down"],
                },
                {
                    "GSE": row["GSE"],
                    "mode": "native",
                    "up": row["consensus_native_up"],
                    "down": row["consensus_native_down"],
                },
            ]
        )
    data = pd.DataFrame(rows)
    x = np.arange(len(COHORT_ORDER))
    width = 0.32
    for offset, mode in [(-width / 2, "strict"), (width / 2, "native")]:
        sub = data[data["mode"] == mode].set_index("GSE").loc[COHORT_ORDER]
        xpos = x + offset
        ax.bar(xpos, sub["down"], width=width, color=COLORS["down"], alpha=0.88 if mode == "strict" else 0.55)
        ax.bar(
            xpos,
            sub["up"],
            width=width,
            bottom=sub["down"],
            color=COLORS["up"],
            alpha=0.88 if mode == "strict" else 0.55,
        )
        for i, total in enumerate(sub["up"] + sub["down"]):
            ax.text(xpos[i], total + 75, f"{int(total)}", ha="center", va="bottom", fontsize=7.8)
    ax.set_xticks(x)
    ax.set_xticklabels([COHORT_LABELS[g] for g in COHORT_ORDER])
    ax.set_ylabel("Consensus DMPs")
    ax.text(0.02, 0.98, "dark = strict, light = native", transform=ax.transAxes, ha="left", va="top", fontsize=7.8)
    ax.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, color=COLORS["up"], label="Up"),
            plt.Rectangle((0, 0), 1, 1, color=COLORS["down"], label="Down"),
        ],
        frameon=False,
        fontsize=8,
        loc="upper right",
    )
    set_panel_title(ax, "B. Dual-pipeline consensus DMPs")
    style_axis(ax)
    return data


def panel_c(ax: plt.Axes, branch: pd.DataFrame) -> pd.DataFrame:
    data = branch.copy()
    data["label"] = data["GSE"].map(COHORT_LABELS).str.replace("\n", " ", regex=False) + " " + data["consensus_mode"]
    x = np.arange(len(data))
    ax.scatter(x, data["logFC_correlation"], s=52, color=COLORS["strict"], label="logFC correlation", zorder=3)
    ax.scatter(x, data["jaccard_overlap"].fillna(0), s=52, color=COLORS["native"], marker="s", label="Jaccard overlap", zorder=3)
    for i, row in data.iterrows():
        if pd.isna(row["jaccard_overlap"]):
            ax.text(i, 0.04, "no\nDMP", ha="center", va="bottom", fontsize=7, color="#5C626B")
    ax.set_ylim(-0.03, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(data["label"], rotation=35, ha="right")
    ax.set_ylabel("Metric value")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    set_panel_title(ax, "C. Minfi/SeSAMe branch concordance")
    style_axis(ax)
    return data[["GSE", "consensus_mode", "logFC_correlation", "jaccard_overlap", "n_both"]]


def panel_d(ax: plt.Axes, cohort: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in cohort.iterrows():
        for branch, label in [
            ("minfi", "Minfi"),
            ("sesame_strict", "SeSAMe strict"),
            ("sesame_native", "SeSAMe native"),
        ]:
            rows.append({"GSE": row["GSE"], "branch": label, "dmr_rows": row[f"{branch}_dmr_rows"]})
    data = pd.DataFrame(rows)
    x = np.arange(len(COHORT_ORDER))
    width = 0.22
    branches = ["Minfi", "SeSAMe strict", "SeSAMe native"]
    colors = [COLORS["minfi"], COLORS["sesame_strict"], COLORS["sesame_native"]]
    for j, (branch, color) in enumerate(zip(branches, colors)):
        sub = data[data["branch"] == branch].set_index("GSE").loc[COHORT_ORDER]
        xpos = x + (j - 1) * width
        ax.bar(xpos, sub["dmr_rows"], width=width, color=color, label=branch)
    ax.set_xticks(x)
    ax.set_xticklabels([COHORT_LABELS[g] for g in COHORT_ORDER])
    ax.set_ylabel("DMR rows")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    set_panel_title(ax, "D. Regional methylation calls")
    style_axis(ax)
    return data


def panel_e(ax: plt.Axes, overlap: pd.DataFrame) -> pd.DataFrame:
    data = overlap.copy()
    data["plotted"] = data["overlap_n"] > 0
    plot_data = data.loc[data["plotted"]].copy()
    if plot_data.empty:
        plot_data = data.copy()
    data["pair"] = data["left_GSE"] + " vs " + data["right_GSE"]
    plot_data["pair"] = plot_data["left_GSE"] + " vs " + plot_data["right_GSE"]
    x = np.arange(len(plot_data))
    colors = [COLORS["strict"] if m == "strict" else COLORS["native"] for m in plot_data["consensus_mode"]]
    ax.bar(x, plot_data["overlap_n"], color=colors, width=0.66)
    for i, (_, row) in enumerate(plot_data.iterrows()):
        ax.text(i, row["overlap_n"] + 0.14, str(int(row["overlap_n"])), ha="center", va="bottom", fontsize=7.8)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_data["pair"] + "\n" + plot_data["consensus_mode"], rotation=30, ha="right")
    ax.set_ylabel("Shared consensus CpGs")
    set_panel_title(ax, "E. Cross-cohort CpG overlap")
    ax.text(
        0.02,
        0.95,
        "Only nonzero pairwise overlaps are plotted; zero-overlap pairs remain in source data.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.8,
        color="#47515E",
    )
    style_axis(ax)
    return data[
        [
            "consensus_mode",
            "left_GSE",
            "right_GSE",
            "left_n",
            "right_n",
            "overlap_n",
            "jaccard",
            "overlap_cpgs",
            "plotted",
        ]
    ]


def panel_f(ax: plt.Axes, shared: pd.DataFrame) -> pd.DataFrame:
    ax.axis("off")
    cols = ["CpG", "Gene", "brain_450k_logFC_mean", "pfc_epic_logFC_mean"]
    data = shared[cols].copy()
    data["brain_450k_logFC_mean"] = data["brain_450k_logFC_mean"].map(lambda x: f"{x:.2f}")
    data["pfc_epic_logFC_mean"] = data["pfc_epic_logFC_mean"].map(lambda x: f"{x:.2f}")
    if data.empty:
        data = pd.DataFrame([["No shared strict CpGs", "", "", ""]], columns=cols)
    table = ax.table(
        cellText=data.values,
        colLabels=["CpG", "Gene", "450K brain\nlogFC", "EPIC PFC\nlogFC"],
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=[0.28, 0.18, 0.22, 0.22],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.2)
    table.scale(1, 1.45)
    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#CBD1D8")
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#4B5563")
        else:
            cell.set_facecolor("#F8FAFC" if row % 2 else "#FFFFFF")
    set_panel_title(ax, "F. Brain cross-platform shared CpGs")
    shared_count = len(shared)
    if shared_count == 0:
        shared_note = "No strict-consensus CpGs are shared by GSE125895 and GSE284764."
    elif shared_count == 1:
        shared_note = "One strict-consensus CpG is shared by GSE125895 and GSE284764."
    else:
        shared_note = f"{shared_count} strict-consensus CpGs are shared by GSE125895 and GSE284764."
    ax.text(
        0.0,
        0.06,
        shared_note,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color="#47515E",
    )
    return shared


def save_figure(fig: plt.Figure, stem: Path) -> list[str]:
    outputs = []
    for ext in ("png", "pdf", "svg"):
        path = stem.with_suffix(f".{ext}")
        fig.savefig(path, dpi=320, bbox_inches="tight")
        if ext == "svg":
            strip_trailing_whitespace(path)
        outputs.append(str(path.relative_to(ROOT)))
    return outputs


def render_panel_only(name: str, draw_fn, *args) -> list[str]:
    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    draw_fn(ax, *args)
    fig.tight_layout()
    outputs = save_figure(fig, FIG_DIR / f"{FIG_ID}_{name}")
    plt.close(fig)
    return outputs


def write_legend(path: Path) -> None:
    doc = Document()
    doc.add_paragraph(
        "Figure 2. Cross-cohort IlluMeta evidence from completed public dementia methylation IDAT analyses. "
        "A, sample composition for the six completed cohorts. "
        "B, strict and native Minfi/SeSAMe consensus DMP counts, split by direction. "
        "C, branch concordance measured by logFC correlation and significant-set Jaccard overlap; zero-consensus cohorts are annotated as no DMP. "
        "D, DMR row counts from Minfi, strict SeSAMe, and native SeSAMe branches. "
        "E, nonzero pairwise consensus-CpG overlaps across cohorts; zero-overlap pairs are retained in the source workbook. "
        "F, strict-consensus CpGs shared between the 450K brain validation cohort and the EPIC PFC extension, with mean logFC values from each cohort."
    )
    doc.save(path)


def write_source_data(path: Path, sources: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(path) as writer:
        for name, df in sources.items():
            df.to_excel(writer, sheet_name=name[:31], index=False)
        pd.DataFrame(
            [
                {
                    "figure_id": FIG_ID,
                    "generated_at": datetime.now(timezone.utc).isoformat(),
                    "source_summary": "All plotted values are read from full IlluMeta result summaries and cross-cohort tables.",
                }
            ]
        ).to_excel(writer, sheet_name="metadata", index=False)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    tables = load_inputs()
    shared = shared_cpg_table()

    fig = plt.figure(figsize=(17.2, 11.6), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.95])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    ax_e = fig.add_subplot(grid[2, 0])
    ax_f = fig.add_subplot(grid[2, 1])

    sources = {
        "panel_A_samples": panel_a(ax_a, tables["cohort"]),
        "panel_B_consensus_counts": panel_b(ax_b, tables["cohort"]),
        "panel_C_branch_metrics": panel_c(ax_c, tables["branch"]),
        "panel_D_dmr_counts": panel_d(ax_d, tables["cohort"]),
        "panel_E_pairwise_overlap": panel_e(ax_e, tables["overlap"]),
        "panel_F_shared_cpgs": panel_f(ax_f, shared),
        "adjustment_summary": tables["adjustment"],
    }
    combined_outputs = save_figure(fig, FIG_DIR / FIG_ID)
    plt.close(fig)

    panel_outputs = []
    panel_outputs.extend(render_panel_only("panel_a", panel_a, tables["cohort"]))
    panel_outputs.extend(render_panel_only("panel_b", panel_b, tables["cohort"]))
    panel_outputs.extend(render_panel_only("panel_c", panel_c, tables["branch"]))
    panel_outputs.extend(render_panel_only("panel_d", panel_d, tables["cohort"]))
    panel_outputs.extend(render_panel_only("panel_e", panel_e, tables["overlap"]))
    panel_outputs.extend(render_panel_only("panel_f", panel_f, shared))

    source_path = FIG_DIR / f"{FIG_ID}_source_data.xlsx"
    legend_path = FIG_DIR / f"{FIG_ID}_legend.docx"
    manifest_path = FIG_DIR / f"{FIG_ID}_manifest.json"
    write_source_data(source_path, sources)
    write_legend(legend_path)

    manifest = {
        "figure_id": FIG_ID,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "script": str((BASE / "make_full_cross_cohort_figure.py").relative_to(ROOT)),
        "inputs": [
            str((CROSS / "cohort_summary.tsv").relative_to(ROOT)),
            str((CROSS / "adjustment_summary.tsv").relative_to(ROOT)),
            str((CROSS / "branch_metrics.tsv").relative_to(ROOT)),
            str((CROSS / "consensus_pairwise_overlap.tsv").relative_to(ROOT)),
            str((CROSS / "top_consensus_cpgs.tsv").relative_to(ROOT)),
            str((CROSS / "top_dmr_regions.tsv").relative_to(ROOT)),
        ],
        "outputs": combined_outputs + panel_outputs + [str(legend_path.relative_to(ROOT))],
        "source_data": [str(source_path.relative_to(ROOT))],
        "panel_claims": [
            "GSE208623, GSE125895, GSE134379, GSE284764, GSE105109, and GSE66351 completed full IlluMeta analyses from public IDAT files.",
            "GSE208623, GSE125895, GSE284764, and native GSE66351 contain nonzero dual-pipeline consensus DMP evidence; GSE134379 and GSE105109 have zero consensus DMPs.",
            "GSE125895 shows high Minfi/SeSAMe significant-set overlap; GSE134379 and GSE105109 have no consensus DMP set.",
            "All six completed cohorts contain DMR outputs in Minfi, strict SeSAMe, and native SeSAMe branches.",
            f"{len(shared)} strict consensus CpGs are shared between GSE125895 and GSE284764 in the current source tables.",
        ],
        "transformations": [
            "Consensus counts were read from IlluMeta summary.json files and cross-cohort summary tables.",
            "DMP overlap was computed on CpG identifiers from Intersection_Consensus_DMPs.csv and Intersection_Native_Consensus_DMPs.csv.",
            "DMR counts are row counts from branch-specific DMR CSV files.",
        ],
        "validation": {
            "script_regenerated_outputs": True,
            "formats": ["png", "pdf", "svg"],
            "source_workbook": str(source_path.relative_to(ROOT)),
        },
        "known_limitations": [
            "GSE134379 is a completed stress-test/extension cohort but is not a clean positive DMP replication result because consensus DMP count is zero and the lambda guard triggered.",
            "GSE284764 contributes a recent EPIC prefrontal-cortex brain extension with nonzero consensus DMPs but also lambda-guard triggering, so it should be presented as guarded supporting evidence.",
            "GSE105109 is a completed BS-only brain extension with DMR outputs but zero strict/native consensus DMPs.",
            "GSE66351 contributes native-consensus DMPs only; its native SeSAMe branch has lambda-guard triggering and should be interpreted conservatively.",
            "Cross-cohort CpG overlap is conservative because cohorts differ by tissue and array platform.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
