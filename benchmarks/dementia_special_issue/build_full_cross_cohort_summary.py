#!/usr/bin/env python3
"""Summarize completed full IlluMeta dementia analyses across cohorts."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "benchmarks" / "dementia_special_issue" / "cross_cohort"

COHORTS = {
    "GSE208623": {
        "label": "Peripheral blood proof-of-work",
        "contrast": "AD_vs_Control",
        "platform": "EPIC",
        "tissue": "peripheral blood",
        "result_dir": ROOT / "projects" / "GSE208623" / "AD_vs_Control_full_long_results",
        "dashboard": ROOT / "projects" / "GSE208623" / "AD_vs_Control_full_long_results_index.html",
    },
    "GSE125895": {
        "label": "Brain validation",
        "contrast": "AD_vs_Control_all_regions",
        "platform": "450K",
        "tissue": "brain; ERC/DLPFC/CRB/HIPPO",
        "result_dir": ROOT / "projects" / "GSE125895" / "AD_vs_Control_all_regions_results_lcC",
        "dashboard": ROOT / "projects" / "GSE125895" / "AD_vs_Control_all_regions_results_lcC_index.html",
    },
    "GSE134379": {
        "label": "Large brain extension",
        "contrast": "AD_vs_Control_all_regions",
        "platform": "450K",
        "tissue": "brain; CBL/MTG",
        "result_dir": ROOT / "projects" / "GSE134379" / "AD_vs_Control_all_regions_results_lcC",
        "dashboard": ROOT / "projects" / "GSE134379" / "AD_vs_Control_all_regions_results_lcC_index.html",
    },
    "GSE284764": {
        "label": "EPIC PFC brain extension",
        "contrast": "AD_vs_Ctl_PFC",
        "platform": "EPIC",
        "tissue": "brain; PFC",
        "result_dir": ROOT / "projects" / "GSE284764" / "AD_vs_Ctl_PFC_results_lcC_run2",
        "dashboard": ROOT / "projects" / "GSE284764" / "AD_vs_Ctl_PFC_results_lcC_run2_index.html",
    },
    "GSE105109": {
        "label": "BS-only brain 450K extension",
        "contrast": "AD_vs_Control_BS_only",
        "platform": "450K",
        "tissue": "brain; entorhinal cortex/cerebellum; BS-only",
        "result_dir": ROOT / "projects" / "GSE105109" / "AD_vs_Control_BS_only_results_lcC",
        "dashboard": ROOT / "projects" / "GSE105109" / "AD_vs_Control_BS_only_results_lcC_index.html",
    },
    "GSE66351": {
        "label": "Brain bulk 450K extension",
        "contrast": "AD_vs_CTRL_bulk",
        "platform": "450K",
        "tissue": "brain; frontal/temporal cortex bulk",
        "result_dir": ROOT / "projects" / "GSE66351" / "AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare",
        "dashboard": ROOT / "projects" / "GSE66351" / "AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare_index.html",
    },
}

DMP_TABLES = {
    "minfi": "Minfi_DMPs_full.csv",
    "sesame_strict": "Sesame_DMPs_full.csv",
    "sesame_native": "Sesame_Native_DMPs_full.csv",
    "consensus_strict": "Intersection_Consensus_DMPs.csv",
    "consensus_native": "Intersection_Native_Consensus_DMPs.csv",
}

DMR_TABLES = {
    "minfi": "Minfi_DMRs.csv",
    "sesame_strict": "Sesame_DMRs.csv",
    "sesame_native": "Sesame_Native_DMRs.csv",
}

METRIC_TABLES = {
    "strict": "Intersection_Comparison_Metrics.csv",
    "native": "Intersection_Native_Comparison_Metrics.csv",
}

PIPELINE_METRIC_TABLES = {
    "minfi": "Minfi_Metrics.csv",
    "sesame_strict": "Sesame_Metrics.csv",
    "sesame_native": "Sesame_Native_Metrics.csv",
}

PIPELINE_LABELS = {
    "minfi": "Minfi",
    "sesame_strict": "SeSAMe strict",
    "sesame_native": "SeSAMe native",
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    return max(0, int(sum(1 for _ in path.open("rb")) - 1))


def read_metric_table(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    if {"metric", "value"} - set(df.columns):
        return {}
    return {str(row.metric): row.value for row in df.itertuples(index=False)}


def clean_value(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none", "null"}:
        return ""
    return text


def split_terms(value: object) -> list[str]:
    text = clean_value(value)
    if not text:
        return []
    return [term for term in (part.strip() for part in text.split(";")) if term]


def join_terms(terms: list[str]) -> str:
    return ";".join(terms)


def cell_terms(terms: list[str]) -> list[str]:
    return [term for term in terms if term.startswith("Cell_")]


def sv_terms(terms: list[str]) -> list[str]:
    return [term for term in terms if term.startswith("SV")]


def read_cell_deconvolution(result_dir: Path) -> dict[str, dict[str, str]]:
    path = result_dir / "Cell_Deconvolution_Summary.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    out: dict[str, dict[str, str]] = {}
    for row in df.to_dict(orient="records"):
        pipeline = clean_value(row.get("Pipeline"))
        if pipeline == "Sesame":
            key = "sesame_strict"
        elif pipeline == "Minfi":
            key = "minfi"
        else:
            key = pipeline.lower().replace(" ", "_")
        out[key] = {
            "method": clean_value(row.get("Method")),
            "tissue": clean_value(row.get("Tissue")),
            "k": clean_value(row.get("K")),
            "sample_count": clean_value(row.get("Sample_Count")),
            "cell_types": clean_value(row.get("Cell_Types")),
        }
    if "sesame_strict" in out and "sesame_native" not in out:
        out["sesame_native"] = dict(out["sesame_strict"])
    return out


def read_cell_group_association(result_dir: Path) -> dict[str, str]:
    path = result_dir / "Cell_Group_Association.csv"
    if not path.exists():
        return {"max_cell_eta2": "", "high_eta_cell_terms": "", "top_cell_group_association": ""}
    df = pd.read_csv(path)
    if df.empty or "Eta_Squared" not in df.columns:
        return {"max_cell_eta2": "", "high_eta_cell_terms": "", "top_cell_group_association": ""}
    df["Eta_Squared"] = pd.to_numeric(df["Eta_Squared"], errors="coerce")
    df = df.dropna(subset=["Eta_Squared"]).sort_values("Eta_Squared", ascending=False)
    if df.empty:
        return {"max_cell_eta2": "", "high_eta_cell_terms": "", "top_cell_group_association": ""}
    high = df.loc[df["Eta_Squared"] >= 0.5]
    top = df.head(3)
    return {
        "max_cell_eta2": f"{float(df['Eta_Squared'].iloc[0]):.6g}",
        "high_eta_cell_terms": ";".join(high["CellType"].astype(str)) if not high.empty else "",
        "top_cell_group_association": ";".join(
            f"{row.CellType}:eta2={float(row.Eta_Squared):.3g}" for row in top.itertuples(index=False)
        ),
    }


def summarize_branch_adjustment(
    gse: str,
    label: str,
    result_dir: Path,
    branch: str,
    filename: str,
    deconv: dict[str, dict[str, str]],
    cell_assoc: dict[str, str],
) -> dict[str, object]:
    metrics = read_metric_table(result_dir / filename)
    used = split_terms(metrics.get("covariates_used"))
    dropped = split_terms(metrics.get("dropped_covariates"))
    sv_used = split_terms(metrics.get("sv_used"))
    cells_used = cell_terms(used)
    dropped_cells = cell_terms(dropped)
    metadata_used = [term for term in used if not term.startswith("Cell_")]
    dropped_sv = sv_terms(dropped)
    deconv_row = deconv.get(branch, {})
    method = deconv_row.get("method", "")
    if method:
        cell_mode = "reference_free_latent" if "RefFree" in method else "reference_based_or_custom"
    elif cells_used:
        cell_mode = "cell_terms_from_metrics"
    else:
        cell_mode = "not_retained_or_unavailable"
    return {
        "GSE": gse,
        "label": label,
        "branch": PIPELINE_LABELS.get(branch, branch),
        "cell_adjustment_mode": cell_mode,
        "cell_deconvolution_method": method,
        "cell_deconvolution_tissue": deconv_row.get("tissue", ""),
        "cell_deconvolution_k": deconv_row.get("k", ""),
        "cell_deconvolution_sample_count": deconv_row.get("sample_count", ""),
        "cell_terms_used": join_terms(cells_used),
        "n_cell_terms_used": len(cells_used),
        "metadata_covariates_used": join_terms(metadata_used),
        "n_metadata_covariates_used": len(metadata_used),
        "sv_used": join_terms(sv_used),
        "n_sv_used": len(sv_used),
        "batch_method_applied": clean_value(metrics.get("batch_method_applied")),
        "batch_candidate": clean_value(metrics.get("batch_candidate")),
        "tier3_batch": clean_value(metrics.get("tier3_batch")),
        "primary_result_mode": clean_value(metrics.get("primary_result_mode")),
        "lambda_guard_status": clean_value(metrics.get("lambda_guard_status")),
        "dropped_cell_terms": join_terms(dropped_cells),
        "dropped_sv_terms": join_terms(dropped_sv),
        "dropped_covariates": join_terms(dropped),
        "max_cell_eta2": cell_assoc["max_cell_eta2"],
        "high_eta_cell_terms": cell_assoc["high_eta_cell_terms"],
        "top_cell_group_association": cell_assoc["top_cell_group_association"],
    }


def compact_by_branch(rows: list[dict[str, object]], key: str) -> str:
    parts = []
    for row in rows:
        value = clean_value(row.get(key))
        if value:
            parts.append(f"{row['branch']}={value}")
    return " | ".join(parts)


def build_adjustment_summary() -> pd.DataFrame:
    rows = []
    for gse, meta in COHORTS.items():
        result_dir = meta["result_dir"]
        deconv = read_cell_deconvolution(result_dir)
        cell_assoc = read_cell_group_association(result_dir)
        for branch, filename in PIPELINE_METRIC_TABLES.items():
            rows.append(
                summarize_branch_adjustment(
                    gse,
                    meta["label"],
                    result_dir,
                    branch,
                    filename,
                    deconv,
                    cell_assoc,
                )
            )
    return pd.DataFrame(rows)


def consensus_set(result_dir: Path, native: bool) -> set[str]:
    name = "Intersection_Native_Consensus_DMPs.csv" if native else "Intersection_Consensus_DMPs.csv"
    path = result_dir / name
    if not path.exists() or count_rows(path) == 0:
        return set()
    df = pd.read_csv(path, usecols=["CpG"])
    return set(df["CpG"].dropna().astype(str))


def build_cohort_summary(adjustment_summary: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    for gse, meta in COHORTS.items():
        out = meta["result_dir"]
        summary = read_json(out / "summary.json")
        adjustment_rows = []
        if adjustment_summary is not None and not adjustment_summary.empty:
            adjustment_rows = adjustment_summary.loc[adjustment_summary["GSE"] == gse].to_dict(orient="records")
        row = {
            "GSE": gse,
            "label": meta["label"],
            "contrast": meta["contrast"],
            "platform": meta["platform"],
            "tissue": meta["tissue"],
            "n_control": summary.get("n_con"),
            "n_AD": summary.get("n_test"),
            "minfi_sig_up": summary.get("minfi_up"),
            "minfi_sig_down": summary.get("minfi_down"),
            "sesame_strict_sig_up": summary.get("sesame_up"),
            "sesame_strict_sig_down": summary.get("sesame_down"),
            "sesame_native_sig_up": summary.get("sesame_native_up"),
            "sesame_native_sig_down": summary.get("sesame_native_down"),
            "consensus_strict_up": summary.get("intersect_up"),
            "consensus_strict_down": summary.get("intersect_down"),
            "consensus_native_up": summary.get("intersect_native_up"),
            "consensus_native_down": summary.get("intersect_native_down"),
            "primary_branch": summary.get("primary_branch"),
            "primary_result_mode": summary.get("primary_result_mode"),
            "lambda_guard_status": summary.get("primary_lambda_guard_status"),
            "crf_sample_tier": summary.get("crf_sample_tier"),
            "cell_adjustment_modes": compact_by_branch(adjustment_rows, "cell_adjustment_mode"),
            "cell_deconvolution_methods": compact_by_branch(adjustment_rows, "cell_deconvolution_method"),
            "cell_terms_used_by_branch": compact_by_branch(adjustment_rows, "cell_terms_used"),
            "dropped_cell_terms_by_branch": compact_by_branch(adjustment_rows, "dropped_cell_terms"),
            "sv_used_by_branch": compact_by_branch(adjustment_rows, "sv_used"),
            "batch_methods_by_branch": compact_by_branch(adjustment_rows, "batch_method_applied"),
            "high_eta_cell_terms": clean_value(adjustment_rows[0].get("high_eta_cell_terms")) if adjustment_rows else "",
            "max_cell_eta2": clean_value(adjustment_rows[0].get("max_cell_eta2")) if adjustment_rows else "",
            "dashboard": str(meta["dashboard"].relative_to(ROOT)),
        }
        for branch, filename in DMP_TABLES.items():
            row[f"{branch}_dmp_rows"] = count_rows(out / filename)
        for branch, filename in DMR_TABLES.items():
            row[f"{branch}_dmr_rows"] = count_rows(out / filename)
        for branch, filename in PIPELINE_METRIC_TABLES.items():
            metrics = read_metric_table(out / filename)
            row[f"{branch}_lambda"] = metrics.get("lambda")
            row[f"{branch}_lambda_guard_status"] = metrics.get("lambda_guard_status")
            row[f"{branch}_lambda_guard_action"] = metrics.get("lambda_guard_action")
        rows.append(row)
    return pd.DataFrame(rows)


def build_branch_metrics() -> pd.DataFrame:
    rows = []
    for gse, meta in COHORTS.items():
        out = meta["result_dir"]
        for mode, filename in METRIC_TABLES.items():
            metrics = read_metric_table(out / filename)
            rows.append(
                {
                    "GSE": gse,
                    "consensus_mode": mode,
                    "logFC_correlation": metrics.get("logFC_correlation"),
                    "jaccard_overlap": metrics.get("jaccard_overlap"),
                    "n_both": metrics.get("n_both"),
                }
            )
    return pd.DataFrame(rows)


def build_pairwise_overlap(native: bool) -> pd.DataFrame:
    mode = "native" if native else "strict"
    sets = {gse: consensus_set(meta["result_dir"], native=native) for gse, meta in COHORTS.items()}
    rows = []
    for left, right in combinations(COHORTS, 2):
        a = sets[left]
        b = sets[right]
        union = a | b
        inter = a & b
        rows.append(
            {
                "consensus_mode": mode,
                "left_GSE": left,
                "right_GSE": right,
                "left_n": len(a),
                "right_n": len(b),
                "overlap_n": len(inter),
                "jaccard": len(inter) / len(union) if union else 0.0,
                "overlap_cpgs": ";".join(sorted(inter)[:100]),
            }
        )
    return pd.DataFrame(rows)


def build_top_consensus() -> pd.DataFrame:
    rows = []
    for gse, meta in COHORTS.items():
        out = meta["result_dir"]
        for mode, filename in {
            "strict": "Intersection_Consensus_DMPs.csv",
            "native": "Intersection_Native_Consensus_DMPs.csv",
        }.items():
            path = out / filename
            if not path.exists() or count_rows(path) == 0:
                continue
            df = pd.read_csv(path)
            sort_cols = [c for c in ["adj.P.Val.selection", "adj.P.Val", "P.Value.selection", "P.Value"] if c in df.columns]
            if sort_cols:
                df = df.sort_values(sort_cols, na_position="last")
            for rank, (_, item) in enumerate(df.head(50).iterrows(), start=1):
                rows.append(
                    {
                        "GSE": gse,
                        "consensus_mode": mode,
                        "rank": rank,
                        "CpG": item.get("CpG"),
                        "Gene": item.get("Gene"),
                        "chr": item.get("chr"),
                        "pos": item.get("pos"),
                        "logFC_mean": item.get("logFC_mean"),
                        "Delta_Beta_Minfi": item.get("Delta_Beta.Minfi"),
                        "Delta_Beta_Sesame": item.get("Delta_Beta.Sesame"),
                        "adj_P_selection": item.get("adj.P.Val.selection"),
                        "adj_P_fisher": item.get("adj.P.Fisher"),
                    }
                )
    return pd.DataFrame(rows)


def build_top_dmrs() -> pd.DataFrame:
    rows = []
    for gse, meta in COHORTS.items():
        out = meta["result_dir"]
        for branch, filename in DMR_TABLES.items():
            path = out / filename
            if not path.exists() or count_rows(path) == 0:
                continue
            df = pd.read_csv(path)
            sort_cols = [col for col in ["p.adjust", "p.value"] if col in df.columns]
            if sort_cols:
                df = df.sort_values(sort_cols, na_position="last")
            for rank, (_, item) in enumerate(df.head(25).iterrows(), start=1):
                rows.append(
                    {
                        "GSE": gse,
                        "branch": branch,
                        "rank": rank,
                        "chr": item.get("chr"),
                        "start": item.get("start"),
                        "end": item.get("end"),
                        "n_probes": item.get("n"),
                        "Genes": item.get("Genes"),
                        "Regions": item.get("Regions"),
                        "Delta_Beta": item.get("Delta_Beta"),
                        "p_value": item.get("p.value"),
                        "p_adjust": item.get("p.adjust"),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adjustment_summary = build_adjustment_summary()
    cohort_summary = build_cohort_summary(adjustment_summary)
    branch_metrics = build_branch_metrics()
    overlap = pd.concat([build_pairwise_overlap(False), build_pairwise_overlap(True)], ignore_index=True)
    top_consensus = build_top_consensus()
    top_dmrs = build_top_dmrs()

    tables = {
        "cohort_summary": cohort_summary,
        "adjustment_summary": adjustment_summary,
        "branch_metrics": branch_metrics,
        "consensus_pairwise_overlap": overlap,
        "top_consensus_cpgs": top_consensus,
        "top_dmr_regions": top_dmrs,
    }
    for name, df in tables.items():
        df.to_csv(OUT_DIR / f"{name}.tsv", sep="\t", index=False)
        df.to_csv(OUT_DIR / f"{name}.csv", index=False)

    workbook = OUT_DIR / "dementia_full_cross_cohort_source_data.xlsx"
    with pd.ExcelWriter(workbook) as writer:
        for name, df in tables.items():
            df.to_excel(writer, sheet_name=name[:31], index=False)

    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "cohorts": {gse: str(meta["result_dir"].relative_to(ROOT)) for gse, meta in COHORTS.items()},
        "outputs": [str((OUT_DIR / f"{name}.tsv").relative_to(ROOT)) for name in tables]
        + [str(workbook.relative_to(ROOT))],
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
