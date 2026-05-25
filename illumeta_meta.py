#!/usr/bin/env python3
"""Cross-cohort CpG meta-analysis for completed IlluMeta result folders."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_BRANCH_FILES = {
    "minfi": "Minfi_DMPs_full.csv",
    "sesame_strict": "Sesame_DMPs_full.csv",
    "sesame_native": "Sesame_Native_DMPs_full.csv",
}

BRANCH_ALIASES = {
    "sesame": "sesame_strict",
    "strict": "sesame_strict",
    "native": "sesame_native",
}

ANNOTATION_COLUMNS = ("Gene", "chr", "pos", "Region", "Island_Context")


@dataclass(frozen=True)
class MetaThresholds:
    min_cohorts: int = 3
    meta_fdr: float = 0.05
    min_direction_fraction: float = 0.70
    min_loo_direction_fraction: float = 0.80
    max_i2: float = 60.0
    min_abs_delta_beta: float = 0.01
    partial_conjunction_r: int = 3
    top_n: int = 1000


@dataclass(frozen=True)
class MetaCohort:
    cohort_id: str
    result_dir: Path
    label: str = ""
    platform: str = ""
    tissue: str = ""
    n_con: int = 0
    n_test: int = 0
    primary_result_mode: str = ""
    primary_lambda_guard_status: str = ""

    @property
    def sample_weight(self) -> float:
        n = int(self.n_con or 0) + int(self.n_test or 0)
        return float(n if n > 0 else 1)


def _log(message: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"na", "nan", "null", "none"}


def _first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        if not _is_missing(value):
            return str(value).strip()
    return ""


def _safe_float(value: object) -> float:
    if _is_missing(value):
        return math.nan
    try:
        out = float(str(value).strip())
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def _safe_int(value: object) -> int:
    val = _safe_float(value)
    return int(val) if math.isfinite(val) else 0


def _format_float(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    if isinstance(value, (float, int)):
        val = float(value)
        if not math.isfinite(val):
            return ""
        if abs(val) >= 1e4 or (0 < abs(val) < 1e-4):
            return f"{val:.6g}"
        return f"{val:.10g}"
    return str(value)


def _sniff_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        return "\t" if "\t" in sample else ","


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def _read_summary(result_dir: Path) -> dict[str, object]:
    summary_path = result_dir / "summary.json"
    if not summary_path.exists():
        return {}
    try:
        with summary_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _resolve_path(path_text: str, root: Path) -> Path:
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def _truthy_include(value: object) -> bool:
    if _is_missing(value):
        return True
    return str(value).strip().lower() not in {"0", "false", "no", "n", "exclude", "skip"}


def load_meta_manifest(manifest: Path, project_root: Path) -> list[MetaCohort]:
    delimiter = _sniff_delimiter(manifest)
    rows: list[dict[str, str]] = []
    with manifest.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Manifest has no header: {manifest}")
        for row in reader:
            if _truthy_include(row.get("include")):
                rows.append(row)
    cohorts: list[MetaCohort] = []
    for i, row in enumerate(rows, start=1):
        result_dir_text = _first_nonempty(
            row.get(col) for col in ("result_dir", "results_dir", "path", "dir")
        )
        if not result_dir_text:
            raise ValueError("Manifest must contain a result_dir/path column.")
        result_dir = _resolve_path(result_dir_text, project_root)
        summary = _read_summary(result_dir)
        cohort_id = _first_nonempty(
            row.get(col) for col in ("cohort", "cohort_id", "gse", "id", "name")
        ) or result_dir.parent.name or f"cohort_{i}"
        cohorts.append(
            MetaCohort(
                cohort_id=cohort_id,
                result_dir=result_dir,
                label=str(row.get("label") or ""),
                platform=str(row.get("platform") or ""),
                tissue=str(row.get("tissue") or ""),
                n_con=_safe_int(summary.get("n_con")),
                n_test=_safe_int(summary.get("n_test")),
                primary_result_mode=str(summary.get("primary_result_mode") or ""),
                primary_lambda_guard_status=str(summary.get("primary_lambda_guard_status") or ""),
            )
        )
    return cohorts


def load_positional_cohorts(result_dirs: list[str], project_root: Path) -> list[MetaCohort]:
    cohorts: list[MetaCohort] = []
    for i, result_dir_text in enumerate(result_dirs, start=1):
        result_dir = _resolve_path(result_dir_text, project_root)
        summary = _read_summary(result_dir)
        cohorts.append(
            MetaCohort(
                cohort_id=result_dir.parent.name or f"cohort_{i}",
                result_dir=result_dir,
                n_con=_safe_int(summary.get("n_con")),
                n_test=_safe_int(summary.get("n_test")),
                primary_result_mode=str(summary.get("primary_result_mode") or ""),
                primary_lambda_guard_status=str(summary.get("primary_lambda_guard_status") or ""),
            )
        )
    return cohorts


def parse_branch_selection(branches: str, branch_file_overrides: list[str] | None = None) -> dict[str, str]:
    selected: dict[str, str] = {}
    for raw in (branches or "").split(","):
        name = raw.strip().lower()
        if not name:
            continue
        name = BRANCH_ALIASES.get(name, name)
        if name not in DEFAULT_BRANCH_FILES:
            raise ValueError(f"Unknown branch '{raw}'. Valid branches: {', '.join(DEFAULT_BRANCH_FILES)}")
        selected[name] = DEFAULT_BRANCH_FILES[name]
    if not selected:
        selected = dict(DEFAULT_BRANCH_FILES)
    for override in branch_file_overrides or []:
        if "=" not in override:
            raise ValueError("--branch-file values must be NAME=FILENAME")
        name, filename = override.split("=", 1)
        name = BRANCH_ALIASES.get(name.strip().lower(), name.strip().lower())
        if name not in selected:
            raise ValueError(f"--branch-file override names an unselected branch: {name}")
        if not filename.strip():
            raise ValueError("--branch-file filename cannot be empty")
        selected[name] = filename.strip()
    return selected


def _normal_two_sided_p(z: float) -> float:
    if not math.isfinite(z):
        return math.nan
    return math.erfc(abs(z) / math.sqrt(2.0))


def _chi2_sf_even(stat: float, df: int) -> float:
    if not math.isfinite(stat) or stat < 0 or df <= 0 or df % 2 != 0:
        return math.nan
    half = stat / 2.0
    m = df // 2
    term = 1.0
    total = 1.0
    for j in range(1, m):
        term *= half / j
        total += term
    sf = math.exp(-half) * total
    return min(1.0, max(0.0, sf))


def bh_fdr(p_values: list[float]) -> list[float]:
    out = [math.nan] * len(p_values)
    valid = [(p, i) for i, p in enumerate(p_values) if math.isfinite(p)]
    if not valid:
        return out
    valid.sort(key=lambda item: item[0])
    n = len(valid)
    adjusted = [0.0] * n
    running = 1.0
    for rank_from_zero in range(n - 1, -1, -1):
        p, _ = valid[rank_from_zero]
        rank = rank_from_zero + 1
        running = min(running, p * n / rank)
        adjusted[rank_from_zero] = min(1.0, max(0.0, running))
    for value, (_, original_idx) in zip(adjusted, valid):
        out[original_idx] = value
    return out


def _random_effect_meta_one(effects: list[float], ses: list[float], valid: list[bool]) -> dict[str, float]:
    weights = [1.0 / (se * se) if ok and se > 0 else 0.0 for ok, se in zip(valid, ses)]
    sum_w = sum(weights)
    k = sum(1 for ok in valid if ok)
    if k < 1 or sum_w <= 0:
        return {
            "fixed_effect": math.nan,
            "fixed_se": math.nan,
            "fixed_p": math.nan,
            "random_effect": math.nan,
            "random_se": math.nan,
            "random_p": math.nan,
            "Q": math.nan,
            "I2": math.nan,
            "tau2": math.nan,
            "k": k,
        }
    fixed_effect = sum(w * e for w, e in zip(weights, effects) if math.isfinite(e)) / sum_w
    fixed_se = math.sqrt(1.0 / sum_w)
    fixed_p = _normal_two_sided_p(fixed_effect / fixed_se)
    q_stat = sum(w * ((e - fixed_effect) ** 2) for w, e, ok in zip(weights, effects, valid) if ok)
    df = k - 1
    c_term = sum_w - (sum(w * w for w in weights) / sum_w)
    tau2 = (q_stat - df) / c_term if c_term > 0 and q_stat > df else 0.0
    if not math.isfinite(tau2) or tau2 < 0:
        tau2 = 0.0
    re_weights = [
        1.0 / ((se * se) + tau2) if ok and se > 0 else 0.0
        for ok, se in zip(valid, ses)
    ]
    sum_re_w = sum(re_weights)
    random_effect = (
        sum(w * e for w, e in zip(re_weights, effects) if math.isfinite(e)) / sum_re_w
        if sum_re_w > 0
        else math.nan
    )
    random_se = math.sqrt(1.0 / sum_re_w) if sum_re_w > 0 else math.nan
    random_p = _normal_two_sided_p(random_effect / random_se)
    i2 = max(0.0, (q_stat - df) / q_stat) * 100.0 if q_stat > 0 and k >= 2 else 0.0
    if k < 2:
        q_stat = math.nan
        i2 = math.nan
        tau2 = math.nan
    return {
        "fixed_effect": fixed_effect,
        "fixed_se": fixed_se,
        "fixed_p": fixed_p,
        "random_effect": random_effect,
        "random_se": random_se,
        "random_p": random_p,
        "Q": q_stat,
        "I2": i2,
        "tau2": tau2,
        "k": k,
    }


def _fisher_partial_conjunction(p_values: list[float], valid: list[bool], r: int) -> float:
    vals = sorted(max(1e-300, min(1.0, p)) for p, ok in zip(p_values, valid) if ok and math.isfinite(p))
    k = len(vals)
    if k < r or r < 1:
        return math.nan
    tail = vals[r - 1 :]
    stat = -2.0 * sum(math.log(p) for p in tail)
    return _chi2_sf_even(stat, 2 * len(tail))


def _directional_pc_p(effects: list[float], p_values: list[float], valid: list[bool], r: int) -> float:
    p_up = []
    p_down = []
    for effect, p in zip(effects, p_values):
        if not math.isfinite(effect) or not math.isfinite(p):
            p_up.append(math.nan)
            p_down.append(math.nan)
            continue
        p_up.append(p / 2.0 if effect >= 0 else 1.0 - p / 2.0)
        p_down.append(p / 2.0 if effect <= 0 else 1.0 - p / 2.0)
    pc_up = _fisher_partial_conjunction(p_up, valid, r)
    pc_down = _fisher_partial_conjunction(p_down, valid, r)
    candidates = [p for p in (pc_up, pc_down) if math.isfinite(p)]
    return min(candidates) if candidates else math.nan


def _read_branch_records(
    cohorts: list[MetaCohort],
    branch: str,
    filename: str,
    allow_missing_branches: bool,
) -> tuple[dict[str, dict[str, object]], list[str]]:
    records: dict[str, dict[str, object]] = {}
    warnings: list[str] = []
    n = len(cohorts)
    for cohort_idx, cohort in enumerate(cohorts):
        table_path = cohort.result_dir / filename
        if not table_path.exists():
            message = f"{branch}: missing {filename} in {cohort.result_dir}"
            if allow_missing_branches:
                warnings.append(message)
                continue
            raise FileNotFoundError(message)
        delimiter = _sniff_delimiter(table_path)
        with _open_text(table_path) as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if not reader.fieldnames:
                raise ValueError(f"No header found in {table_path}")
            missing_required = [col for col in ("CpG", "logFC", "P.Value") if col not in reader.fieldnames]
            if missing_required:
                raise ValueError(f"{table_path} missing required columns: {', '.join(missing_required)}")
            if "SE" not in reader.fieldnames and "t" not in reader.fieldnames:
                raise ValueError(f"{table_path} must include either SE or t for meta-analysis")
            seen_cpgs: set[str] = set()
            for row in reader:
                cpg = str(row.get("CpG") or "").strip()
                if not cpg:
                    continue
                if cpg in seen_cpgs:
                    continue
                seen_cpgs.add(cpg)
                rec = records.get(cpg)
                if rec is None:
                    rec = {
                        "CpG": cpg,
                        "annotations": {col: "" for col in ANNOTATION_COLUMNS},
                        "effects": [math.nan] * n,
                        "ses": [math.nan] * n,
                        "p_values": [math.nan] * n,
                        "deltas": [math.nan] * n,
                    }
                    records[cpg] = rec
                annotations = rec["annotations"]
                assert isinstance(annotations, dict)
                for col in ANNOTATION_COLUMNS:
                    if not annotations.get(col) and not _is_missing(row.get(col)):
                        annotations[col] = str(row.get(col)).strip()
                effect = _safe_float(row.get("logFC"))
                se = _safe_float(row.get("SE"))
                if not math.isfinite(se):
                    t_stat = _safe_float(row.get("t"))
                    if math.isfinite(effect) and math.isfinite(t_stat) and t_stat != 0:
                        se = abs(effect / t_stat)
                if not math.isfinite(se) or se <= 0:
                    se = math.nan
                rec["effects"][cohort_idx] = effect
                rec["ses"][cohort_idx] = se
                rec["p_values"][cohort_idx] = _safe_float(row.get("P.Value"))
                rec["deltas"][cohort_idx] = _safe_float(row.get("Delta_Beta"))
    return records, warnings


def _pooled_delta(deltas: list[float], valid: list[bool], weights: list[float]) -> float:
    total_w = 0.0
    total = 0.0
    for delta, ok, weight in zip(deltas, valid, weights):
        if ok and math.isfinite(delta):
            total_w += weight
            total += weight * delta
    return total / total_w if total_w > 0 else math.nan


def _loo_direction_fraction(effects: list[float], ses: list[float], valid: list[bool], full_sign: int, min_cohorts: int) -> tuple[int, int, float]:
    same = 0
    possible = 0
    if full_sign == 0:
        return same, possible, 0.0
    for idx in range(len(valid)):
        if not valid[idx]:
            continue
        loo_valid = list(valid)
        loo_valid[idx] = False
        if sum(1 for ok in loo_valid if ok) < min_cohorts:
            continue
        possible += 1
        meta = _random_effect_meta_one(effects, ses, loo_valid)
        effect = meta["random_effect"]
        sign = 1 if effect > 0 else -1 if effect < 0 else 0
        if math.isfinite(effect) and sign == full_sign:
            same += 1
    return same, possible, (same / possible if possible > 0 else 0.0)


def _analyze_branch(
    cohorts: list[MetaCohort],
    branch: str,
    filename: str,
    thresholds: MetaThresholds,
    allow_missing_branches: bool,
) -> tuple[list[dict[str, object]], dict[str, object], list[str]]:
    records, warnings = _read_branch_records(cohorts, branch, filename, allow_missing_branches)
    sample_weights = [cohort.sample_weight for cohort in cohorts]
    rows: list[dict[str, object]] = []
    random_p_values: list[float] = []
    fixed_p_values: list[float] = []
    pc_fisher_values: list[float] = []
    pc_directional_values: list[float] = []

    for rec in records.values():
        effects = list(rec["effects"])
        ses = list(rec["ses"])
        p_values = list(rec["p_values"])
        deltas = list(rec["deltas"])
        valid = [
            math.isfinite(effect) and math.isfinite(se) and se > 0
            for effect, se in zip(effects, ses)
        ]
        meta = _random_effect_meta_one(effects, ses, valid)
        k = int(meta["k"])
        enough = k >= thresholds.min_cohorts
        n_up = sum(1 for effect, ok in zip(effects, valid) if ok and effect > 0)
        n_down = sum(1 for effect, ok in zip(effects, valid) if ok and effect < 0)
        majority = max(n_up, n_down)
        direction_fraction = majority / k if k > 0 else 0.0
        direction = "up" if n_up > n_down else "down" if n_down > n_up else "tie"
        full_sign = 1 if meta["random_effect"] > 0 else -1 if meta["random_effect"] < 0 else 0
        loo_same, loo_valid_count, loo_fraction = _loo_direction_fraction(
            effects, ses, valid, full_sign, thresholds.min_cohorts
        )
        pc_fisher = _fisher_partial_conjunction(p_values, valid, thresholds.partial_conjunction_r)
        pc_directional = _directional_pc_p(effects, p_values, valid, thresholds.partial_conjunction_r)
        row = {
            "CpG": rec["CpG"],
            **rec["annotations"],
            "branch": branch,
            "n_valid_cohorts": k,
            "fixed_effect_logFC": meta["fixed_effect"],
            "fixed_se": meta["fixed_se"],
            "fixed_p": meta["fixed_p"] if enough else math.nan,
            "random_effect_logFC": meta["random_effect"],
            "random_se": meta["random_se"],
            "random_p": meta["random_p"] if enough else math.nan,
            "Q": meta["Q"],
            "I2": meta["I2"],
            "tau2": meta["tau2"],
            "direction": direction,
            "direction_fraction": direction_fraction,
            "n_up": n_up,
            "n_down": n_down,
            "pooled_delta_beta_weighted": _pooled_delta(deltas, valid, sample_weights),
            "loo_same_direction_count": loo_same,
            "loo_valid_count": loo_valid_count,
            "loo_direction_fraction": loo_fraction,
            "pc_fisher_p": pc_fisher if enough else math.nan,
            "pc_directional_p": pc_directional if enough else math.nan,
        }
        for cohort, effect, p_val, delta in zip(cohorts, effects, p_values, deltas):
            safe_id = cohort.cohort_id.replace(" ", "_")
            row[f"effect_{safe_id}"] = effect
            row[f"p_{safe_id}"] = p_val
            row[f"delta_{safe_id}"] = delta
        rows.append(row)
        fixed_p_values.append(row["fixed_p"])
        random_p_values.append(row["random_p"])
        pc_fisher_values.append(row["pc_fisher_p"])
        pc_directional_values.append(row["pc_directional_p"])

    fixed_fdr = bh_fdr(fixed_p_values)
    random_fdr = bh_fdr(random_p_values)
    pc_fisher_fdr = bh_fdr(pc_fisher_values)
    pc_directional_fdr = bh_fdr(pc_directional_values)
    for row, fixed_adj, random_adj, pc_fisher_adj, pc_directional_adj in zip(
        rows, fixed_fdr, random_fdr, pc_fisher_fdr, pc_directional_fdr
    ):
        row["fixed_fdr"] = fixed_adj
        row["random_fdr"] = random_adj
        row["pc_fisher_fdr"] = pc_fisher_adj
        row["pc_directional_fdr"] = pc_directional_adj
        row["core_candidate"] = (
            int(row["n_valid_cohorts"]) >= thresholds.min_cohorts
            and math.isfinite(row["random_fdr"])
            and row["random_fdr"] < thresholds.meta_fdr
            and row["direction_fraction"] >= thresholds.min_direction_fraction
            and row["loo_valid_count"] > 0
            and row["loo_direction_fraction"] >= thresholds.min_loo_direction_fraction
            and math.isfinite(row["I2"])
            and row["I2"] <= thresholds.max_i2
            and math.isfinite(row["pooled_delta_beta_weighted"])
            and abs(row["pooled_delta_beta_weighted"]) >= thresholds.min_abs_delta_beta
        )
        row["replicable_pc_candidate"] = (
            int(row["n_valid_cohorts"]) >= thresholds.min_cohorts
            and math.isfinite(row["pc_directional_fdr"])
            and row["pc_directional_fdr"] < thresholds.meta_fdr
            and row["direction_fraction"] >= thresholds.min_direction_fraction
        )

    analyzable_i2 = [
        row["I2"] for row in rows
        if int(row["n_valid_cohorts"]) >= thresholds.min_cohorts and math.isfinite(row["I2"])
    ]
    core_loo = [row["loo_direction_fraction"] for row in rows if row["core_candidate"]]
    summary = {
        "branch": branch,
        "input_table": filename,
        "n_total_cpgs": len(rows),
        "n_meta_analyzable_cpgs": sum(
            1 for row in rows if int(row["n_valid_cohorts"]) >= thresholds.min_cohorts
        ),
        "n_random_fdr_lt_0_05": sum(
            1 for row in rows if math.isfinite(row["random_fdr"]) and row["random_fdr"] < thresholds.meta_fdr
        ),
        "n_core_candidates": sum(1 for row in rows if row["core_candidate"]),
        "n_pc_directional_candidates": sum(1 for row in rows if row["replicable_pc_candidate"]),
        "median_i2_analyzable": _median(analyzable_i2),
        "median_loo_direction_fraction_core": _median(core_loo),
    }
    return rows, summary, warnings


def _median(values: list[float]) -> float:
    vals = sorted(v for v in values if math.isfinite(v))
    if not vals:
        return math.nan
    n = len(vals)
    mid = n // 2
    if n % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2.0


def _sort_key(row: dict[str, object]) -> tuple[float, float, float]:
    def finite_or_inf(value: object) -> float:
        val = float(value) if isinstance(value, (float, int)) else math.nan
        return val if math.isfinite(val) else math.inf
    return (
        finite_or_inf(row.get("random_fdr")),
        finite_or_inf(row.get("I2")),
        finite_or_inf(row.get("random_p")),
    )


def _write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str], gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gzip_output else open
    mode = "wt" if gzip_output else "w"
    with opener(path, mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format_float(row.get(key)) for key in fieldnames})


def _summarize_branch_concordance(candidate_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in candidate_rows:
        grouped.setdefault(str(row["CpG"]), []).append(row)
    out: list[dict[str, object]] = []
    for cpg, group in grouped.items():
        directions = sorted(set(str(row.get("direction") or "") for row in group))
        direction = directions[0] if len(directions) == 1 else "mixed"
        branches = sorted(set(str(row.get("branch") or "") for row in group))
        out.append(
            {
                "CpG": cpg,
                "Gene": _first_nonempty(row.get("Gene") for row in group),
                "chr": _first_nonempty(row.get("chr") for row in group),
                "pos": _first_nonempty(row.get("pos") for row in group),
                "Region": _first_nonempty(row.get("Region") for row in group),
                "Island_Context": _first_nonempty(row.get("Island_Context") for row in group),
                "n_candidate_branches": len(branches),
                "branches": ";".join(branches),
                "direction": direction,
                "mean_random_effect_logFC": _mean(row.get("random_effect_logFC") for row in group),
                "mean_pooled_delta_beta": _mean(row.get("pooled_delta_beta_weighted") for row in group),
                "max_random_fdr": _max(row.get("random_fdr") for row in group),
                "max_I2": _max(row.get("I2") for row in group),
                "min_direction_fraction": _min(row.get("direction_fraction") for row in group),
                "min_loo_direction_fraction": _min(row.get("loo_direction_fraction") for row in group),
            }
        )
    out.sort(key=lambda row: (-int(row["n_candidate_branches"]), row["max_random_fdr"], row["max_I2"]))
    return out


def _mean(values: Iterable[object]) -> float:
    vals = [_safe_float(v) for v in values]
    vals = [v for v in vals if math.isfinite(v)]
    return sum(vals) / len(vals) if vals else math.nan


def _max(values: Iterable[object]) -> float:
    vals = [_safe_float(v) for v in values]
    vals = [v for v in vals if math.isfinite(v)]
    return max(vals) if vals else math.nan


def _min(values: Iterable[object]) -> float:
    vals = [_safe_float(v) for v in values]
    vals = [v for v in vals if math.isfinite(v)]
    return min(vals) if vals else math.nan


def _base_fieldnames(cohorts: list[MetaCohort]) -> list[str]:
    fields = [
        "CpG",
        "Gene",
        "chr",
        "pos",
        "Region",
        "Island_Context",
        "branch",
        "n_valid_cohorts",
        "fixed_effect_logFC",
        "fixed_se",
        "fixed_p",
        "random_effect_logFC",
        "random_se",
        "random_p",
        "Q",
        "I2",
        "tau2",
        "direction",
        "direction_fraction",
        "n_up",
        "n_down",
        "pooled_delta_beta_weighted",
        "loo_same_direction_count",
        "loo_valid_count",
        "loo_direction_fraction",
        "pc_fisher_p",
        "pc_directional_p",
        "fixed_fdr",
        "random_fdr",
        "pc_fisher_fdr",
        "pc_directional_fdr",
        "core_candidate",
        "replicable_pc_candidate",
    ]
    for cohort in cohorts:
        safe_id = cohort.cohort_id.replace(" ", "_")
        fields.extend([f"effect_{safe_id}", f"p_{safe_id}", f"delta_{safe_id}"])
    return fields


def _write_report(
    out_dir: Path,
    cohorts: list[MetaCohort],
    branch_summaries: list[dict[str, object]],
    concordant: list[dict[str, object]],
    thresholds: MetaThresholds,
    warnings: list[str],
) -> None:
    same_ge2 = [
        row for row in concordant
        if int(row["n_candidate_branches"]) >= 2 and row["direction"] != "mixed"
    ]
    lines = [
        "# IlluMeta Cross-cohort CpG Meta-analysis Report",
        "",
        "This run combines completed IlluMeta branch-level DMP outputs across cohorts.",
        "Branches are analyzed separately; cross-branch concordance is a preprocessing robustness filter, not an independent-study meta-analysis.",
        "",
        "## Inputs",
        "",
        "| Cohort | Result directory | N control | N test | Primary mode |",
        "| --- | --- | ---: | ---: | --- |",
    ]
    for cohort in cohorts:
        lines.append(
            f"| {cohort.cohort_id} | `{cohort.result_dir}` | {cohort.n_con} | {cohort.n_test} | {cohort.primary_result_mode or 'NA'} |"
        )
    lines.extend(
        [
            "",
            "## Candidate Rule",
            "",
            f"- Random-effects FDR < {thresholds.meta_fdr:g}",
            f"- At least {thresholds.min_cohorts} valid cohorts",
            f"- Direction fraction >= {thresholds.min_direction_fraction:g}",
            f"- Leave-one-cohort-out direction fraction >= {thresholds.min_loo_direction_fraction:g}",
            f"- I2 <= {thresholds.max_i2:g}%",
            f"- Absolute sample-size-weighted pooled delta beta >= {thresholds.min_abs_delta_beta:g}",
            f"- Directional partial-conjunction screen uses r={thresholds.partial_conjunction_r}",
            "",
            "## Branch Summary",
            "",
            "| Branch | Analyzable CpGs | Random FDR hits | Core candidates | Directional PC candidates | Median I2 | Median core LOO direction |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for summary in branch_summaries:
        lines.append(
            "| {branch} | {n_meta_analyzable_cpgs:,} | {n_random_fdr_lt_0_05:,} | "
            "{n_core_candidates:,} | {n_pc_directional_candidates:,} | {median_i2_analyzable:.2f} | "
            "{median_loo_direction_fraction_core:.2f} |".format(**summary)
        )
    lines.extend(
        [
            "",
            "## Branch-concordant Core Candidates",
            "",
            f"Core candidates present in at least two branches with the same direction: {len(same_ge2):,}",
            "",
        ]
    )
    if same_ge2:
        lines.extend(
            [
                "| CpG | Gene | Branches | Direction | Mean delta beta | Max FDR | Max I2 |",
                "| --- | --- | --- | --- | ---: | ---: | ---: |",
            ]
        )
        for row in same_ge2[:20]:
            lines.append(
                f"| {row['CpG']} | {row.get('Gene', '')} | {row['branches']} | {row['direction']} | "
                f"{_format_float(row['mean_pooled_delta_beta'])} | {_format_float(row['max_random_fdr'])} | "
                f"{_format_float(row['max_I2'])} |"
            )
    if warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in warnings:
            lines.append(f"- {warning}")
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "- Treat this as effect-size meta-analysis and robustness prioritization, not a per-cohort FDR intersection.",
            "- Branch concordance checks preprocessing robustness because branches reuse the same samples.",
            "- High-I2 partial-conjunction hits should be reported as heterogeneous directional replication, not universal markers.",
        ]
    )
    (out_dir / "meta_analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    methods = [
        "# IlluMeta Meta-analysis Methods",
        "",
        "IlluMeta cross-cohort meta-analysis used branch-level DMP tables from completed IlluMeta runs.",
        "For each preprocessing branch, cohort-level logFC estimates and standard errors were combined per CpG.",
        "When an SE column was absent, SE was reconstructed as abs(logFC / t) from the limma moderated t-statistic.",
        "Fixed-effect estimates used inverse-variance weights. Random-effects estimates used a DerSimonian-Laird tau2 estimator.",
        "Benjamini-Hochberg FDR was applied within each branch. Directional consistency, I2, sample-size-weighted delta beta, and leave-one-cohort-out directional stability were used as robustness filters.",
        "Minfi, Sesame strict, and Sesame native branches were not pooled as independent cohorts; branch concordance was used only as a preprocessing sensitivity criterion.",
        "",
    ]
    (out_dir / "meta_methods.md").write_text("\n".join(methods), encoding="utf-8")


def run_meta_analysis(
    cohorts: list[MetaCohort],
    branch_files: dict[str, str],
    out_dir: Path,
    thresholds: MetaThresholds,
    allow_missing_branches: bool = False,
) -> dict[str, object]:
    if len(cohorts) < thresholds.min_cohorts:
        raise ValueError(
            f"Need at least {thresholds.min_cohorts} cohorts for this meta-analysis; got {len(cohorts)}"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    all_summaries: list[dict[str, object]] = []
    all_candidate_rows: list[dict[str, object]] = []
    warnings: list[str] = []
    fieldnames = _base_fieldnames(cohorts)

    for branch, filename in branch_files.items():
        _log(f"Running {branch} meta-analysis from {filename}...")
        rows, summary, branch_warnings = _analyze_branch(
            cohorts, branch, filename, thresholds, allow_missing_branches
        )
        warnings.extend(branch_warnings)
        rows.sort(key=_sort_key)
        all_summaries.append(summary)
        _write_tsv(out_dir / f"{branch}_meta_full.tsv.gz", rows, fieldnames, gzip_output=True)
        _write_tsv(out_dir / f"{branch}_meta_top{thresholds.top_n}.tsv", rows[: thresholds.top_n], fieldnames)
        core = [row for row in rows if row["core_candidate"]]
        pc_hits = [row for row in rows if row["replicable_pc_candidate"]]
        _write_tsv(out_dir / f"{branch}_core_candidates.tsv", core, fieldnames)
        _write_tsv(out_dir / f"{branch}_pc_candidates.tsv", pc_hits, fieldnames)
        all_candidate_rows.extend(core)
        _log(f"  {branch}: {summary['n_core_candidates']} core candidates.")

    summary_fields = [
        "branch",
        "input_table",
        "n_total_cpgs",
        "n_meta_analyzable_cpgs",
        "n_random_fdr_lt_0_05",
        "n_core_candidates",
        "n_pc_directional_candidates",
        "median_i2_analyzable",
        "median_loo_direction_fraction_core",
    ]
    _write_tsv(out_dir / "meta_branch_summary.tsv", all_summaries, summary_fields)

    concordant = _summarize_branch_concordance(all_candidate_rows)
    concordance_fields = [
        "CpG",
        "Gene",
        "chr",
        "pos",
        "Region",
        "Island_Context",
        "n_candidate_branches",
        "branches",
        "direction",
        "mean_random_effect_logFC",
        "mean_pooled_delta_beta",
        "max_random_fdr",
        "max_I2",
        "min_direction_fraction",
        "min_loo_direction_fraction",
    ]
    _write_tsv(out_dir / "branch_concordant_core_candidates.tsv", concordant, concordance_fields)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "cohorts": [
            {
                "cohort_id": cohort.cohort_id,
                "result_dir": str(cohort.result_dir),
                "label": cohort.label,
                "platform": cohort.platform,
                "tissue": cohort.tissue,
                "n_con": cohort.n_con,
                "n_test": cohort.n_test,
                "primary_result_mode": cohort.primary_result_mode,
                "primary_lambda_guard_status": cohort.primary_lambda_guard_status,
            }
            for cohort in cohorts
        ],
        "branches": branch_files,
        "thresholds": thresholds.__dict__,
        "branch_summaries": all_summaries,
        "warnings": warnings,
    }
    (out_dir / "meta_input_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    _write_report(out_dir, cohorts, all_summaries, concordant, thresholds, warnings)
    return manifest


def run_meta_cli(args) -> int:
    project_root = Path(args.project_root or os.getcwd()).expanduser().resolve()
    out_dir = _resolve_path(args.output, project_root)
    branch_files = parse_branch_selection(args.branches, args.branch_file)
    thresholds = MetaThresholds(
        min_cohorts=args.min_cohorts,
        meta_fdr=args.meta_fdr,
        min_direction_fraction=args.min_direction_fraction,
        min_loo_direction_fraction=args.min_loo_direction_fraction,
        max_i2=args.max_i2,
        min_abs_delta_beta=args.min_abs_delta_beta,
        partial_conjunction_r=args.partial_conjunction_r,
        top_n=args.top_n,
    )
    if thresholds.min_cohorts < 2:
        raise ValueError("--min-cohorts must be at least 2")
    if thresholds.partial_conjunction_r < 1:
        raise ValueError("--partial-conjunction-r must be at least 1")
    cohorts: list[MetaCohort] = []
    if args.manifest:
        cohorts.extend(load_meta_manifest(_resolve_path(args.manifest, project_root), project_root))
    if args.result_dirs:
        cohorts.extend(load_positional_cohorts(args.result_dirs, project_root))
    if not cohorts:
        raise ValueError("Provide result directories or --manifest.")
    seen: set[str] = set()
    unique_cohorts: list[MetaCohort] = []
    for cohort in cohorts:
        key = str(cohort.result_dir)
        if key in seen:
            continue
        seen.add(key)
        if not cohort.result_dir.exists():
            raise FileNotFoundError(f"Result directory not found: {cohort.result_dir}")
        unique_cohorts.append(cohort)
    _log(f"Starting cross-cohort meta-analysis: {len(unique_cohorts)} cohorts, {len(branch_files)} branch(es).")
    run_meta_analysis(
        unique_cohorts,
        branch_files,
        out_dir,
        thresholds,
        allow_missing_branches=args.allow_missing_branches,
    )
    _log(f"Meta-analysis complete: {out_dir}")
    return 0
