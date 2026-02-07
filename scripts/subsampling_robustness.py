#!/usr/bin/env python3
"""
Subsampling / imbalance robustness benchmark for IlluMeta.

Goal:
- Repeatedly subsample a "large" dataset under balanced and unbalanced designs
  and quantify stability + correction adequacy, producing one figure + one table.

Typical workflow:
1) Plan runs (creates subset configure.tsv files + jobs TSV):
   python3 scripts/subsampling_robustness.py plan \
     --config projects/GSEXXXX/configure.tsv \
     --idat-dir projects/GSEXXXX/idat \
     --group-con Control --group-test Case \
     --out-root benchmarks/subsampling_robustness/GSEXXXX

2) Run them (sequential; logs per run):
   python3 scripts/subsampling_robustness.py run \
     --jobs benchmarks/subsampling_robustness/GSEXXXX/jobs.tsv

3) Summarize (writes table + multi-panel figure):
   python3 scripts/subsampling_robustness.py summarize \
     --plan benchmarks/subsampling_robustness/GSEXXXX/plan.tsv \
     --out-dir benchmarks/paper_figures
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import re
import shlex
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


REQUIRED_JOB_FIELDS = ("config", "group_con", "group_test", "name", "output")


def _sniff_delimiter(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            head = fh.readline()
    except OSError:
        return "\t"
    if "\t" in head:
        return "\t"
    if "," in head:
        return ","
    return "\t"


def read_table(path: Path) -> Tuple[List[Dict[str, str]], List[str], str]:
    delim = _sniff_delimiter(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        rows = [{k: (v or "") for k, v in row.items()} for row in reader]
    return rows, (reader.fieldnames or []), delim


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def normalize_text(s: str) -> str:
    return (s or "").strip().lower()


def detect_id_col(headers: Sequence[str], override: str = "") -> str:
    if override:
        if override not in headers:
            raise ValueError(f"--id-column '{override}' not found in configure.tsv")
        return override
    for col in headers:
        if re.search(r"(GSM|geo_accession)", col, flags=re.IGNORECASE):
            return col
    if "Basename" in headers:
        return "Basename"
    if headers:
        return headers[0]
    raise ValueError("configure.tsv missing headers")


def safe_int_round(x: float) -> int:
    if not math.isfinite(x):
        return 0
    return int(round(x))


def clamp(n: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, n))


@dataclass(frozen=True)
class Scenario:
    scenario: str
    frac_con: float
    frac_test: float
    rep: int
    seed: int

    @property
    def name(self) -> str:
        # Stable filesystem-friendly run name
        def pct(x: float) -> str:
            return f"{int(round(100 * x)):02d}"

        return f"{self.scenario}_c{pct(self.frac_con)}_t{pct(self.frac_test)}_rep{self.rep:02d}"


def build_scenarios(
    balanced_fracs: Sequence[float],
    balanced_reps: int,
    unbalanced_pairs: Sequence[Tuple[float, float]],
    unbalanced_reps: int,
    include_full: bool,
    full_reps: int,
    seed: int,
) -> List[Scenario]:
    out: List[Scenario] = []
    cur_seed = seed

    if include_full:
        for rep in range(1, full_reps + 1):
            out.append(Scenario("full", 1.0, 1.0, rep=rep, seed=cur_seed))
            cur_seed += 1

    for frac in balanced_fracs:
        for rep in range(1, balanced_reps + 1):
            out.append(Scenario("balanced", frac, frac, rep=rep, seed=cur_seed))
            cur_seed += 1

    for (fc, ft) in unbalanced_pairs:
        for rep in range(1, unbalanced_reps + 1):
            out.append(Scenario("unbalanced", fc, ft, rep=rep, seed=cur_seed))
            cur_seed += 1

    return out


def sample_rows(
    rows: Sequence[Dict[str, str]],
    headers: Sequence[str],
    id_col: str,
    group_col: str,
    group_con: str,
    group_test: str,
    frac_con: float,
    frac_test: float,
    min_per_group: int,
    seed: int,
) -> Tuple[List[Dict[str, str]], Dict[str, int], List[Tuple[str, str]]]:
    con_norm = normalize_text(group_con)
    test_norm = normalize_text(group_test)

    con_rows = [r for r in rows if normalize_text(r.get(group_col, "")) == con_norm]
    test_rows = [r for r in rows if normalize_text(r.get(group_col, "")) == test_norm]

    if not con_rows or not test_rows:
        raise ValueError(
            f"Could not find both groups in configure.tsv: control='{group_con}' (n={len(con_rows)}), "
            f"test='{group_test}' (n={len(test_rows)})."
        )

    if len(con_rows) < min_per_group or len(test_rows) < min_per_group:
        raise ValueError(
            f"Not enough samples to subsample with min_per_group={min_per_group}: "
            f"control n={len(con_rows)}, test n={len(test_rows)}."
        )

    rng = random.Random(seed)

    n_con = clamp(safe_int_round(frac_con * len(con_rows)), min_per_group, len(con_rows))
    n_test = clamp(safe_int_round(frac_test * len(test_rows)), min_per_group, len(test_rows))
    con_sel = rng.sample(con_rows, n_con) if n_con < len(con_rows) else list(con_rows)
    test_sel = rng.sample(test_rows, n_test) if n_test < len(test_rows) else list(test_rows)

    subset = list(con_sel) + list(test_sel)

    # Record selected IDs (for reproducibility and debugging)
    chosen: List[Tuple[str, str]] = []
    for r in subset:
        chosen.append((r.get(id_col, "").strip(), r.get(group_col, "").strip()))

    counts = {"n_con": len(con_sel), "n_test": len(test_sel), "n_total": len(subset)}
    return subset, counts, chosen


def plan_cmd(args: argparse.Namespace) -> int:
    base_config = Path(args.config).expanduser().resolve()
    if not base_config.exists():
        raise SystemExit(f"Missing --config: {base_config}")

    idat_dir = Path(args.idat_dir).expanduser().resolve()
    if not idat_dir.exists():
        raise SystemExit(f"Missing --idat-dir: {idat_dir}")

    out_root = Path(args.out_root).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    rows, headers, _ = read_table(base_config)
    if not headers:
        raise SystemExit(f"Empty config: {base_config}")
    if "primary_group" not in headers:
        raise SystemExit("configure.tsv missing required column: primary_group")

    id_col = detect_id_col(headers, override=args.id_column or "")
    group_col = "primary_group"

    unbalanced_pairs: List[Tuple[float, float]] = []
    for pair in args.unbalanced or []:
        unbalanced_pairs.append((pair[0], pair[1]))

    scenarios = build_scenarios(
        balanced_fracs=args.balanced_fracs,
        balanced_reps=args.balanced_reps,
        unbalanced_pairs=unbalanced_pairs,
        unbalanced_reps=args.unbalanced_reps,
        include_full=not args.no_full,
        full_reps=args.full_reps,
        seed=args.seed,
    )

    config_dir = out_root / "configs"
    jobs_path = out_root / "jobs.tsv"
    plan_path = out_root / "plan.tsv"

    plan_rows: List[Dict[str, str]] = []
    jobs_rows: List[Dict[str, str]] = []

    extra_args = (args.illumeta_extra_args or "").strip()
    # Enforce Tier3 non-fatal behavior for small subsamples.
    if "--tier3-on-fail" not in extra_args and "--tier3_on_fail" not in extra_args:
        extra_args = ("--tier3-on-fail skip " + extra_args).strip()

    # Always point to the same IDAT directory.
    extra_args = f"--idat-dir {shlex.quote(str(idat_dir))} {extra_args}".strip()

    for sc in scenarios:
        subset, counts, chosen = sample_rows(
            rows=rows,
            headers=headers,
            id_col=id_col,
            group_col=group_col,
            group_con=args.group_con,
            group_test=args.group_test,
            frac_con=sc.frac_con,
            frac_test=sc.frac_test,
            min_per_group=args.min_per_group,
            seed=sc.seed,
        )

        run_cfg_dir = config_dir / sc.name
        run_cfg_path = run_cfg_dir / "configure.tsv"
        run_samples_path = run_cfg_dir / "samples.tsv"
        write_tsv(run_cfg_path, headers, subset)
        write_tsv(run_samples_path, ["sample_id", "primary_group"], [{"sample_id": sid, "primary_group": grp} for sid, grp in chosen])

        out_dir = out_root / "runs" / sc.name

        plan_rows.append(
            {
                "name": sc.name,
                "scenario": sc.scenario,
                "frac_con": f"{sc.frac_con:.4g}",
                "frac_test": f"{sc.frac_test:.4g}",
                "rep": str(sc.rep),
                "seed": str(sc.seed),
                "n_con": str(counts["n_con"]),
                "n_test": str(counts["n_test"]),
                "n_total": str(counts["n_total"]),
                "config": str(run_cfg_path),
                "output": str(out_dir),
            }
        )

        jobs_rows.append(
            {
                "name": sc.name,
                "config": str(run_cfg_path),
                "group_con": args.group_con,
                "group_test": args.group_test,
                "output": str(out_dir),
                "extra_args": extra_args,
            }
        )

    write_tsv(plan_path, ["name", "scenario", "frac_con", "frac_test", "rep", "seed", "n_con", "n_test", "n_total", "config", "output"], plan_rows)
    write_tsv(jobs_path, ["config", "group_con", "group_test", "name", "output", "extra_args"], jobs_rows)

    print(f"Wrote plan: {plan_path}")
    print(f"Wrote jobs: {jobs_path}")
    return 0


def load_jobs(path: Path) -> List[Dict[str, str]]:
    rows, headers, _ = read_table(path)
    if not headers:
        raise ValueError("Job file is empty or missing headers.")
    for field in REQUIRED_JOB_FIELDS:
        if field not in headers:
            raise ValueError(f"Missing required column: {field}")
    jobs: List[Dict[str, str]] = []
    for row in rows:
        job = {k: (v or "").strip() for k, v in row.items()}
        if not job.get("config"):
            continue
        jobs.append(job)
    if not jobs:
        raise ValueError("No valid jobs found.")
    return jobs


def run_one_job(job: Dict[str, str], python_exe: str, illumeta_path: str, log_dir: Path, force: bool) -> Dict[str, str]:
    name = job["name"]
    config = job["config"]
    output = job["output"]
    extra_args = job.get("extra_args", "")

    out_dir = Path(output)
    out_dir.parent.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "summary.json"
    dashboard_path = Path(f"{output}_index.html")
    if summary_path.exists() and not force:
        return {
            "name": name,
            "output": output,
            "exit_code": "0",
            "status": "skipped",
            "duration_sec": "0.0",
            "summary_exists": "True",
            "dashboard_exists": str(dashboard_path.is_file()),
            "log_path": "",
        }

    cmd = [
        python_exe,
        illumeta_path,
        "analysis",
        "--config",
        config,
        "--group_con",
        job["group_con"],
        "--group_test",
        job["group_test"],
        "-o",
        output,
    ]
    if extra_args:
        cmd.extend(shlex.split(extra_args))

    env = os.environ.copy()
    env["ILLUMETA_TRACEBACK"] = "1"

    log_path = log_dir / f"{name}.log"
    start = time.time()
    with log_path.open("w", encoding="utf-8") as logf:
        proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, env=env)
    duration = time.time() - start

    return {
        "name": name,
        "output": output,
        "exit_code": str(proc.returncode),
        "status": "ok" if proc.returncode == 0 else "failed",
        "duration_sec": f"{duration:.1f}",
        "summary_exists": str(summary_path.is_file()),
        "dashboard_exists": str(dashboard_path.is_file()),
        "log_path": str(log_path),
    }


def run_cmd(args: argparse.Namespace) -> int:
    jobs_path = Path(args.jobs).expanduser().resolve()
    jobs = load_jobs(jobs_path)

    log_dir = Path(args.log_dir).expanduser().resolve() if args.log_dir else (jobs_path.parent / "logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    report_rows: List[Dict[str, str]] = []
    for job in jobs:
        row = run_one_job(
            job=job,
            python_exe=args.python,
            illumeta_path=args.illumeta,
            log_dir=log_dir,
            force=args.force,
        )
        report_rows.append(row)
        print(f"[{row['status']}] {row['name']} (exit={row['exit_code']}, sec={row['duration_sec']})")
        if args.stop_on_failure and row["exit_code"] != "0":
            break

    report_path = jobs_path.parent / "run_report.tsv"
    write_tsv(
        report_path,
        ["name", "output", "status", "exit_code", "duration_sec", "summary_exists", "dashboard_exists", "log_path"],
        report_rows,
    )
    print(f"Wrote run report: {report_path}")
    return 0


def _read_metric_kv_csv(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        out = {}
        for row in reader:
            k = (row.get("metric") or "").strip().strip('"')
            v = (row.get("value") or "").strip()
            if k:
                out[k] = v
        return out


def _to_float(v: str) -> Optional[float]:
    if v is None:
        return None
    s = str(v).strip().strip('"')
    if not s or s.lower() == "na":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _fmt_mean_sd(vals: List[float], digits: int = 3) -> str:
    if not vals:
        return "NA"
    if len(vals) == 1:
        return f"{vals[0]:.{digits}g}"
    m = statistics.mean(vals)
    sd = statistics.stdev(vals)
    return f"{m:.{digits}g}±{sd:.{digits}g}"


def _pearsonr(xs: List[float], ys: List[float]) -> Optional[float]:
    if len(xs) < 2 or len(ys) < 2 or len(xs) != len(ys):
        return None
    # statistics.correlation is available in 3.11, but compute explicitly for robustness.
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    num = 0.0
    dx = 0.0
    dy = 0.0
    for x, y in zip(xs, ys):
        num += (x - mx) * (y - my)
        dx += (x - mx) ** 2
        dy += (y - my) ** 2
    if dx <= 0 or dy <= 0:
        return None
    return num / math.sqrt(dx * dy)


def _load_consensus(path: Path, top_n: int) -> Tuple[List[str], Dict[str, float], int]:
    if not path.exists():
        return [], {}, 0
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = []
        for row in reader:
            cpg = (row.get("CpG") or "").strip().strip('"')
            if not cpg:
                continue
            adj = _to_float(row.get("adj.P.Val", ""))
            pval = _to_float(row.get("P.Value", ""))
            # prefer mean logFC if present
            lfc = _to_float(row.get("logFC_mean", "")) if "logFC_mean" in row else _to_float(row.get("logFC", ""))
            rows.append((cpg, adj if adj is not None else 1.0, pval if pval is not None else 1.0, lfc))
    total = len(rows)
    rows.sort(key=lambda t: (t[1], t[2], t[0]))
    top = rows[: max(0, top_n)]
    ids = [t[0] for t in top]
    lfc_map = {t[0]: (t[3] if t[3] is not None else float("nan")) for t in top}
    return ids, lfc_map, total


def summarize_cmd(args: argparse.Namespace) -> int:
    plan_path = Path(args.plan).expanduser().resolve()
    plan_rows, headers, _ = read_table(plan_path)
    if not plan_rows:
        raise SystemExit(f"Empty plan: {plan_path}")

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve reference run (full dataset) for overlap/correlation.
    ref_name = (args.reference or "").strip()
    if not ref_name:
        full = [r for r in plan_rows if (r.get("scenario") or "") == "full"]
        if full:
            ref_name = (full[0].get("name") or "").strip()
        else:
            # fallback: max n_total
            best = max(plan_rows, key=lambda r: int((r.get("n_total") or "0") or "0"))
            ref_name = (best.get("name") or "").strip()
    if not ref_name:
        raise SystemExit("Could not determine a reference run. Provide --reference.")

    by_name = {r.get("name", ""): r for r in plan_rows}
    if ref_name not in by_name:
        raise SystemExit(f"Reference run not found in plan.tsv: {ref_name}")

    top_n = int(args.top_n)

    # Load reference consensus list (Intersection Native) for comparisons.
    ref_out = Path(by_name[ref_name]["output"])
    ref_cons_path = ref_out / "Intersection_Native_Consensus_DMPs.csv"
    ref_ids, ref_lfc, ref_total = _load_consensus(ref_cons_path, top_n=top_n)
    ref_set = set(ref_ids)

    run_summ_rows: List[Dict[str, str]] = []
    scenario_groups: Dict[str, List[Dict[str, str]]] = {}

    for row in plan_rows:
        name = (row.get("name") or "").strip()
        scenario = (row.get("scenario") or "").strip()
        out = Path(row.get("output") or "")
        sj_path = out / "summary.json"
        cai_path = out / "Correction_Adequacy_Summary.csv"
        cons_path = out / "Intersection_Native_Consensus_DMPs.csv"

        status = "ok" if sj_path.exists() else "missing"
        sj = {}
        if sj_path.exists():
            try:
                sj = json.loads(sj_path.read_text(encoding="utf-8"))
            except Exception:
                status = "bad_summary"

        cai = _read_metric_kv_csv(cai_path)
        cai_f = {k: _to_float(v) for k, v in cai.items()}

        ids, lfc_map, total_cons = _load_consensus(cons_path, top_n=top_n)
        ids_set = set(ids)
        inter = ref_set & ids_set
        union = ref_set | ids_set
        jaccard = (len(inter) / len(union)) if union else None

        # logFC correlation on overlapping CpGs (top-N overlap)
        xs = []
        ys = []
        for cpg in sorted(inter):
            x = ref_lfc.get(cpg, float("nan"))
            y = lfc_map.get(cpg, float("nan"))
            if x is None or y is None:
                continue
            if not (math.isfinite(x) and math.isfinite(y)):
                continue
            xs.append(float(x))
            ys.append(float(y))
        logfc_r = _pearsonr(xs, ys) if xs and ys else None

        out_row = {
            "name": name,
            "scenario": scenario,
            "frac_con": row.get("frac_con", ""),
            "frac_test": row.get("frac_test", ""),
            "rep": row.get("rep", ""),
            "n_con": row.get("n_con", ""),
            "n_test": row.get("n_test", ""),
            "n_total": row.get("n_total", ""),
            "status": status,
            "crf_tier": str(sj.get("crf_sample_tier", "")),
            "primary_branch": str(sj.get("primary_branch", "")),
            "primary_mode": str(sj.get("primary_result_mode", "")),
            "primary_caf": str(sj.get("primary_caf_score", "")),
            "cai": str(cai_f.get("cai") if cai_f.get("cai") is not None else ""),
            "calibration_score": str(cai_f.get("calibration_score") if cai_f.get("calibration_score") is not None else ""),
            "preservation_score": str(cai_f.get("preservation_score") if cai_f.get("preservation_score") is not None else ""),
            "batch_removal_score": str(cai_f.get("batch_removal_score") if cai_f.get("batch_removal_score") is not None else ""),
            "null_lambda": str(cai_f.get("null_lambda") if cai_f.get("null_lambda") is not None else ""),
            "null_ks_p": str(cai_f.get("null_ks_p") if cai_f.get("null_ks_p") is not None else ""),
            "consensus_total": str(total_cons),
            "jaccard_topN_vs_ref": f"{jaccard:.6g}" if jaccard is not None else "",
            "logfc_r_topN_vs_ref": f"{logfc_r:.6g}" if logfc_r is not None else "",
            "output": str(out),
        }

        run_summ_rows.append(out_row)
        scenario_groups.setdefault(f"{scenario}:{row.get('frac_con','')}:{row.get('frac_test','')}", []).append(out_row)

    runs_path = out_dir / "subsampling_runs.tsv"
    write_tsv(
        runs_path,
        [
            "name",
            "scenario",
            "frac_con",
            "frac_test",
            "rep",
            "n_con",
            "n_test",
            "n_total",
            "status",
            "crf_tier",
            "primary_branch",
            "primary_mode",
            "primary_caf",
            "cai",
            "calibration_score",
            "preservation_score",
            "batch_removal_score",
            "null_lambda",
            "null_ks_p",
            "consensus_total",
            "jaccard_topN_vs_ref",
            "logfc_r_topN_vs_ref",
            "output",
        ],
        run_summ_rows,
    )

    # Aggregate per scenario bucket (mean±sd across reps).
    def _collect_float(rows_in: List[Dict[str, str]], key: str) -> List[float]:
        vals: List[float] = []
        for r in rows_in:
            if r.get("status") != "ok":
                continue
            v = _to_float(r.get(key, ""))
            if v is None or not math.isfinite(v):
                continue
            vals.append(float(v))
        return vals

    scenario_rows: List[Dict[str, str]] = []
    for bucket, items in scenario_groups.items():
        scenario, fc, ft = bucket.split(":", 2)
        n_total = [int(r.get("n_total") or "0") for r in items if r.get("status") == "ok"]
        n_con = [int(r.get("n_con") or "0") for r in items if r.get("status") == "ok"]
        n_test = [int(r.get("n_test") or "0") for r in items if r.get("status") == "ok"]

        scenario_rows.append(
            {
                "scenario": scenario,
                "frac_con": fc,
                "frac_test": ft,
                "n_runs": str(len(items)),
                "n_ok": str(sum(1 for r in items if r.get("status") == "ok")),
                "n_total": _fmt_mean_sd([float(x) for x in n_total], digits=4),
                "n_con": _fmt_mean_sd([float(x) for x in n_con], digits=4),
                "n_test": _fmt_mean_sd([float(x) for x in n_test], digits=4),
                "CAI": _fmt_mean_sd(_collect_float(items, "cai"), digits=4),
                "calibration": _fmt_mean_sd(_collect_float(items, "calibration_score"), digits=4),
                "preservation": _fmt_mean_sd(_collect_float(items, "preservation_score"), digits=4),
                "batch_removal": _fmt_mean_sd(_collect_float(items, "batch_removal_score"), digits=4),
                "null_lambda": _fmt_mean_sd(_collect_float(items, "null_lambda"), digits=4),
                f"jaccard_top{top_n}": _fmt_mean_sd(_collect_float(items, "jaccard_topN_vs_ref"), digits=4),
                f"logFC_r_top{top_n}": _fmt_mean_sd(_collect_float(items, "logfc_r_topN_vs_ref"), digits=4),
            }
        )

    # Stable order
    def _order_key(r: Dict[str, str]) -> Tuple[int, float, float]:
        scen = r.get("scenario", "")
        priority = {"full": 0, "balanced": 1, "unbalanced": 2}.get(scen, 9)
        fc = _to_float(r.get("frac_con", "")) or 0.0
        ft = _to_float(r.get("frac_test", "")) or 0.0
        return (priority, fc, ft)

    scenario_rows.sort(key=_order_key)

    table_path = out_dir / "Table_Subsampling_Robustness.tsv"
    table_fields = [
        "scenario",
        "frac_con",
        "frac_test",
        "n_ok",
        "n_runs",
        "n_total",
        "n_con",
        "n_test",
        "CAI",
        "calibration",
        "preservation",
        "batch_removal",
        "null_lambda",
        f"jaccard_top{top_n}",
        f"logFC_r_top{top_n}",
    ]
    write_tsv(table_path, table_fields, scenario_rows)

    # Write a simple Markdown table for quick copy/paste.
    md_path = out_dir / "Table_Subsampling_Robustness.md"
    col_widths = {k: max(len(k), max((len(r.get(k, "")) for r in scenario_rows), default=0)) for k in table_fields}
    with md_path.open("w", encoding="utf-8") as fh:
        fh.write("| " + " | ".join(table_fields) + " |\n")
        fh.write("| " + " | ".join("-" * col_widths[k] for k in table_fields) + " |\n")
        for r in scenario_rows:
            fh.write("| " + " | ".join((r.get(k, "") or "") for k in table_fields) + " |\n")

    print(f"Wrote run-level summary: {runs_path}")
    print(f"Wrote table: {table_path}")
    print(f"Wrote table (md): {md_path}")

    # Plot (optional; requires matplotlib). We keep plotting minimal and robust.
    fig_prefix = out_dir / "Figure_Subsampling_Robustness"
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as e:
        print(f"[!] matplotlib not available; skipping figure: {e}")
        return 0

    # Prepare per-run points (only OK)
    points = [r for r in run_summ_rows if r.get("status") == "ok"]
    # Build x-axis buckets consistent with scenario_rows order
    buckets = [(r["scenario"], r["frac_con"], r["frac_test"]) for r in scenario_rows]
    bucket_to_x = {b: i for i, b in enumerate(buckets)}
    x_labels = []
    for scen, fc, ft in buckets:
        if scen == "full":
            x_labels.append("Full\n(100/100)")
            continue
        x_labels.append(f"{scen}\n({int(round(100*float(fc))):d}/{int(round(100*float(ft))):d})")

    def jitter(x: float, seed: int) -> float:
        rng = random.Random(seed)
        return x + rng.uniform(-0.12, 0.12)

    metrics = [
        ("cai", "CAI (higher is better)", (0.0, 1.0), 1.0),
        ("null_lambda", "Null λ (ideal ≈ 1.0)", None, 1.0),
        ("jaccard_topN_vs_ref", f"Top-{top_n} consensus Jaccard vs ref", (0.0, 1.0), 1.0),
        ("logfc_r_topN_vs_ref", f"Top-{top_n} logFC correlation vs ref", (-1.0, 1.0), 1.0),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    axes = axes.flatten()

    for ax, (key, title, ylim, ref_line) in zip(axes, metrics):
        xs = []
        ys = []
        for p in points:
            b = (p["scenario"], p["frac_con"], p["frac_test"])
            if b not in bucket_to_x:
                continue
            v = _to_float(p.get(key, ""))
            if v is None or not math.isfinite(v):
                continue
            x = bucket_to_x[b]
            # deterministic jitter per run name
            seed_i = sum(ord(c) for c in p["name"]) % 100000
            xs.append(jitter(float(x), seed_i))
            ys.append(float(v))
        ax.scatter(xs, ys, s=28, alpha=0.8, color="#2563EB")

        # Mean±sd per bucket
        for b, x in bucket_to_x.items():
            vals = []
            for p in points:
                if (p["scenario"], p["frac_con"], p["frac_test"]) != b:
                    continue
                v = _to_float(p.get(key, ""))
                if v is None or not math.isfinite(v):
                    continue
                vals.append(float(v))
            if not vals:
                continue
            m = statistics.mean(vals)
            sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
            ax.errorbar([x], [m], yerr=[sd], fmt="o", color="#111827", capsize=4, markersize=5)

        if ref_line is not None:
            ax.axhline(ref_line, color="#16A34A", linestyle="--", linewidth=1.5, alpha=0.9)

        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xticks(list(range(len(x_labels))))
        ax.set_xticklabels(x_labels, rotation=0, fontsize=9)
        if ylim:
            ax.set_ylim(*ylim)
        ax.grid(axis="y", alpha=0.25)

    fig.suptitle("IlluMeta robustness under subsampling and group imbalance", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    fig.savefig(str(fig_prefix) + ".png", dpi=300)
    fig.savefig(str(fig_prefix) + ".pdf")
    plt.close(fig)

    print(f"Wrote figure: {fig_prefix}.png/.pdf")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Subsampling robustness benchmark (plan/run/summarize).")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("plan", help="Create subset configure.tsv files and a jobs TSV.")
    sp.add_argument("--config", required=True, help="Base configure.tsv for the full dataset.")
    sp.add_argument("--idat-dir", required=True, help="IDAT directory containing the full dataset IDAT pairs.")
    sp.add_argument("--group-con", required=True, help="Control group label (as in primary_group).")
    sp.add_argument("--group-test", required=True, help="Test group label (as in primary_group).")
    sp.add_argument("--out-root", required=True, help="Output root directory for configs/jobs/runs.")
    sp.add_argument("--id-column", default="", help="Optional sample ID column override (must exist in configure.tsv).")
    sp.add_argument("--min-per-group", type=int, default=6, help="Minimum samples per group for any planned run (default: 6).")
    sp.add_argument("--seed", type=int, default=42, help="Base seed for reproducible sampling (default: 42).")
    sp.add_argument("--balanced-fracs", type=float, nargs="+", default=[0.25, 0.50, 0.75], help="Balanced sampling fractions (default: 0.25 0.5 0.75).")
    sp.add_argument("--balanced-reps", type=int, default=3, help="Replicates per balanced fraction (default: 3).")
    sp.add_argument("--unbalanced", type=float, nargs=2, action="append", default=[[0.25, 0.75]], help="Unbalanced pair: frac_con frac_test (repeatable). Default: 0.25 0.75.")
    sp.add_argument("--unbalanced-reps", type=int, default=3, help="Replicates per unbalanced pair (default: 3).")
    sp.add_argument("--no-full", action="store_true", help="Do not include the full dataset run (100%/100%).")
    sp.add_argument("--full-reps", type=int, default=1, help="Replicates for full dataset run (default: 1).")
    sp.add_argument("--illumeta-extra-args", default="", help="Extra args appended to every illumeta.py analysis invocation.")

    sr = sub.add_parser("run", help="Run jobs TSV sequentially (writes logs and a run report).")
    sr.add_argument("--jobs", required=True, help="Jobs TSV created by `plan`.")
    sr.add_argument("--python", default=sys.executable, help="Python executable (default: current).")
    sr.add_argument("--illumeta", default="illumeta.py", help="Path to illumeta.py (default: illumeta.py).")
    sr.add_argument("--log-dir", default="", help="Log directory (default: sibling 'logs' next to jobs.tsv).")
    sr.add_argument("--force", action="store_true", help="Re-run even if summary.json already exists.")
    sr.add_argument("--stop-on-failure", action="store_true", help="Stop after first failed run.")

    ss = sub.add_parser("summarize", help="Summarize planned runs into one table + one figure.")
    ss.add_argument("--plan", required=True, help="Plan TSV created by `plan`.")
    ss.add_argument("--out-dir", required=True, help="Output directory for table/figure.")
    ss.add_argument("--reference", default="", help="Reference run name (default: first full run, else max n_total).")
    ss.add_argument("--top-n", type=int, default=1000, help="Top-N consensus DMPs for overlap/correlation (default: 1000).")

    return p


def main() -> int:
    args = build_parser().parse_args()
    if args.cmd == "plan":
        return plan_cmd(args)
    if args.cmd == "run":
        return run_cmd(args)
    if args.cmd == "summarize":
        return summarize_cmd(args)
    raise SystemExit("Unknown command")


if __name__ == "__main__":
    raise SystemExit(main())

