#!/usr/bin/env python3
"""Region-level meta-analysis across cohorts, calibrated by its own null.

Why this exists
---------------
Single-CpG cross-cohort meta-analysis on array data runs into a hard ceiling that is not
statistical but physical: the pooled effects that survive multiple-testing correction sit
at the array's reproducibility limit. Measured on this workflow's own output, the pooled
absolute beta difference of the calls a real analysis makes is indistinguishable from that
of the calls a sign-flipped null makes -- the probability that a randomly chosen real call
exceeds a randomly chosen false one is 0.49 to 0.52 across three preprocessing routes. So
effect size, at the single-CpG level, carries no information about which calls are real.
Counts separate; magnitudes do not.

Neighbouring CpGs are co-methylated, so a real regulatory difference tends to move a run of
them together while noise does not. Aggregating a run therefore improves the signal-to-noise
ratio in a way that thresholding a single probe harder cannot. This is the standard argument
for region-level testing (bumphunter, comb-p, DMRcate, dmrff), and dmrff is the method built
to work from effect estimates and standard errors, which is exactly what a meta-analysis
produces.

What is different here
----------------------
The same co-methylation that makes regions worth testing also makes their significance hard
to compute: the CpGs in a region are correlated, so combining them as if they were
independent is anti-conservative, and the size of that anti-conservatism depends on a
correlation structure that a summary-statistic meta-analysis does not have access to.
Methods that need the correlation matrix have to go back to the methylation matrices.

This module does not estimate the correlation. It calls regions with an explicit,
independence-assuming statistic -- stated as such, never presented as a calibrated P-value --
and then calibrates that statistic empirically by running the identical procedure on
sign-flipped replicates of the same cohorts. Flipping a whole cohort's sign preserves every
probe-probe correlation inside that cohort, so the null replicates carry the same
co-methylation structure as the observed data and the empirical null absorbs it without
anyone having to model it. What comes out is a region-level empirical P-value and a false
discovery count that were measured rather than assumed.

The cost is honest and worth stating: with k cohorts only 2^(k-1) - 1 non-trivial
whole-cohort sign patterns exist, so at k = 5 the null has 15 replicates and the resulting
tail resolution is coarse. Region P-values below about 1/15 are reported as bounded, not as
point estimates.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Sequence


DEFAULT_MAX_GAP = 500
DEFAULT_SEED_P = 0.05
DEFAULT_MIN_CPGS = 2


@dataclass(frozen=True)
class CpGRecord:
    """One CpG's pooled cross-cohort result, positioned on the genome."""

    cpg: str
    chrom: str
    pos: int
    effect: float
    se: float
    p: float
    delta_beta: float = float("nan")
    gene: str = ""

    @property
    def usable(self) -> bool:
        return (
            math.isfinite(self.effect)
            and math.isfinite(self.se)
            and self.se > 0
            and math.isfinite(self.p)
            and self.pos >= 0
        )


@dataclass
class Region:
    """A run of same-direction CpGs and the statistic combining them."""

    chrom: str
    start: int
    end: int
    n_cpgs: int
    cpgs: tuple[str, ...]
    genes: tuple[str, ...]
    direction: str
    estimate: float
    se: float
    z: float
    p_independent: float
    mean_delta_beta: float
    min_cpg_p: float
    p_empirical: float = float("nan")
    p_empirical_is_bound: bool = False
    fdr_empirical: float = float("nan")

    @property
    def width(self) -> int:
        return self.end - self.start + 1

    def as_row(self) -> dict[str, object]:
        return {
            "chr": self.chrom,
            "start": self.start,
            "end": self.end,
            "width": self.width,
            "n_cpgs": self.n_cpgs,
            "genes": ";".join(self.genes),
            "direction": self.direction,
            "estimate": self.estimate,
            "se": self.se,
            "z": self.z,
            "p_independent": self.p_independent,
            "p_empirical": self.p_empirical,
            "p_empirical_is_bound": self.p_empirical_is_bound,
            "fdr_empirical": self.fdr_empirical,
            "mean_delta_beta": self.mean_delta_beta,
            "min_cpg_p": self.min_cpg_p,
            "cpgs": ";".join(self.cpgs),
        }


def _norm_sf(z: float) -> float:
    """Two-sided standard-normal tail probability.

    ``math.erfc`` is used rather than a series expansion because the tail is where region
    statistics live and a truncating approximation would silently flatten it.
    """
    if not math.isfinite(z):
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2.0))


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    base = chrom[3:] if chrom.lower().startswith("chr") else chrom
    try:
        return (int(base), "")
    except ValueError:
        return (10_000, base)


def combine_region(members: Sequence[CpGRecord]) -> tuple[float, float, float]:
    """Combine member CpGs assuming independence.

    Returns ``(estimate, se, z)`` for the summed-effect statistic used by dmrff:
    ``S = sum(b_i)`` with ``SE(S) = sqrt(sum(s_i^2))`` under independence. The independence
    assumption is wrong by construction -- the CpGs are adjacent and co-methylated -- which
    is why the caller is expected to calibrate ``z`` against sign-flipped replicates rather
    than read ``p_independent`` as a P-value.
    """
    total = sum(m.effect for m in members)
    var = sum(m.se * m.se for m in members)
    if var <= 0 or not math.isfinite(var) or not math.isfinite(total):
        return (float("nan"), float("nan"), float("nan"))
    se = math.sqrt(var)
    return (total, se, total / se)


def call_regions(
    records: Iterable[CpGRecord],
    *,
    max_gap: int = DEFAULT_MAX_GAP,
    seed_p: float = DEFAULT_SEED_P,
    min_cpgs: int = DEFAULT_MIN_CPGS,
) -> list[Region]:
    """Group CpGs into candidate regions and compute each region's statistic.

    A candidate region is a maximal run of usable CpGs on one chromosome that share the sign
    of their pooled effect, each pass ``seed_p``, and are separated by no more than
    ``max_gap`` base pairs from the previous member. Runs shorter than ``min_cpgs`` are
    dropped: a "region" of one CpG is the single-CpG analysis under a different name, and
    reporting it here would let single-probe noise re-enter through the region table.
    """
    usable = [r for r in records if r.usable]
    usable.sort(key=lambda r: (_chrom_sort_key(r.chrom), r.pos))

    regions: list[Region] = []
    run: list[CpGRecord] = []

    def flush() -> None:
        if len(run) < min_cpgs:
            run.clear()
            return
        estimate, se, z = combine_region(run)
        if not math.isfinite(z):
            run.clear()
            return
        deltas = [m.delta_beta for m in run if math.isfinite(m.delta_beta)]
        genes = tuple(
            dict.fromkeys(
                part
                for m in run
                for part in str(m.gene).split(";")
                if part and part != "NA"
            )
        )
        regions.append(
            Region(
                chrom=run[0].chrom,
                start=run[0].pos,
                end=run[-1].pos,
                n_cpgs=len(run),
                cpgs=tuple(m.cpg for m in run),
                genes=genes,
                direction="up" if estimate > 0 else "down",
                estimate=estimate,
                se=se,
                z=z,
                p_independent=_norm_sf(z),
                mean_delta_beta=(sum(deltas) / len(deltas)) if deltas else float("nan"),
                min_cpg_p=min(m.p for m in run),
            )
        )
        run.clear()

    for rec in usable:
        if rec.p > seed_p:
            flush()
            continue
        if run:
            same_chrom = rec.chrom == run[-1].chrom
            same_sign = (rec.effect > 0) == (run[-1].effect > 0)
            close = rec.pos - run[-1].pos <= max_gap
            if not (same_chrom and same_sign and close):
                flush()
        run.append(rec)
    flush()

    return regions


def calibrate_regions(
    observed: Sequence[Region],
    null_statistics: Sequence[Sequence[float]],
) -> list[Region]:
    """Attach empirical P-values and an FDR estimate from sign-flipped replicates.

    ``null_statistics`` holds one sequence of absolute region ``z`` values per null
    replicate. The empirical P-value of an observed region is the fraction of null region
    statistics at least as extreme, pooled over replicates, with the usual ``(x + 1) /
    (n + 1)`` correction so that no region is assigned probability zero on finite
    resampling. Regions whose empirical P-value equals that floor are marked as bounded:
    the null cannot resolve past it, and reporting the floor as a point estimate would be a
    precision the design does not have.

    The FDR estimate is the expected number of null regions at or beyond the region's
    statistic -- averaged over replicates -- divided by the number of observed regions at
    or beyond it. It is a direct empirical analogue of the Benjamini-Hochberg ratio and
    needs no independence assumption, because the replicates carry the dependence.
    """
    if not observed:
        return []
    pooled = sorted(
        abs(z) for rep in null_statistics for z in rep if math.isfinite(z)
    )
    n_reps = max(len(null_statistics), 1)
    n_null = len(pooled)
    if n_null == 0:
        return list(observed)

    obs_sorted = sorted(
        (abs(r.z) for r in observed if math.isfinite(r.z)), reverse=True
    )
    floor = 1.0 / (n_null + 1)

    import bisect

    for region in observed:
        if not math.isfinite(region.z):
            continue
        stat = abs(region.z)
        n_ge_null = n_null - bisect.bisect_left(pooled, stat)
        region.p_empirical = (n_ge_null + 1) / (n_null + 1)
        region.p_empirical_is_bound = n_ge_null == 0
        n_ge_obs = len(obs_sorted) - bisect.bisect_left(
            sorted(obs_sorted), stat
        )
        expected_false = n_ge_null / n_reps
        region.fdr_empirical = (
            min(1.0, expected_false / n_ge_obs) if n_ge_obs else float("nan")
        )
    return list(observed)


def summarize(regions: Sequence[Region], fdr_threshold: float = 0.05) -> dict[str, object]:
    """Counts a reader needs before believing any individual region."""
    called = [r for r in regions if math.isfinite(r.fdr_empirical) and r.fdr_empirical < fdr_threshold]
    widths = [r.width for r in called]
    sizes = [r.n_cpgs for r in called]
    deltas = [abs(r.mean_delta_beta) for r in called if math.isfinite(r.mean_delta_beta)]
    return {
        "n_candidate_regions": len(regions),
        "n_regions_fdr": len(called),
        "median_n_cpgs": _median(sizes),
        "median_width_bp": _median(widths),
        "median_abs_mean_delta_beta": _median(deltas),
        "n_up": sum(1 for r in called if r.direction == "up"),
        "n_down": sum(1 for r in called if r.direction == "down"),
    }


def _median(values: Sequence[float]) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2.0


# ---------------------------------------------------------------------------
# Command-line entry point
# ---------------------------------------------------------------------------


def _pool_pattern(records, signs, min_cohorts):
    """Pool every CpG under one sign pattern and return positioned CpGRecords."""
    from illumeta_meta import _random_effect_meta_one  # local: avoids a circular import

    out: list[CpGRecord] = []
    for cpg, rec in records.items():
        ann = rec["annotations"]
        chrom = str(ann.get("chr") or "").strip()
        pos_raw = str(ann.get("pos") or "").strip()
        if not chrom or not pos_raw:
            continue
        try:
            pos = int(float(pos_raw))
        except ValueError:
            continue
        effects = [e * s for e, s in zip(rec["effects"], signs)]
        ses = list(rec["ses"])
        valid = [
            math.isfinite(e) and math.isfinite(v) and v > 0
            for e, v in zip(effects, ses)
        ]
        if sum(valid) < min_cohorts:
            continue
        meta = _random_effect_meta_one(effects, ses, valid)
        deltas = [d * s for d, s in zip(rec["deltas"], signs)]
        tau2 = meta.get("tau2") or 0.0
        num = den = 0.0
        for d, v, ok in zip(deltas, ses, valid):
            if ok and math.isfinite(d):
                w = 1.0 / (v * v + tau2)
                num += w * d
                den += w
        out.append(
            CpGRecord(
                cpg=cpg,
                chrom=chrom,
                pos=pos,
                effect=meta["random_effect"],
                se=meta["random_se"],
                p=meta["random_p"],
                delta_beta=(num / den) if den > 0 else float("nan"),
                gene=str(ann.get("Gene") or ""),
            )
        )
    return out


def run_regions_cli(args) -> int:
    """Region-level meta-analysis with a sign-flipped empirical null.

    The null is not optional. The region statistic assumes independence between adjacent
    CpGs, which is false by construction, so without the empirical calibration there is no
    quantity here a reader could act on. The command therefore always runs the flips and
    refuses rather than emitting an uncalibrated table.
    """
    import csv as _csv
    import itertools
    import json
    from pathlib import Path as _Path

    from illumeta_meta import (  # noqa: E402
        _log,
        _read_branch_records,
        load_meta_manifest,
        load_positional_cohorts,
        _resolve_path,
    )

    project_root = _Path(getattr(args, "project_root", ".")).resolve()
    allow_missing_summary = bool(getattr(args, "allow_missing_summary", False))

    cohorts = []
    if getattr(args, "manifest", None):
        cohorts.extend(
            load_meta_manifest(
                _resolve_path(args.manifest, project_root), project_root, allow_missing_summary
            )
        )
    if getattr(args, "result_dirs", None):
        cohorts.extend(
            load_positional_cohorts(args.result_dirs, project_root, allow_missing_summary)
        )
    if not cohorts:
        raise ValueError("Provide result directories or --manifest.")

    k = len(cohorts)
    n_patterns = 2 ** (k - 1) - 1
    if n_patterns < args.min_null_patterns:
        raise ValueError(
            f"{k} cohorts give only {n_patterns} non-trivial whole-cohort sign patterns, "
            f"below the {args.min_null_patterns} this command requires to estimate a null. "
            "Region P-values would be bounded so coarsely as to be uninformative."
        )
    if n_patterns > args.max_null_patterns:
        raise ValueError(
            f"{k} cohorts give {n_patterns} sign patterns, above --max-null-patterns "
            f"({args.max_null_patterns}). Each pattern re-pools every CpG, so raise the "
            "limit deliberately rather than by accident."
        )

    out_dir = _Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    branches = [b.strip() for b in args.branches.split(",") if b.strip()]
    branch_files = {
        "minfi": "Minfi_DMPs_full.csv",
        "sesame_strict": "Sesame_DMPs_full.csv",
        "sesame_native": "Sesame_Native_DMPs_full.csv",
    }

    summaries = []
    for branch in branches:
        filename = branch_files.get(branch)
        if filename is None:
            raise ValueError(f"Unknown branch: {branch}")
        records, warnings, _tables, _w = _read_branch_records(
            cohorts, branch, filename,
            bool(getattr(args, "allow_missing_branches", False)),
            bool(getattr(args, "tier3_primary", True)),
            bool(getattr(args, "allow_missing_tier3_primary", False)),
        )
        for warning in warnings:
            _log(f"[warn] {warning}")
        _log(f"{branch}: {len(records):,} CpGs, {n_patterns} null patterns")

        patterns = [(1,) + p for p in itertools.product((1, -1), repeat=k - 1)]
        observed = None
        nulls: list[list[float]] = []
        for index, signs in enumerate(patterns):
            regions = call_regions(
                _pool_pattern(records, signs, args.min_cohorts),
                max_gap=args.max_gap,
                seed_p=args.seed_p,
                min_cpgs=args.min_cpgs,
            )
            if index == 0:
                observed = regions
                _log(f"  observed: {len(regions):,} candidate regions")
            else:
                nulls.append([r.z for r in regions])
        calibrate_regions(observed, nulls)

        called = sorted(
            (r for r in observed
             if math.isfinite(r.fdr_empirical) and r.fdr_empirical < args.region_fdr),
            key=lambda r: r.fdr_empirical,
        )
        rows = [r.as_row() for r in (called or observed[:1])]
        dest = out_dir / f"{branch}_regions.tsv"
        with dest.open("w", encoding="utf-8", newline="") as fh:
            writer = _csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            if called:
                writer.writerows(r.as_row() for r in called)
        summary = summarize(observed, args.region_fdr)
        counts = [len(z) for z in nulls]
        summary.update(
            branch=branch,
            n_null_patterns=len(nulls),
            mean_null_candidate_regions=round(sum(counts) / len(counts), 1),
            fdr_threshold=args.region_fdr,
        )
        summaries.append(summary)
        _log(
            f"  {summary['n_regions_fdr']:,} regions at empirical FDR < {args.region_fdr} "
            f"({summary['n_candidate_regions']:,} candidates; null mean "
            f"{summary['mean_null_candidate_regions']:,})"
        )

    dest = out_dir / "region_summary.tsv"
    with dest.open("w", encoding="utf-8", newline="") as fh:
        writer = _csv.DictWriter(fh, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    (out_dir / "region_manifest.json").write_text(
        json.dumps(
            {
                "cohorts": [str(c.result_dir) for c in cohorts],
                "n_cohorts": k,
                "n_null_patterns": n_patterns,
                "max_gap": args.max_gap,
                "seed_p": args.seed_p,
                "min_cpgs": args.min_cpgs,
                "min_cohorts": args.min_cohorts,
                "region_fdr": args.region_fdr,
                "branch_summaries": summaries,
                "null": "whole-cohort sign flipping; region P-values are empirical and "
                        "bounded below by 1/(n_null_regions + 1)",
            },
            indent=2,
        )
    )
    _log(f"Region analysis complete: {out_dir}")
    return 0
