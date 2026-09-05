#!/usr/bin/env python3
"""Convert another tool's per-cohort EWAS tables into IlluMeta meta-analysis inputs.

IlluMeta's cross-cohort layer -- DerSimonian-Laird pooling with tau2 propagation,
leave-one-cohort-out direction stability, directional partial conjunction, I2 and
Knapp-Hartung sensitivity -- does not depend on IlluMeta having produced the
per-cohort estimates. ChAMP, RnBeads, meffil, limma run by hand, or a published
supplementary table all yield the same three quantities the pooling needs: an
effect, its standard error, and the CpG it belongs to.

This module normalises such a table into the layout `illumeta.py meta` already
reads, so the meta layer is usable as a standalone tool without re-running anyone's
preprocessing. The conversion is deliberately a separate step that writes files on
disk rather than an in-memory shim inside the meta reader: the reader is covered by
the metafor equivalence tests and the release checksums, and an external-input path
that reached into it would put unvalidated parsing upstream of validated arithmetic.

The scale check is the part that matters scientifically. Inverse-variance pooling
assumes every cohort's effect is measured on one scale. IlluMeta reports logFC on
the M-value scale; several tools report a difference in beta. Those differ by
roughly a factor of 4 near beta = 0.5 and diverge further at the extremes, so
pooling them together produces a number that is not an effect size in any scale.
`effect_scale` is therefore a required manifest column with no default.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

# The three branch names the meta layer knows, and the file each expects. External
# results are single-pipeline by nature -- ChAMP has no "sesame_native" -- so a
# converted cohort populates one branch and the meta run selects that branch alone.
from illumeta_meta import ANNOTATION_COLUMNS, BRANCH_ALIASES, DEFAULT_BRANCH_FILES

# Manifest columns naming the source table's own column headers. Only cpg, effect
# and pvalue are structurally required; se/t are required as a pair-of-which-one
# (checked below) because SE can be reconstructed from effect/t exactly.
COLUMN_MAP_FIELDS = {
    "col_cpg": "CpG",
    "col_effect": "logFC",
    "col_pvalue": "P.Value",
    "col_se": "SE",
    "col_t": "t",
    "col_delta_beta": "Delta_Beta",
    "col_gene": "Gene",
    "col_chr": "chr",
    "col_pos": "pos",
    "col_region": "Region",
    "col_island_context": "Island_Context",
}

# Header names to try when a mapping column is absent, so the common case needs no
# mapping at all. Compared after _norm_header strips every separator, which is why
# "coefficient.se", "coefficient_se" and "coefficientSE" all reduce to one candidate --
# and why a table carrying two spellings of the same field is an ambiguity to report
# rather than a race between candidate orderings.
#
# The entries are the real column names of the tools this is meant to accept:
# ChAMP returns limma's logFC/t/P.Value; meffil returns coefficient, coefficient.se and
# t.statistic; RnBeads returns mean.diff and diffmeth.p.val but no SE or t at all, so a
# RnBeads site table cannot be imported without one -- that is a property of the tool's
# output, not a gap in this list.
DEFAULT_HEADER_CANDIDATES = {
    "CpG": ("cpg", "probeid", "probe", "id", "name", "cgid", "illuminaid", "site"),
    "logFC": ("logfc", "effect", "estimate", "coefficient", "coef", "beta", "logratio"),
    "P.Value": ("pvalue", "p", "pval", "praw", "diffmethpval", "pvalues"),
    "SE": ("se", "stderror", "stderr", "standarderror", "coefficientse", "sebeta", "sd"),
    "t": ("t", "tvalue", "tstat", "tstatistic", "statistic", "z", "zscore", "zvalue"),
    "Delta_Beta": ("deltabeta", "db", "betadiff", "betafc", "meandiff", "deltab"),
    "Gene": ("gene", "genesymbol", "symbol", "ucscrefgenename"),
    "chr": ("chr", "chrom", "chromosome", "cpgchrm", "seqnames"),
    "pos": ("pos", "position", "start", "mapinfo", "cpgbeg"),
    "Region": ("region", "ucscrefgenegroup", "feature", "generegion"),
    "Island_Context": ("islandcontext", "relationtoisland", "cgi", "relationtoucsccpgisland"),
}

VALID_EFFECT_SCALES = {"m", "beta"}


@dataclass
class ExternalCohortSpec:
    cohort_id: str
    table: Path
    n_con: int
    n_test: int
    effect_scale: str
    branch: str = "minfi"
    tool: str = ""
    tool_version: str = ""
    platform: str = ""
    tissue: str = ""
    label: str = ""
    column_map: dict[str, str] = field(default_factory=dict)


@dataclass
class ConversionResult:
    cohort_id: str
    out_dir: Path
    branch: str
    table_written: Path
    n_rows_in: int
    n_rows_written: int
    n_dropped_no_cpg: int
    n_dropped_duplicate: int
    n_dropped_unusable: int
    se_source: str
    resolved_columns: dict[str, str]
    warnings: list[str]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def _sniff_delimiter(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        return "\t" if "\t" in sample else ","


def _safe_float(value: object) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NULL", "NONE", "."}:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def _safe_int(value: object, field_name: str, cohort_id: str) -> int:
    text = str(value if value is not None else "").strip()
    if not text:
        raise ValueError(
            f"{cohort_id}: '{field_name}' is required. Group sizes cannot be inferred "
            "from a results table, and the meta layer weights the pooled delta beta by them."
        )
    try:
        parsed = int(float(text))
    except ValueError as exc:
        raise ValueError(f"{cohort_id}: '{field_name}' must be an integer, got {text!r}") from exc
    if parsed <= 0:
        raise ValueError(f"{cohort_id}: '{field_name}' must be positive, got {parsed}")
    return parsed


def _norm_header(name: str) -> str:
    """Reduce a header to letters and digits only.

    Separators are dropped rather than kept, so "log_FC", "logFC" and "log.FC" collapse
    to one key. That is what makes a table carrying two of those spellings detectable as
    ambiguous; keeping the underscore would let the first candidate in the list win
    silently, and which of two effect columns got pooled would depend on the order of a
    tuple in this file.
    """
    return "".join(ch for ch in name.strip().lower() if ch.isalnum())


def resolve_columns(
    header: list[str], column_map: dict[str, str], cohort_id: str
) -> tuple[dict[str, str], list[str]]:
    """Map canonical field -> the source table's actual header name.

    An explicit mapping is honoured verbatim and must name a header that exists;
    a typo that silently fell through to auto-detection would be worse than a stop,
    because the run would succeed while pooling the wrong column.
    """
    warnings: list[str] = []
    by_norm: dict[str, list[str]] = {}
    for col in header:
        by_norm.setdefault(_norm_header(col), []).append(col)

    resolved: dict[str, str] = {}
    for manifest_field, canonical in COLUMN_MAP_FIELDS.items():
        explicit = (column_map.get(manifest_field) or "").strip()
        if explicit:
            if explicit not in header:
                raise ValueError(
                    f"{cohort_id}: {manifest_field}={explicit!r} names a column that is not "
                    f"in the table header. Available: {', '.join(header)}"
                )
            resolved[canonical] = explicit
            continue
        for candidate in DEFAULT_HEADER_CANDIDATES.get(canonical, ()):
            matches = by_norm.get(candidate)
            if not matches:
                continue
            if len(matches) > 1:
                raise ValueError(
                    f"{cohort_id}: header contains {len(matches)} columns that normalise to "
                    f"{candidate!r} ({', '.join(matches)}); set {manifest_field} explicitly."
                )
            resolved[canonical] = matches[0]
            warnings.append(f"{cohort_id}: auto-detected {canonical} <- {matches[0]!r}")
            break

    for required in ("CpG", "logFC", "P.Value"):
        if required not in resolved:
            raise ValueError(
                f"{cohort_id}: could not identify the {required} column. Set "
                f"{[k for k, v in COLUMN_MAP_FIELDS.items() if v == required][0]} in the manifest. "
                f"Header was: {', '.join(header)}"
            )
    if "SE" not in resolved and "t" not in resolved:
        raise ValueError(
            f"{cohort_id}: the table needs either a standard-error column or a t/z statistic. "
            "Inverse-variance pooling has no way to weight a cohort without one, and a "
            "P-value alone cannot supply it without assuming the very distribution being tested. "
            f"Header was: {', '.join(header)}"
        )
    return resolved, warnings


def convert_cohort(spec: ExternalCohortSpec, out_root: Path) -> ConversionResult:
    if spec.effect_scale not in VALID_EFFECT_SCALES:
        raise ValueError(
            f"{spec.cohort_id}: effect_scale must be one of {sorted(VALID_EFFECT_SCALES)}, "
            f"got {spec.effect_scale!r}. This column has no default on purpose: pooling an "
            "M-value logFC with a beta-scale difference yields a quantity that is neither."
        )
    branch = BRANCH_ALIASES.get(spec.branch.strip().lower(), spec.branch.strip().lower())
    if branch not in DEFAULT_BRANCH_FILES:
        raise ValueError(
            f"{spec.cohort_id}: branch must be one of {', '.join(DEFAULT_BRANCH_FILES)}, got {spec.branch!r}"
        )
    if not spec.table.exists():
        raise FileNotFoundError(f"{spec.cohort_id}: table not found: {spec.table}")

    delimiter = _sniff_delimiter(spec.table)
    warnings: list[str] = []
    rows_out: list[dict[str, str]] = []
    seen: set[str] = set()
    n_in = n_dup = n_no_cpg = n_unusable = 0
    se_direct = se_from_t = 0

    with _open_text(spec.table) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"{spec.cohort_id}: no header in {spec.table}")
        resolved, resolve_warnings = resolve_columns(list(reader.fieldnames), spec.column_map, spec.cohort_id)
        warnings.extend(resolve_warnings)

        for row in reader:
            n_in += 1
            cpg = str(row.get(resolved["CpG"]) or "").strip()
            if not cpg:
                n_no_cpg += 1
                continue
            if cpg in seen:
                # The meta reader raises on a duplicate CpG. Collapsing here would hide a
                # real problem in the upstream table (two annotations for one probe, a
                # concatenated file), so drop and count instead of guessing which to keep.
                n_dup += 1
                continue
            seen.add(cpg)

            effect = _safe_float(row.get(resolved["logFC"]))
            se = _safe_float(row.get(resolved["SE"])) if "SE" in resolved else math.nan
            if math.isfinite(se) and se > 0:
                se_direct += 1
            else:
                # Exact for any t defined as effect/SE, which covers limma, plain lm and
                # every tool that reports a Wald statistic. Recorded in the provenance so
                # a reader knows the SE was derived rather than reported.
                t_stat = _safe_float(row.get(resolved["t"])) if "t" in resolved else math.nan
                if math.isfinite(effect) and math.isfinite(t_stat) and t_stat != 0:
                    se = abs(effect / t_stat)
                    se_from_t += 1
                else:
                    se = math.nan
            if not (math.isfinite(effect) and math.isfinite(se) and se > 0):
                n_unusable += 1
                continue

            out: dict[str, str] = {
                "CpG": cpg,
                "logFC": repr(effect),
                "SE": repr(se),
                "P.Value": str(row.get(resolved["P.Value"]) or "").strip(),
                "t": repr(effect / se),
            }
            delta = _safe_float(row.get(resolved["Delta_Beta"])) if "Delta_Beta" in resolved else math.nan
            out["Delta_Beta"] = repr(delta) if math.isfinite(delta) else ""
            for col in ANNOTATION_COLUMNS:
                src = resolved.get(col)
                out[col] = str(row.get(src) or "").strip() if src else ""
            rows_out.append(out)

    if not rows_out:
        raise ValueError(
            f"{spec.cohort_id}: no usable rows in {spec.table} "
            f"({n_in} read, {n_no_cpg} without a CpG, {n_dup} duplicates, {n_unusable} without a finite effect and SE)"
        )
    if "Delta_Beta" not in resolved:
        warnings.append(
            f"{spec.cohort_id}: no Delta_Beta column, so the pooled delta-beta magnitude gate "
            "(--min-abs-delta-beta) cannot be evaluated for this cohort. Pass "
            "--min-abs-delta-beta 0 or supply col_delta_beta."
        )
    if spec.effect_scale == "beta":
        warnings.append(
            f"{spec.cohort_id}: effects are on the beta scale. Every cohort in the meta run must "
            "use the same scale; do not mix this cohort with M-value logFC inputs."
        )

    out_dir = out_root / spec.cohort_id / "external_results"
    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / DEFAULT_BRANCH_FILES[branch]
    fieldnames = ["CpG", "logFC", "SE", "t", "P.Value", "Delta_Beta", *ANNOTATION_COLUMNS]
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    se_source = "reported" if se_from_t == 0 else ("derived_from_t" if se_direct == 0 else "mixed")
    summary = {
        # Consumed by the meta reader.
        "n_con": spec.n_con,
        "n_test": spec.n_test,
        # Not a tier3 mode, so the tier3 fixed-effect variance warning correctly stays
        # silent; the value still tells a reader these estimates are not IlluMeta's.
        "primary_result_mode": "external_import",
        "primary_lambda_guard_status": "not_evaluated",
        "primary_tier3_meta_method": "",
        # Provenance. IlluMeta did not preprocess these samples and must not imply it did.
        "external_import": {
            "tool": spec.tool,
            "tool_version": spec.tool_version,
            "source_table": str(spec.table),
            "effect_scale": spec.effect_scale,
            "branch": branch,
            "se_source": se_source,
            "se_reported": se_direct,
            "se_derived_from_t": se_from_t,
            "rows_read": n_in,
            "rows_written": len(rows_out),
            "rows_dropped_no_cpg": n_no_cpg,
            "rows_dropped_duplicate": n_dup,
            "rows_dropped_unusable": n_unusable,
            "resolved_columns": resolved,
            "converted_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "warnings": warnings,
        },
        "platform": spec.platform,
        "tissue": spec.tissue,
        "label": spec.label,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    return ConversionResult(
        cohort_id=spec.cohort_id,
        out_dir=out_dir,
        branch=branch,
        table_written=table_path,
        n_rows_in=n_in,
        n_rows_written=len(rows_out),
        n_dropped_no_cpg=n_no_cpg,
        n_dropped_duplicate=n_dup,
        n_dropped_unusable=n_unusable,
        se_source=se_source,
        resolved_columns=resolved,
        warnings=warnings,
    )


def load_external_manifest(manifest: Path, project_root: Path) -> list[ExternalCohortSpec]:
    delimiter = _sniff_delimiter(manifest)
    specs: list[ExternalCohortSpec] = []
    with _open_text(manifest) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"External manifest has no header: {manifest}")
        for i, row in enumerate(reader, start=1):
            include = str(row.get("include") or "").strip().lower()
            if include in {"0", "false", "no", "n", "exclude", "skip"}:
                continue
            cohort_id = str(row.get("cohort") or row.get("cohort_id") or "").strip()
            if not cohort_id:
                raise ValueError(f"{manifest}: row {i} has no 'cohort' value")
            table_text = str(row.get("table") or row.get("path") or "").strip()
            if not table_text:
                raise ValueError(f"{manifest}: row {i} ({cohort_id}) has no 'table' value")
            table = Path(table_text).expanduser()
            if not table.is_absolute():
                table = project_root / table
            specs.append(
                ExternalCohortSpec(
                    cohort_id=cohort_id,
                    table=table.resolve(),
                    n_con=_safe_int(row.get("n_con"), "n_con", cohort_id),
                    n_test=_safe_int(row.get("n_test"), "n_test", cohort_id),
                    effect_scale=str(row.get("effect_scale") or "").strip().lower(),
                    branch=str(row.get("branch") or "minfi").strip(),
                    tool=str(row.get("tool") or "").strip(),
                    tool_version=str(row.get("tool_version") or "").strip(),
                    platform=str(row.get("platform") or "").strip(),
                    tissue=str(row.get("tissue") or "").strip(),
                    label=str(row.get("label") or "").strip(),
                    column_map={k: str(row.get(k) or "").strip() for k in COLUMN_MAP_FIELDS},
                )
            )
    if not specs:
        raise ValueError(f"{manifest}: no included rows")
    return specs


def run_import_external(
    manifest: Path, out_root: Path, project_root: Path, log=print
) -> tuple[list[ConversionResult], Path]:
    specs = load_external_manifest(manifest, project_root)
    branches = {BRANCH_ALIASES.get(s.branch.lower(), s.branch.lower()) for s in specs}
    scales = {s.effect_scale for s in specs}
    if len(scales) > 1:
        raise ValueError(
            "All cohorts in one meta-analysis must share an effect scale; this manifest mixes "
            f"{sorted(scales)}. Convert to a single scale before importing."
        )

    out_root.mkdir(parents=True, exist_ok=True)
    results = [convert_cohort(spec, out_root) for spec in specs]
    for res in results:
        log(f"  {res.cohort_id}: {res.n_rows_written:,}/{res.n_rows_in:,} rows -> {res.table_written}")
        for warning in res.warnings:
            log(f"    ! {warning}")

    # The meta manifest is written for the user rather than described in prose, because
    # the next command is `illumeta.py meta -m <this file>` and a hand-built manifest is
    # where a path typo silently drops a cohort from the pooling.
    meta_manifest = out_root / "external_meta_manifest.tsv"
    with meta_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cohort", "result_dir", "label", "platform", "tissue"])
        for spec, res in zip(specs, results):
            writer.writerow([res.cohort_id, str(res.out_dir), spec.label, spec.platform, spec.tissue])

    log(f"\nWrote {meta_manifest}")
    log("Next:")
    log(f"  illumeta.py meta -m {meta_manifest} --branches {','.join(sorted(branches))} -o <output_dir>")
    if len(branches) == 1:
        log(
            "  (single-pipeline inputs: the dual-route consensus columns are not meaningful here, "
            "so select only the branch that was imported)"
        )
    return results, meta_manifest
