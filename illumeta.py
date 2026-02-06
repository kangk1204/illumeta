#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import csv
import html
import json
import re
import statistics
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime
from urllib import request as _urllib_request
from urllib.error import URLError, HTTPError

__version__ = "1.0.0"

# Configuration for R script paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def safe_path_component(name: str, fallback: str = "results") -> str:
    """Return an ASCII-only path component safe across filesystems/locales."""
    if not name:
        return fallback
    normalized = unicodedata.normalize("NFKD", name)
    ascii_name = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_name = re.sub(r"\s+", "_", ascii_name)
    ascii_name = re.sub(r"[^A-Za-z0-9._-]+", "_", ascii_name)
    ascii_name = re.sub(r"_+", "_", ascii_name).strip("._-")
    return ascii_name or fallback

def sniff_delimiter(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            first = handle.readline()
    except OSError:
        return "\t"
    return "\t" if "\t" in first else ","

def load_config_rows(path: str):
    delim = sniff_delimiter(path)
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        rows = [row for row in reader]
    return rows, reader.fieldnames or [], delim

def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False

def profile_columns(rows, headers):
    profiles = []
    for col in headers:
        raw_vals = [row.get(col, "") for row in rows]
        vals = []
        for v in raw_vals:
            if v is None:
                vals.append("")
            elif isinstance(v, str):
                vals.append(v.strip())
            else:
                vals.append(str(v).strip())
        missing = [v == "" for v in vals]
        non_missing = [v for v, m in zip(vals, missing) if not m]
        unique_vals = set(non_missing)
        numeric_count = sum(1 for v in non_missing if is_number(v))
        numeric_frac = numeric_count / len(non_missing) if non_missing else 0.0
        dtype = "numeric" if numeric_frac >= 0.9 else "categorical"
        profiles.append({
            "column": col,
            "missing_count": sum(missing),
            "missing_frac": (sum(missing) / len(vals)) if vals else 0.0,
            "unique_non_missing": len(unique_vals),
            "dtype_guess": dtype,
            "example_values": non_missing[:3],
        })
    return profiles

def resolve_id_column(headers, override=None):
    if override:
        if override not in headers:
            raise ValueError(f"Specified --id-column '{override}' not found in configure.tsv.")
        return override
    for col in headers:
        if re.search(r"(GSM|geo_accession)", col, flags=re.IGNORECASE):
            return col
    if "Basename" in headers:
        return "Basename"
    raise ValueError("Could not identify a sample ID column. Provide --id-column or include GSM/geo_accession/Basename in configure.tsv.")

GROUP_COLUMN_HINTS = (
    "group",
    "condition",
    "disease",
    "phenotype",
    "status",
    "case",
    "control",
    "treatment",
    "response",
    "diagnosis",
)
EXCLUDE_GROUP_HINTS = (
    "id",
    "gsm",
    "geo_accession",
    "sample",
    "title",
    "source",
    "description",
    "protocol",
    "characteristics",
    "relation",
    "contact",
    "submission",
    "release",
    "update",
    "date",
    "platform",
    "series",
    "supplementary",
    "file",
    "url",
)

TIER1_KEYWORDS = {
    "diagnosis", "disease", "cancer", "tumor", "subtype", "response", "treatment",
    "genotype", "mutation", "histology", "grade", "stage", "phenotype", "case", "control",
}
TIER2_KEYWORDS = {
    "sex", "gender", "age", "agegroup", "age_group", "race", "ethnicity", "smoker",
    "smoking", "bmi", "gestationalage", "gestational_age",
}
TIER3_KEYWORDS = {
    "batch", "plate", "slide", "array", "sentrix", "position", "chip", "run",
    "center", "hospital", "site", "date",
}

ROMAN_MAP = {
    "i": 1, "ii": 2, "iii": 3, "iv": 4, "v": 5, "vi": 6, "vii": 7, "viii": 8, "ix": 9, "x": 10,
}
ORDERED_SETS = [
    ["low", "medium", "high"],
    ["mild", "moderate", "severe"],
    ["early", "mid", "late"],
]
CONTROL_SYNONYMS = {
    "control",
    "controls",
    "normal",
    "healthy",
    "reference",
    "wildtype",
    "wt",
    "baseline",
    "vehicle",
    "mock",
    "untreated",
}
TEST_SYNONYMS = {
    "case",
    "cases",
    "tumor",
    "disease",
    "affected",
    "patient",
    "treated",
    "exposed",
    "stimulated",
    "ko",
    "knockout",
}

def normalize_group_value(value: str) -> str:
    if value is None:
        return ""
    if not isinstance(value, str):
        value = str(value)
    value = value.strip().lower()
    value = re.sub(r"[\s_/-]+", "", value)
    value = re.sub(r"[^a-z0-9]+", "", value)
    return value

def parse_group_map(map_str: str, group_con: str, group_test: str) -> dict:
    mapping = {}
    if not map_str:
        return mapping
    for chunk in re.split(r"[;,]", map_str):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError("Invalid --group-map entry (expected key=value): " + chunk)
        raw, mapped = chunk.split("=", 1)
        mapped = mapped.strip()
        mapped_norm = normalize_group_value(mapped)
        if mapped_norm in {"control", "con", "ctrl"}:
            mapped = group_con
        elif mapped_norm in {"test", "case"}:
            mapped = group_test
        for raw_item in re.split(r"[|/]", raw):
            raw_item = raw_item.strip()
            if not raw_item:
                continue
            mapping[normalize_group_value(raw_item)] = mapped
    return mapping

def profile_values(values):
    cleaned = []
    for v in values:
        if v is None:
            cleaned.append("")
        elif isinstance(v, str):
            cleaned.append(v.strip())
        else:
            cleaned.append(str(v).strip())
    missing = [v == "" for v in cleaned]
    non_missing = [v for v in cleaned if v != ""]
    unique_vals = sorted(set(non_missing))
    non_missing_count = len(non_missing)
    numeric_count = sum(1 for v in non_missing if is_number(v))
    numeric_frac = numeric_count / non_missing_count if non_missing_count else 0.0
    return {
        "missing_frac": (sum(missing) / len(cleaned)) if cleaned else 0.0,
        "unique_values": unique_vals,
        "unique_count": len(unique_vals),
        "numeric_frac": numeric_frac,
        "non_missing_count": non_missing_count,
        "non_missing_frac": (non_missing_count / len(cleaned)) if cleaned else 0.0,
    }

def guess_tier(name_lower: str):
    if any(k in name_lower for k in TIER1_KEYWORDS):
        return "primary", 100
    if any(k in name_lower for k in TIER3_KEYWORDS):
        return "technical", -100
    if any(k in name_lower for k in TIER2_KEYWORDS):
        return "demographic", -5
    return "unknown", 0

def ordinal_rank(value: str):
    norm = normalize_group_value(value)
    if not norm:
        return None
    if norm.isdigit():
        return int(norm)
    for prefix in ("stage", "grade"):
        if norm.startswith(prefix):
            tail = norm[len(prefix):]
            if tail.isdigit():
                return int(tail)
            if tail in ROMAN_MAP:
                return ROMAN_MAP[tail]
    if norm in ROMAN_MAP:
        return ROMAN_MAP[norm]
    for seq in ORDERED_SETS:
        for idx, token in enumerate(seq, start=1):
            if norm == token:
                return idx
    return None

def is_ordinal(values):
    non_missing = [v for v in values if v not in ("", None)]
    if len(non_missing) < 3:
        return False, {}
    ranks = {}
    valid = 0
    for v in non_missing:
        rank = ordinal_rank(v)
        if rank is not None:
            valid += 1
            ranks.setdefault(v, rank)
    if valid / max(1, len(non_missing)) < 0.7:
        return False, {}
    if len(set(ranks.values())) < 3:
        return False, {}
    return True, ranks

def score_group_candidate(name: str, values):
    prof = profile_values(values)
    unique_count = prof["unique_count"]
    name_lower = name.lower()
    if any(h in name_lower for h in EXCLUDE_GROUP_HINTS):
        if "source_name" not in name_lower:
            return None
    tier, tier_score = guess_tier(name_lower)
    keyword_score = 2 if any(k in name_lower for k in GROUP_COLUMN_HINTS) else 0
    value_norms = {normalize_group_value(v) for v in prof["unique_values"]}
    synonym_score = 0
    if value_norms & CONTROL_SYNONYMS:
        synonym_score += 2
    if value_norms & TEST_SYNONYMS:
        synonym_score += 2
    substring_score = 0
    for val in value_norms:
        if any(s in val for s in CONTROL_SYNONYMS):
            substring_score += 1
        if any(s in val for s in TEST_SYNONYMS):
            substring_score += 1
    substring_score = 1 if substring_score > 0 else 0
    if unique_count < 2:
        return None
    ordinal_flag, ordinal_map = is_ordinal(prof["unique_values"])
    if prof["numeric_frac"] >= 0.8 and unique_count >= 7 and not ordinal_flag:
        return None
    signal_score = keyword_score + synonym_score + substring_score
    max_unique_allow = max(6, int(len(values) * 0.25))
    if unique_count > max_unique_allow and signal_score < 2:
        return None
    counts = [values.count(v) for v in prof["unique_values"]]
    min_ct = min(counts) if counts else 0
    max_ct = max(counts) if counts else 0
    balance_ratio = (min_ct / max_ct) if max_ct else 0
    balance_score = 2 if balance_ratio >= 0.5 else 1 if balance_ratio >= 0.25 else 0
    unique_score = 2 if unique_count == 2 else 1 if unique_count == 3 else 0
    unique_penalty = -1 if unique_count >= 5 and not ordinal_flag else 0
    if unique_count >= 7 and not ordinal_flag:
        unique_penalty -= 1
    missing_score = 1 if prof["missing_frac"] < 0.1 else 0
    coverage = prof["non_missing_frac"]
    if coverage >= 0.9:
        coverage_score = 2
    elif coverage >= 0.75:
        coverage_score = 1
    elif coverage >= 0.6:
        coverage_score = 0
    else:
        coverage_score = -1
    min_count_penalty = -1 if (min_ct < 2 and len(values) >= 6) else 0
    max_frac = (max_ct / sum(counts)) if counts else 0
    imbalance_penalty = -1 if max_frac >= 0.85 else 0
    total = (
        tier_score
        + keyword_score
        + synonym_score
        + substring_score
        + unique_score
        + unique_penalty
        + missing_score
        + balance_score
        + coverage_score
        + min_count_penalty
        + imbalance_penalty
    )
    return {
        "name": name,
        "values": values,
        "score": total,
        "tier": tier,
        "tier_score": tier_score,
        "missing_frac": prof["missing_frac"],
        "unique_count": unique_count,
        "ordinal": ordinal_flag,
        "keyword_score": keyword_score,
        "synonym_score": synonym_score,
        "substring_score": substring_score,
        "balance_ratio": balance_ratio,
        "balance_score": balance_score,
        "coverage": coverage,
        "non_missing_count": prof["non_missing_count"],
        "numeric_frac": prof["numeric_frac"],
        "ordinal_map": ordinal_map,
    }

def extract_characteristics_keys(rows, headers):
    key_map = {}
    for idx, row in enumerate(rows):
        for col in headers:
            if "characteristics" not in col.lower():
                continue
            val = (row.get(col) or "").strip()
            if ":" not in val:
                continue
            key, value = val.split(":", 1)
            key = key.strip()
            value = value.strip()
            if not key:
                continue
            key_norm = normalize_group_value(key)
            if not key_norm:
                continue
            entry = key_map.setdefault(key_norm, {"label": key, "values": [""] * len(rows)})
            entry["values"][idx] = value
    return key_map

def collect_group_candidates(rows, headers, id_column=None):
    candidates = []
    for col in headers:
        if col in {"primary_group", "Basename"}:
            continue
        if id_column and col == id_column:
            continue
        cand = score_group_candidate(col, [row.get(col, "") for row in rows])
        if cand:
            cand["source"] = f"column:{col}"
            candidates.append(cand)
    for entry in extract_characteristics_keys(rows, headers).values():
        cand = score_group_candidate(entry["label"], entry["values"])
        if cand:
            cand["source"] = f"characteristics:{entry['label']}"
            candidates.append(cand)
    return candidates

def write_group_candidates(path, candidates):
    if not candidates:
        return
    fieldnames = [
        "rank", "name", "source", "score", "tier", "tier_score",
        "unique_count", "ordinal", "balance_ratio", "coverage", "missing_frac", "numeric_frac",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, cand in enumerate(candidates, start=1):
            writer.writerow(
                {
                    "rank": idx,
                    "name": cand.get("name", ""),
                    "source": cand.get("source", ""),
                    "score": cand.get("score", ""),
                    "tier": cand.get("tier", ""),
                    "tier_score": cand.get("tier_score", ""),
                    "unique_count": cand.get("unique_count", ""),
                    "ordinal": cand.get("ordinal", ""),
                    "balance_ratio": f"{cand.get('balance_ratio', 0):.3f}",
                    "coverage": f"{cand.get('coverage', 0):.3f}",
                    "missing_frac": f"{cand.get('missing_frac', 0):.3f}",
                    "numeric_frac": f"{cand.get('numeric_frac', 0):.3f}",
                }
            )

def infer_group_source(rows, headers, id_column=None, allow_technical=False):
    candidates = collect_group_candidates(rows, headers, id_column=id_column)
    if not candidates:
        return None
    tier_buckets = {"primary": [], "unknown": [], "demographic": [], "technical": []}
    for cand in candidates:
        tier_buckets.setdefault(cand.get("tier", "unknown"), []).append(cand)

    preferred = []
    if tier_buckets["primary"]:
        preferred = tier_buckets["primary"]
    elif tier_buckets["unknown"]:
        preferred = tier_buckets["unknown"]
    elif tier_buckets["demographic"]:
        preferred = tier_buckets["demographic"]
    elif allow_technical:
        preferred = tier_buckets["technical"]

    if not preferred:
        return None

    preferred.sort(
        key=lambda c: (
            c["score"],
            c.get("coverage", 0),
            c.get("balance_ratio", 0),
            -c.get("missing_frac", 0),
            -c.get("unique_count", 0),
        ),
        reverse=True,
    )
    return preferred[0]

def apply_group_mapping(value, mapping, group_con, group_test):
    raw = value or ""
    raw = raw.strip() if isinstance(raw, str) else str(raw).strip()
    if not raw:
        return ""
    norm = normalize_group_value(raw)
    if mapping and norm in mapping:
        return mapping[norm]
    if norm in CONTROL_SYNONYMS:
        return group_con
    if norm in TEST_SYNONYMS:
        return group_test
    return raw

def write_config_rows(path, headers, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({h: row.get(h, "") for h in headers})

def auto_group_config(
    config_path: str,
    group_con: str,
    group_test: str,
    *,
    group_column: str = None,
    group_key: str = None,
    group_map: str = None,
    output_path: str = None,
    id_column: str = None,
    overwrite: bool = False,
    allow_technical: bool = False,
):
    rows, headers, delim = load_config_rows(config_path)
    if delim != "\t":
        raise ValueError("configure.tsv must be tab-delimited (TSV). CSV/other delimiters are not supported.")
    if not headers:
        raise ValueError("configure.tsv appears empty or unreadable.")
    if "primary_group" not in headers:
        headers = ["primary_group"] + headers
        for row in rows:
            row["primary_group"] = ""
    else:
        headers = ["primary_group"] + [h for h in headers if h != "primary_group"]
    pending_idx = []
    for idx, row in enumerate(rows):
        cur = (row.get("primary_group") or "").strip()
        if overwrite or not cur:
            pending_idx.append(idx)
    if not pending_idx:
        return config_path, {"updated": False, "source": "primary_group"}

    mapping = parse_group_map(group_map, group_con, group_test)
    source_values = None
    source_label = None
    if group_column:
        if group_column not in headers:
            raise ValueError(f"Auto-group column '{group_column}' not found in configure.tsv.")
        source_values = [row.get(group_column, "") for row in rows]
        source_label = f"column:{group_column}"
    elif group_key:
        key_map = extract_characteristics_keys(rows, headers)
        key_norm = normalize_group_value(group_key)
        entry = key_map.get(key_norm)
        if not entry:
            raise ValueError(f"Auto-group key '{group_key}' not found in characteristics fields.")
        source_values = entry["values"]
        source_label = f"characteristics:{entry['label']}"
    else:
        try:
            resolved_id = resolve_id_column(headers, override=id_column)
        except ValueError:
            resolved_id = None
        inferred = infer_group_source(rows, headers, id_column=resolved_id, allow_technical=allow_technical)
        if not inferred:
            raise ValueError(
                "Unable to auto-detect a group column. Provide --group-column or --group-key "
                "to populate primary_group automatically."
            )
        source_values = inferred["values"]
        source_label = inferred["source"]

    # Emit candidate report for transparency/debugging.
    candidates_path = None
    try:
        resolved_id = resolve_id_column(headers, override=id_column)
    except ValueError:
        resolved_id = None
    candidates = collect_group_candidates(rows, headers, id_column=resolved_id)
    if candidates:
        candidates.sort(
            key=lambda c: (
                c.get("score", 0),
                c.get("coverage", 0),
                c.get("balance_ratio", 0),
                -c.get("missing_frac", 0),
                -c.get("unique_count", 0),
            ),
            reverse=True,
        )
        report_path = os.path.join(os.path.dirname(os.path.abspath(config_path)), "auto_group_candidates.tsv")
        write_group_candidates(report_path, candidates)
        candidates_path = report_path

    for idx in pending_idx:
        mapped = apply_group_mapping(source_values[idx], mapping, group_con, group_test)
        rows[idx]["primary_group"] = mapped

    still_missing = [i for i in pending_idx if not (rows[i].get("primary_group") or "").strip()]
    if still_missing:
        preview = ", ".join(str(i + 1) for i in still_missing[:5])
        raise ValueError(f"Auto-group left {len(still_missing)} rows empty (rows: {preview}).")

    out_path = output_path or os.path.join(os.path.dirname(os.path.abspath(config_path)), "configure_autogroup.tsv")
    write_config_rows(out_path, headers, rows)
    counts = Counter((row.get("primary_group") or "").strip() for row in rows if (row.get("primary_group") or "").strip())
    return out_path, {
        "updated": True,
        "source": source_label,
        "output_path": out_path,
        "candidates_path": candidates_path,
        "group_counts": dict(counts),
        "filled_rows": len(pending_idx),
        "total_rows": len(rows),
        "group_map": mapping,
    }

def list_idat_basenames(idat_dir: str):
    if not idat_dir or not os.path.isdir(idat_dir):
        return []
    basenames = []
    for name in os.listdir(idat_dir):
        if name.endswith("_Grn.idat") or name.endswith("_Grn.idat.gz"):
            base = re.sub(r"_Grn\.idat(\.gz)?$", "", name)
            basenames.append(os.path.join(idat_dir, base))
    return basenames

def find_basename_for_id(sample_id: str, basenames):
    if not sample_id:
        return None
    for bn in basenames:
        if os.path.basename(bn) == sample_id:
            return bn
    # Fallback: match basenames that start with sample_id followed by '_' or end of string.
    # This avoids GSM1234 matching GSM12345 (substring collision).
    for bn in basenames:
        bname = os.path.basename(bn)
        if bname.startswith(sample_id) and (len(bname) == len(sample_id) or bname[len(sample_id)] == '_'):
            return bn
    return None

def resolve_basename(candidate: str, sample_id: str, basenames, project_dir: str, idat_dir: str):
    cand = candidate or ""
    if not cand:
        return find_basename_for_id(sample_id, basenames)
    if not os.path.isabs(cand):
        cand_project = os.path.join(project_dir, cand)
        cand_idat = os.path.join(idat_dir, cand)
        cand = cand_project
        if not has_idat_pair(cand) and idat_dir != project_dir and has_idat_pair(cand_idat):
            cand = cand_idat
    if not has_idat_green(cand):
        fallback = find_basename_for_id(sample_id, basenames)
        if fallback:
            return fallback
    return cand

def has_idat_green(basename: str) -> bool:
    if not basename:
        return False
    return os.path.exists(basename + "_Grn.idat") or os.path.exists(basename + "_Grn.idat.gz")

def has_idat_red(basename: str) -> bool:
    if not basename:
        return False
    return os.path.exists(basename + "_Red.idat") or os.path.exists(basename + "_Red.idat.gz")

def has_idat_pair(basename: str) -> bool:
    return has_idat_green(basename) and has_idat_red(basename)

def preflight_analysis(config_path: str, idat_dir: str, group_con: str, group_test: str,
                       min_total_size: int, id_column: str = None, drop_missing_idat: bool = True):
    rows, headers, delim = load_config_rows(config_path)
    if delim != "\t":
        raise ValueError("configure.tsv must be tab-delimited (TSV). CSV/other delimiters are not supported.")
    if not headers:
        raise ValueError("configure.tsv appears empty or unreadable.")
    if "primary_group" not in headers:
        raise ValueError("configure.tsv missing 'primary_group' column.")
    id_col = resolve_id_column(headers, override=id_column)
    project_dir = os.path.dirname(os.path.abspath(config_path))
    idat_dir = idat_dir or os.path.join(project_dir, "idat")
    if not os.path.isdir(idat_dir):
        raise ValueError(f"IDAT directory not found: {idat_dir}")

    warnings = []
    n_rows_raw = len(rows)
    filtered_config_path = config_path
    missing_idat_preview = []

    basenames = list_idat_basenames(idat_dir)
    sample_ids = []
    missing_id = []
    duplicate_ids = []
    missing_pairs = []

    for row in rows:
        sample_id = (row.get(id_col) or "").strip()
        if not sample_id:
            missing_id.append(row)
            continue
        sample_ids.append(sample_id)
        if sample_ids.count(sample_id) > 1 and sample_id not in duplicate_ids:
            duplicate_ids.append(sample_id)

        cand = (row.get("Basename") or "").strip()
        basename = resolve_basename(cand, sample_id, basenames, project_dir, idat_dir)
        if not has_idat_pair(basename):
            missing_pairs.append((sample_id, basename))

    if missing_id:
        raise ValueError("Sample ID column contains missing/empty values. Please fill it before running analysis.")
    if duplicate_ids:
        raise ValueError(f"Duplicate sample IDs detected: {', '.join(duplicate_ids)}")
    if missing_pairs:
        missing_preview = [sid for sid, _ in missing_pairs[:10]]
        if not drop_missing_idat:
            preview_str = ", ".join(missing_preview)
            raise ValueError(f"Missing IDAT pairs for {len(missing_pairs)} samples (e.g. {preview_str}).")
        missing_ids = {sid for sid, _ in missing_pairs}
        filtered_rows = [row for row in rows if (row.get(id_col) or "").strip() not in missing_ids]
        if not filtered_rows:
            preview_str = ", ".join(missing_preview)
            raise ValueError(f"All samples are missing IDAT pairs (e.g. {preview_str}).")
        root, ext = os.path.splitext(config_path)
        filtered_config_path = f"{root}_idat{ext or '.tsv'}"
        write_config_rows(filtered_config_path, headers, filtered_rows)
        rows = filtered_rows
        missing_idat_preview = missing_preview
        warnings.append(
            f"Filtered out {len(missing_ids)} samples with missing IDAT pairs; "
            f"using {filtered_config_path}"
        )

    con_lower = (group_con or "").strip().lower()
    test_lower = (group_test or "").strip().lower()
    if con_lower == test_lower:
        raise ValueError("Control and test group labels are identical; please provide distinct labels.")
    group_vals = [(row.get("primary_group") or "").strip().lower() for row in rows]
    if any(g == "" for g in group_vals):
        raise ValueError("primary_group column contains missing/empty values. Please fill before analysis.")
    n_con = sum(1 for g in group_vals if g == con_lower)
    n_test = sum(1 for g in group_vals if g == test_lower)
    n_total = n_con + n_test
    if n_total == 0:
        raise ValueError(f"No samples found matching groups: {group_con} or {group_test}.")
    if min_total_size and n_total < min_total_size:
        raise ValueError(f"Too few samples after group selection: {n_total} < min_total_size ({min_total_size}).")

    other_groups = sorted(set(g for g in group_vals if g and g not in {con_lower, test_lower}))
    if other_groups:
        preview = ", ".join(other_groups[:5])
        warnings.append(
            f"Found {len(other_groups)} additional group label(s) that will be ignored "
            f"({preview}{'...' if len(other_groups) > 5 else ''})."
        )
    batch_patterns = r"(sentrix|slide|array|plate|chip|batch)"
    batch_candidates = [h for h in headers if re.search(batch_patterns, h, flags=re.IGNORECASE)]
    if not batch_candidates:
        warnings.append(
            "No batch candidate columns detected (e.g., Sentrix_ID/Position/plate). "
            "Batch correction may be limited."
        )
    profiles = profile_columns(rows, headers)
    high_missing = [p["column"] for p in profiles
                    if p["missing_frac"] >= 0.5 and p["column"] not in (id_col, "Basename", "primary_group")]
    if high_missing:
        warnings.append(
            "High missingness columns (>=50% empty): " + ", ".join(high_missing[:8]) +
            ("..." if len(high_missing) > 8 else "")
        )
    constant_cols = [p["column"] for p in profiles
                     if p["unique_non_missing"] <= 1 and p["column"] not in (id_col, "Basename", "primary_group")]
    if constant_cols:
        warnings.append(
            "Constant columns will be ignored: " + ", ".join(constant_cols[:8]) +
            ("..." if len(constant_cols) > 8 else "")
        )

    return {
        "config_path": filtered_config_path,
        "idat_dir": idat_dir,
        "sample_count_raw": n_rows_raw,
        "sample_count": len(rows),
        "group_con": n_con,
        "group_test": n_test,
        "missing_idat_pairs": len(missing_pairs),
        "missing_idat_preview": missing_idat_preview,
        "warnings": warnings,
        "column_profile": profiles,
        "batch_candidates": batch_candidates,
        "group_labels": {
            "control": group_con,
            "test": group_test,
            "other": other_groups,
        },
    }

def detect_r_major_minor(env=None):
    override = (env or {}).get("ILLUMETA_R_LIB_VERSION") or os.environ.get("ILLUMETA_R_LIB_VERSION")
    if override:
        return override
    cmd = [
        "Rscript",
        "-e",
        "cat(paste0(R.version$major, '.', sub('\\\\..*', '', R.version$minor)))",
    ]
    try:
        res = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except FileNotFoundError:
        res = None
    if res and res.returncode == 0:
        version = res.stdout.strip()
        if version:
            return version
    try:
        res = subprocess.run(
            ["R", "--version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except FileNotFoundError:
        return None
    if res.returncode != 0:
        return None
    match = re.search(r"R version (\d+)\.(\d+)", res.stdout)
    if not match:
        return None
    return f"{match.group(1)}.{match.group(2)}"

def default_r_lib_base(base_dir: str, env=None) -> str:
    allow_external = (env or {}).get("ILLUMETA_ALLOW_EXTERNAL_LIB") == "1" or (
        os.environ.get("ILLUMETA_ALLOW_EXTERNAL_LIB") == "1"
    )
    if sys.platform == "darwin" and base_dir.startswith("/Volumes/") and not allow_external:
        return os.path.join(os.path.expanduser("~"), ".illumeta", "r-lib")
    return os.path.join(base_dir, ".r-lib")

def default_r_lib(base_dir: str, env=None) -> str:
    lib_root = default_r_lib_base(base_dir, env)
    r_version = detect_r_major_minor(env)
    if r_version:
        return os.path.join(lib_root, f"R-{r_version}")
    return lib_root

def read_recorded_r_lib(base_dir: str):
    record_path = os.path.join(base_dir, ".illumeta_r_lib_path")
    if not os.path.isfile(record_path):
        return None
    try:
        with open(record_path, "r", encoding="utf-8") as handle:
            value = handle.read().strip()
    except OSError:
        return None
    return value or None
R_SCRIPTS_DIR = os.path.join(BASE_DIR, "r_scripts")
DOWNLOAD_SCRIPT = os.path.join(R_SCRIPTS_DIR, "download.R")
ANALYZE_SCRIPT = os.path.join(R_SCRIPTS_DIR, "analyze.R")
SETUP_MARKER = os.path.join(BASE_DIR, ".r_setup_done")
SETUP_SCRIPT = os.path.join(R_SCRIPTS_DIR, "setup_env.R")
DEFAULT_CONDA_PREFIX = os.environ.get("CONDA_PREFIX")
CORE_R_PACKAGES = [
    "xml2",
    "XML",
    "optparse",
    "GEOquery",
    "minfi",
    "sesame",
    "limma",
    "dmrff",
    "sesameData",
    "sva",
    "variancePartition",
    "pvca",
    "illuminaio",
    "RefFreeEWAS",
    "Biobase",
    "reformulas",
    "ggplot2",
    "plotly",
    "DT",
    "data.table",
    "htmlwidgets",
    "dplyr",
    "stringr",
    "jsonlite",
    "ggrepel",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "planet",
]
EPICV2_R_PACKAGES = [
    "IlluminaHumanMethylationEPICv2manifest",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
]
OPTIONAL_R_PACKAGES = [
    "FlowSorted.Blood.EPIC",
    "FlowSorted.Blood.450k",
    "FlowSorted.CordBlood.EPIC",
    "FlowSorted.CordBlood.450k",
    "FlowSorted.CordBloodCombined.450k",
    "FlowSorted.DLPFC.450k",
    "FlowSorted.Saliva.450k",
    "wateRmelon",
    "methylclock",
    "methylclockData",
]
def add_conda_paths(env: dict) -> dict:
    """Ensure LD_LIBRARY_PATH/PKG_CONFIG_PATH/PATH include conda libs so xml2/xml load correctly."""
    use_conda_libs = (env.get("ILLUMETA_USE_CONDA_LIBS") == "1") or (
        os.environ.get("ILLUMETA_USE_CONDA_LIBS") == "1"
    )
    if not use_conda_libs:
        return env
    prefix = env.get("CONDA_PREFIX") or os.environ.get("CONDA_PREFIX")
    if not prefix or not os.path.isdir(prefix):
        return env
    conda_lib = os.path.join(prefix, "lib")
    conda_pkgconfig = os.path.join(conda_lib, "pkgconfig")
    conda_bin = os.path.join(prefix, "bin")
    def prepend(path_var, new_path):
        cur = env.get(path_var, "")
        parts = cur.split(os.pathsep) if cur else []
        if new_path and os.path.exists(new_path) and new_path not in parts:
            env[path_var] = os.pathsep.join([new_path] + parts if parts else [new_path])
    prepend("LD_LIBRARY_PATH", conda_lib)
    prepend("PKG_CONFIG_PATH", conda_pkgconfig)
    prepend("PATH", conda_bin)
    return env

def log(msg: str):
    """Prints a timestamped info message."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

def log_err(msg: str):
    """Prints a timestamped error message to stderr."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", file=sys.stderr)

def check_r_installation():
    """Checks if R is available in the path."""
    try:
        subprocess.run(["R", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        log_err("Error: R is not installed or not in the PATH.")
        log_err("  - Ubuntu/Debian: sudo apt-get update && sudo apt-get install r-base")
        log_err("  - macOS (Homebrew): brew install --cask r")
        log_err("  - Conda: conda install -c conda-forge r-base")
        log_err("After installation, ensure 'R' and 'Rscript' are in your PATH and re-run illumeta.py.")
        sys.exit(1)

def check_pandoc_installation():
    """Checks if pandoc is available in the path."""
    try:
        subprocess.run(["pandoc", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        log_err("Error: 'pandoc' is not installed or not in the PATH.")
        log_err("Pandoc is required for generating interactive HTML reports.")
        log_err("  - macOS: brew install pandoc")
        log_err("  - Linux: sudo apt-get install pandoc (Debian/Ubuntu) or yum install pandoc (CentOS/RHEL)")
        log_err("  - Windows: choco install pandoc")
        sys.exit(1)

def check_r_package(pkg, env):
    expr = (
        'lib <- Sys.getenv("R_LIBS_USER"); '
        'if (nzchar(lib)) .libPaths(c(lib, .Library, .Library.site)); '
        f'ok <- requireNamespace("{pkg}", quietly=TRUE); '
        'cat(if (isTRUE(ok)) "OK" else "MISSING")'
    )
    try:
        res = subprocess.run(
            ["Rscript", "-e", expr],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except FileNotFoundError:
        return False, "Rscript not found"
    if res.returncode == 0 and "OK" in res.stdout:
        return True, ""
    err = (res.stderr + res.stdout).strip()
    return False, err

def check_r_packages(pkgs, env=None):
    if not pkgs:
        return [], {}
    if env is None:
        env = ensure_r_lib_env(os.environ.copy())
        env = add_conda_paths(env)
    missing = []
    errors = {}
    for pkg in pkgs:
        ok, err = check_r_package(pkg, env)
        if not ok:
            missing.append(pkg)
            if err:
                errors[pkg] = err
    return missing, errors

def missing_r_packages(pkgs, env=None):
    """Returns a list of missing or unloadable R packages."""
    missing, _ = check_r_packages(pkgs, env=env)
    return missing

def format_r_error(err: str) -> str:
    if not err:
        return ""
    lower = err.lower()
    if "segfault" in lower:
        return "R crashed while loading this package (segfault)"
    line = err.splitlines()[0].strip()
    return line[:200]

def ensure_r_lib_env(env):
    """
    Ensure R uses a repo-local, writable library by default for reproducibility.
    To respect an externally-set R_LIBS_USER instead, set ILLUMETA_RESPECT_R_LIBS_USER=1.
    """
    respect_existing = (env.get("ILLUMETA_RESPECT_R_LIBS_USER") == "1") or (
        os.environ.get("ILLUMETA_RESPECT_R_LIBS_USER") == "1"
    )
    existing = env.get("R_LIBS_USER")
    if respect_existing and existing:
        return env

    recorded_lib = read_recorded_r_lib(BASE_DIR)
    default_lib_root = default_r_lib_base(BASE_DIR, env)
    default_lib = recorded_lib or default_r_lib(BASE_DIR, env)
    if recorded_lib:
        default_lib_root = os.path.dirname(recorded_lib.rstrip(os.sep))
    if existing and existing != default_lib and not respect_existing:
        log(f"[*] Overriding R_LIBS_USER={existing} -> {default_lib} (set ILLUMETA_RESPECT_R_LIBS_USER=1 to keep)")

    if BASE_DIR.startswith("/Volumes/") and default_lib_root != os.path.join(BASE_DIR, ".r-lib") and not respect_existing:
        log(f"[*] Project is on external volume; using local R library at {default_lib} (set ILLUMETA_ALLOW_EXTERNAL_LIB=1 to use .r-lib)")

    env["R_LIBS_USER"] = default_lib
    if not os.path.exists(default_lib):
        os.makedirs(default_lib, exist_ok=True)
        log(f"[*] Created default R library at {default_lib}")
    return env

def ensure_r_dependencies():
    """Runs the setup_env.R script to install required R packages and cache data."""
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)
    force_setup = os.environ.get("ILLUMETA_FORCE_SETUP") == "1"
    require_epicv2 = os.environ.get("ILLUMETA_REQUIRE_EPICV2") == "1"
    marker_ok = False
    marker_epicv2_required = None
    if os.path.exists(SETUP_MARKER):
        try:
            with open(SETUP_MARKER) as f:
                content = f.read()
            for line in content.splitlines():
                if line.startswith("epicv2_required="):
                    marker_epicv2_required = line.split("=", 1)[1].strip() == "1"
            if marker_epicv2_required is None:
                marker_epicv2_required = False
            marker_ok = marker_epicv2_required == require_epicv2
        except Exception:
            marker_ok = False
    core_pkgs = CORE_R_PACKAGES + (EPICV2_R_PACKAGES if require_epicv2 else [])
    optional_pkgs = OPTIONAL_R_PACKAGES + ([] if require_epicv2 else EPICV2_R_PACKAGES)
    missing_core = missing_r_packages(core_pkgs, env=env) if marker_ok else []
    missing_optional = missing_r_packages(optional_pkgs, env=env) if marker_ok else []
    if marker_ok and not force_setup and not missing_core:
        if missing_optional:
            log(f"[*] Optional R packages missing (analysis will skip related features): {', '.join(missing_optional)}")
            if any(pkg in EPICV2_R_PACKAGES for pkg in missing_optional):
                log("[*] EPIC v2 support is not installed. Use R 4.4+ and set ILLUMETA_REQUIRE_EPICV2=1 to require it.")
        log("[*] R dependencies already set up (skipping). Set ILLUMETA_FORCE_SETUP=1 to force reinstall.")
        return
    if marker_ok and missing_core and not force_setup:
        log(f"[*] R setup marker found but missing packages: {', '.join(missing_core)}. Re-running setup_env.R ...")

    log("[*] Ensuring R dependencies (this may take a few minutes on first run)...")
    try:
        subprocess.run(["Rscript", SETUP_SCRIPT], check=True, env=env)
        with open(SETUP_MARKER, "w") as f:
            f.write(f"setup completed at {datetime.now().isoformat()}\n")
            f.write(f"epicv2_required={1 if require_epicv2 else 0}\n")
        missing_optional_after = missing_r_packages(optional_pkgs, env=env)
        if missing_optional_after:
            log(f"[*] Optional R packages missing (features will be skipped): {', '.join(missing_optional_after)}")
            if any(pkg in EPICV2_R_PACKAGES for pkg in missing_optional_after):
                log("[*] EPIC v2 support is not installed. Use R 4.4+ and set ILLUMETA_REQUIRE_EPICV2=1 to require it.")
    except subprocess.CalledProcessError as e:
        log_err(f"[!] Error while installing R dependencies: {e}")
        log_err("    Hint: If you see 'library path not writable', set a user library first:")
        log_err("      export R_LIBS_USER=\"$HOME/R/library\" && mkdir -p \"$R_LIBS_USER\"")
        log_err("    Then rerun: Rscript r_scripts/setup_env.R")
        if os.path.exists(SETUP_MARKER):
            try:
                os.remove(SETUP_MARKER)
            except OSError:
                pass
        sys.exit(1)

def run_doctor(args):
    """Checks system and R dependencies without installing anything."""
    check_r_installation()
    if not args.skip_pandoc:
        check_pandoc_installation()

    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)

    log("[*] Checking required R packages...")
    require_epicv2 = os.environ.get("ILLUMETA_REQUIRE_EPICV2") == "1"
    core_pkgs = CORE_R_PACKAGES + (EPICV2_R_PACKAGES if require_epicv2 else [])
    optional_pkgs = OPTIONAL_R_PACKAGES + ([] if require_epicv2 else EPICV2_R_PACKAGES)
    missing_core, core_errors = check_r_packages(core_pkgs, env=env)
    missing_optional, optional_errors = check_r_packages(optional_pkgs, env=env)

    if missing_core:
        log_err("[!] Missing required R packages:")
        saw_segfault = False
        for pkg in missing_core:
            err = format_r_error(core_errors.get(pkg, ""))
            if "segfault" in err.lower():
                saw_segfault = True
            if err:
                log_err(f"  - {pkg}: {err}")
            else:
                log_err(f"  - {pkg}")
        if saw_segfault:
            log_err("    Hint: A package crashed R while loading. Rerun setup with `ILLUMETA_CLEAN_MISMATCHED_RLIB=1`.")
        log_err("    Fix: run `ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R` (or run any IlluMeta command once).")
        sys.exit(1)

    log("[*] Core R packages: OK")
    if missing_optional:
        log("[*] Optional R packages missing (features will be skipped): " + ", ".join(missing_optional))
        for pkg in missing_optional:
            err = format_r_error(optional_errors.get(pkg, ""))
            if err:
                log(f"[*]   - {pkg}: {err}")
        if any(pkg in EPICV2_R_PACKAGES for pkg in missing_optional):
            log("[*] EPIC v2 support is not installed. Use R 4.4+ and set ILLUMETA_REQUIRE_EPICV2=1 to require it.")
    else:
        log("[*] Optional R packages: OK")

def run_download(args):
    """Executes the download step."""
    out_dir = args.out_dir if args.out_dir else os.path.abspath(args.gse_id)
    log(f"[*] Starting download pipeline for {args.gse_id}...")
    log(f"[*] Output directory: {out_dir}")
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    cmd = ["Rscript", DOWNLOAD_SCRIPT, "--gse", args.gse_id, "--out", out_dir]
    if args.platform:
        cmd += ["--platform", args.platform]
    
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)

    try:
        subprocess.run(cmd, check=True, env=env)
        log(f"[*] Download complete. Fill 'primary_group' in: {os.path.join(out_dir, 'configure.tsv')}")
        log("[*] Tip: you can auto-fill it during analysis with --auto-group/--group-column/--group-key.")
    except subprocess.CalledProcessError as e:
        log_err(f"[!] Error during download step: {e}")
        sys.exit(1)


# --- GEO Search (IlluMeta Search) ---
SEARCH_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SEARCH_PLATFORM_ACCESSIONS = ["GPL13534", "GPL21145", "GPL34372"]  # 450k, EPIC(850k), EPIC v2 (950k)
SEARCH_PLATFORM_LABELS = {"GPL13534": "450k", "GPL21145": "850k", "GPL34372": "950k"}
SEARCH_IDAT_REGEX = re.compile(r"(\.idat(\.gz)?|RAW\.tar)", re.IGNORECASE)
SEARCH_TOOL_NAME = "illumeta-search"
SEARCH_ESUMMARY_BATCH_SIZE = 150
SEARCH_SUPPL_RETRY = 3
SEARCH_SUPPL_RETRY_DELAY = 0.5
SEARCH_ORGANISM_FILTER = '"Homo sapiens"[Organism]'
SEARCH_EUTILS_RETRY = 3
SEARCH_EUTILS_RETRY_DELAY = 1.0


def build_search_term(keywords: str) -> str:
    if not keywords.strip():
        raise ValueError("keywords must be a non-empty string")
    platform_query = " OR ".join([f"{p}[Accession]" for p in SEARCH_PLATFORM_ACCESSIONS])
    keyword_clause = f"({keywords})"
    term = f"{keyword_clause} AND (idat[All Fields]) AND ({platform_query}) AND gse[Entry Type] AND {SEARCH_ORGANISM_FILTER}"
    return term


def search_eutils_request(path: str, params: dict, email: str, sleep_s: float):
    if email:
        params["email"] = email
    params["tool"] = SEARCH_TOOL_NAME
    url = f"{SEARCH_EUTILS}/{path}"
    last_exc = None
    for attempt in range(1, SEARCH_EUTILS_RETRY + 1):
        try:
            resp = requests.get(url, params=params, timeout=20)
            resp.raise_for_status()
            time.sleep(sleep_s)
            return resp
        except Exception as exc:
            last_exc = exc
            if attempt < SEARCH_EUTILS_RETRY:
                time.sleep(SEARCH_EUTILS_RETRY_DELAY)
    if last_exc:
        raise last_exc
    raise RuntimeError("Unexpected error contacting NCBI E-utilities.")


def search_gse_ids(term: str, email: str, sleep_s: float, retmax: int):
    log(f"[*] Searching GEO GSE with term: {term}")
    params = {"db": "gds", "term": term, "retmax": str(retmax), "retmode": "xml"}
    resp = search_eutils_request("esearch.fcgi", params, email, sleep_s)
    root = ET.fromstring(resp.text)
    ids = [elem.text for elem in root.findall(".//Id") if elem.text]
    log(f"[*] Found {len(ids)} candidate GSE records")
    return ids


def fetch_search_summaries(ids, email: str, sleep_s: float):
    if not ids:
        return []
    log(f"[*] Fetching summaries for {len(ids)} IDs (batch size: {SEARCH_ESUMMARY_BATCH_SIZE})...")
    out = []
    for start in range(0, len(ids), SEARCH_ESUMMARY_BATCH_SIZE):
        chunk = ids[start : start + SEARCH_ESUMMARY_BATCH_SIZE]
        id_str = ",".join(chunk)
        params = {"db": "gds", "id": id_str, "retmode": "xml"}
        resp = search_eutils_request("esummary.fcgi", params, email, sleep_s)
        root = ET.fromstring(resp.text)
        before = len(out)
        for doc in root.findall(".//DocSum"):
            rec = {}
            for item in doc.findall("Item"):
                name = item.get("Name")
                if name:
                    if len(item):
                        children = [(child.text or "").strip() for child in item.findall("Item")]
                        value = ";".join([c for c in children if c])
                    else:
                        value = (item.text or "").strip()
                    rec[name.lower()] = value

            accession = rec.get("accession")
            if accession and accession.startswith("GSE"):
                raw_gpl = rec.get("gpl", "")
                tokens = [t for t in re.split(r"[;,\s]+", raw_gpl) if t]
                gpls = []
                for tok in tokens:
                    tok = tok.upper()
                    if tok.startswith("GPL"):
                        gpls.append(tok)
                    elif tok.isdigit():
                        gpls.append(f"GPL{tok}")
                    else:
                        gpls.append(tok)
                platform = ";".join(gpls)
                platform_types = sorted({SEARCH_PLATFORM_LABELS[g] for g in gpls if g in SEARCH_PLATFORM_LABELS})
                platform_label = ";".join(platform_types)
                pubmed_ids = re.sub(r"\s+", " ", rec.get("pubmedids", "")).strip()
                out.append(
                    {
                        "gse_id": accession,
                        "title": rec.get("title", ""),
                        "platform": platform,
                        "platform_type": platform_label,
                        "sample_count": rec.get("n_samples", ""),
                        "pubmed_id": pubmed_ids,
                    }
                )
        parsed = len(out) - before
        log(f"    - Parsed {parsed} from batch {start // SEARCH_ESUMMARY_BATCH_SIZE + 1}")
    log(f"[*] Parsed {len(out)} GSE summaries total")
    return out


def geo_suppl_url(gse_id: str) -> str:
    prefix = gse_id[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"


def search_fetch_with_retry(url: str, method: str = "get", headers=None, stream: bool = False):
    last_resp = None
    for attempt in range(1, SEARCH_SUPPL_RETRY + 1):
        try:
            resp = requests.request(
                method=method,
                url=url,
                timeout=15,
                allow_redirects=True,
                headers=headers,
                stream=stream,
            )
            last_resp = resp
            if resp.status_code == 200:
                return resp
        except Exception:
            last_resp = None
        if attempt < SEARCH_SUPPL_RETRY:
            time.sleep(SEARCH_SUPPL_RETRY_DELAY)
    return last_resp


def raw_tar_exists(raw_url: str, gse_id: str):
    ftp_statuses = []
    alt_statuses = []

    def probe(url: str, method: str, headers=None, stream: bool = False):
        try:
            resp = requests.request(method=method, url=url, timeout=15, allow_redirects=True, headers=headers, stream=stream)
            status = resp.status_code
            resp.close()
            return status
        except Exception:
            return None

    for attempt in range(1, SEARCH_SUPPL_RETRY + 1):
        status = probe(raw_url, "head")
        ftp_statuses.append(status)
        if status == 200:
            return True
        if attempt < SEARCH_SUPPL_RETRY:
            time.sleep(SEARCH_SUPPL_RETRY_DELAY)

    range_headers = {"Range": "bytes=0-0"}
    for attempt in range(1, SEARCH_SUPPL_RETRY + 1):
        status = probe(raw_url, "get", headers=range_headers, stream=True)
        ftp_statuses.append(status)
        if status in (200, 206):
            return True
        if attempt < SEARCH_SUPPL_RETRY:
            time.sleep(SEARCH_SUPPL_RETRY_DELAY)

    alt_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file"
    for attempt in range(1, 3):
        status = probe(alt_url, "head")
        alt_statuses.append(status)
        if status == 200:
            return True
        if attempt < 2:
            time.sleep(SEARCH_SUPPL_RETRY_DELAY)

    non_none_alt = [s for s in alt_statuses if s is not None]
    if non_none_alt and all(code == 404 for code in non_none_alt):
        return False

    non_none_ftp = [s for s in ftp_statuses if s is not None]
    if non_none_ftp and all(code == 404 for code in non_none_ftp) and not non_none_alt:
        return None
    return None


def has_idat_in_suppl(gse_id: str) -> str:
    url = geo_suppl_url(gse_id)
    resp = search_fetch_with_retry(url, method="get")
    if resp is None:
        return "error"
    resp_text = resp.text
    status = resp.status_code
    resp.close()
    if status != 200:
        return "error"
    if SEARCH_IDAT_REGEX.search(resp_text):
        return "yes"

    raw_url = f"{url}{gse_id}_RAW.tar"
    raw_status = raw_tar_exists(raw_url, gse_id)
    if raw_status is True:
        return "yes"
    if raw_status is False:
        return "no"
    return "error"


def enrich_with_suppl_check(rows, check_suppl: bool, sleep_s: float):
    if not check_suppl:
        return
    log("[*] Checking GEO supplementary directories for actual IDAT/RAW files...")
    for idx, row in enumerate(rows, start=1):
        gse_id = row["gse_id"]
        sys.stdout.write(f"    [{idx}/{len(rows)}] {gse_id}: ")
        sys.stdout.flush()
        found = has_idat_in_suppl(gse_id)
        row["suppl_has_idat"] = found
        if found == "yes":
            print("YES")
        elif found == "no":
            print("no")
        else:
            print("error")
        time.sleep(sleep_s)


def write_search_tsv(rows, path: str) -> None:
    if not rows:
        log("[!] No results to write.")
        return
    fieldnames = ["gse_id", "title", "platform", "platform_type", "sample_count", "pubmed_id", "suppl_has_idat"]
    for row in rows:
        row.setdefault("suppl_has_idat", "unchecked")
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    log(f"[*] Saved {len(rows)} records to {path}")


def run_search(args):
    try:
        import requests as _requests
    except ImportError:
        log_err("Error: Python package 'requests' is missing.")
        log_err("Install with: python3 -m pip install requests")
        return

    globals()["requests"] = _requests

    if args.email:
        log(f"[*] Using contact email for NCBI requests: {args.email}")

    keywords = (args.keywords or "").strip()
    if not keywords:
        log_err("--keywords is required and cannot be empty.")
        return

    term = build_search_term(keywords)
    ids = search_gse_ids(term, args.email, args.sleep, args.retmax)
    summaries = fetch_search_summaries(ids, args.email, args.sleep)
    enrich_with_suppl_check(summaries, args.check_suppl, args.sleep)
    write_search_tsv(summaries, args.output)
    log("[*] Done.")

def write_failure_summary(base_dir: str, stage: str, code: str, message: str, details=None):
    if not base_dir:
        base_dir = os.getcwd()
    os.makedirs(base_dir, exist_ok=True)
    payload = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "stage": stage,
        "code": code,
        "message": message,
        "details": details or {},
    }
    try:
        with open(os.path.join(base_dir, "failure_summary.json"), "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=True)
        with open(os.path.join(base_dir, "failure_reason.txt"), "w", encoding="utf-8") as handle:
            handle.write(f"{code}: {message}\n")
    except Exception as e:
        log_err(f"[!] Failed to write failure summary: {e}")

def run_analysis(args):
    """Executes the analysis step."""
    # Resolve Config Path
    config_path = None
    if args.config:
        config_path = args.config
    elif args.input_dir:
        config_path = os.path.join(args.input_dir, "configure.tsv")
    
    if not config_path:
        log_err("Error: You must provide either --input-dir or --config.")
        write_failure_summary(args.input_dir or os.getcwd(), "config", "CONFIG_MISSING",
                              "Missing --input-dir/--config.")
        sys.exit(1)

    log(f"[*] Starting analysis pipeline using configuration: {config_path}...")
    
    if not os.path.exists(config_path):
        log_err(f"Error: Configuration file {config_path} not found.")
        write_failure_summary(args.input_dir or os.getcwd(), "config", "CONFIG_NOT_FOUND",
                              f"Configuration file not found: {config_path}")
        sys.exit(1)
    config_base_dir = os.path.dirname(os.path.abspath(config_path))
    default_dir_base = args.input_dir if args.input_dir else config_base_dir
    default_folder = safe_path_component(f"{args.group_test}_vs_{args.group_con}_results", fallback="analysis_results")
    planned_output_dir = args.output or os.path.join(default_dir_base, default_folder)

    if args.beginner_safe:
        # Enforce more conservative defaults for novice safety
        if args.min_total_size < 8:
            args.min_total_size = 8
        if args.delta_beta <= 0:
            args.delta_beta = 0.05

    auto_group_info = None
    auto_group_requested = args.auto_group or args.group_column or args.group_key or args.group_map
    if auto_group_requested:
        try:
            config_path, auto_group_info = auto_group_config(
                config_path=config_path,
                group_con=args.group_con,
                group_test=args.group_test,
                group_column=args.group_column,
                group_key=args.group_key,
                group_map=args.group_map,
                output_path=args.auto_group_output,
                id_column=args.id_column,
                overwrite=args.auto_group_overwrite,
                allow_technical=args.auto_group_allow_technical,
            )
            if auto_group_info and auto_group_info.get("updated"):
                log(f"[*] Auto-group: filled primary_group using {auto_group_info.get('source')} -> {config_path}")
            else:
                log("[*] Auto-group: primary_group already filled; using existing configure.tsv")
        except ValueError as e:
            log_err(f"[!] Auto-group failed: {e}")
            write_failure_summary(planned_output_dir, "auto_group", "AUTO_GROUP_FAILED", str(e))
            sys.exit(1)

    if args.beginner_safe:
        try:
            rows, headers, delim = load_config_rows(config_path)
            if delim != "\t":
                raise ValueError("configure.tsv must be tab-delimited (TSV). CSV/other delimiters are not supported.")
            con_norm = normalize_group_value(args.group_con)
            test_norm = normalize_group_value(args.group_test)
            counts = Counter()
            for row in rows:
                val = (row.get("primary_group") or "").strip()
                if val:
                    counts[normalize_group_value(val)] += 1
            extra = [g for g in counts if g not in {con_norm, test_norm}]
            n_con = counts.get(con_norm, 0)
            n_test = counts.get(test_norm, 0)
            if extra:
                raise ValueError(
                    "Beginner-safe mode requires only two groups (control/test). "
                    f"Found extra group labels: {', '.join(sorted(extra)[:5])}. "
                    "Use --group-map to collapse or specify a cleaner group column/key."
                )
            if n_con < 3 or n_test < 3:
                raise ValueError(
                    f"Beginner-safe mode requires >=3 samples per group. "
                    f"Got {args.group_con}={n_con}, {args.group_test}={n_test}."
                )
        except ValueError as e:
            log_err(f"[!] Beginner-safe check failed: {e}")
            write_failure_summary(planned_output_dir, "beginner_safe", "BEGINNER_SAFE_FAILED", str(e))
            sys.exit(1)
        
    # Determine output directory
    if args.output:
        output_dir = args.output
    else:
        # Default: test_vs_con_results (sanitize for filesystem/locale safety)
        dir_base = args.input_dir if args.input_dir else config_base_dir
        folder_name = f"{args.group_test}_vs_{args.group_con}_results"
        safe_folder = safe_path_component(folder_name, fallback="analysis_results")
        if safe_folder != folder_name:
            log(f"[*] Output folder normalized for filesystem safety: {folder_name} -> {safe_folder}")
        output_dir = os.path.join(dir_base, safe_folder)
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    should_prepare_cross = (
        not args.cross_reactive_list
        and not args.disable_cross_reactive
        and not getattr(args, "unsafe_skip_cross_reactive", False)
    )
    if should_prepare_cross:
        cross_dir = resolve_cross_reactive_dir(config_path, args.config_yaml)
        created, status = ensure_cross_reactive_lists(cross_dir, allow_network=True)
        if status == "downloaded":
            log(f"[*] Downloaded cross-reactive probe lists to {cross_dir}")
        elif status == "built_local":
            log(f"[*] Built cross-reactive probe lists from local sources in {cross_dir}")
        elif status == "exists":
            log(f"[*] Cross-reactive probe lists already present: {cross_dir}")
        elif status == "skip_env":
            log("[!] Cross-reactive list auto-download disabled via ILLUMETA_SKIP_XREACT_DOWNLOAD=1.")
        elif status == "missing_sources":
            log_err("[!] Cross-reactive list auto-download failed (missing sources).")
            log_err(f"    Provide --cross-reactive-list, install maxprobes, or place lists under {cross_dir}.")
    
    cmd = [
        "Rscript", ANALYZE_SCRIPT, 
        "--config", config_path, 
        "--out", output_dir, 
        "--group_con", args.group_con,
        "--group_test", args.group_test,
        "--max_plots", str(args.max_plots),
        "--pval", str(args.pval),
        "--lfc", str(args.lfc),
        "--delta_beta", str(args.delta_beta),
        "--permutations", str(args.permutations),
        "--min_total_size", str(args.min_total_size),
        "--qc_intensity_threshold", str(args.qc_intensity_threshold)
    ]
    if args.beginner_safe:
        cmd.append("--beginner_safe")
    if args.marker_list:
        cmd.extend(["--marker_list", args.marker_list])
    if args.cross_reactive_list:
        cmd.extend(["--cross_reactive_list", args.cross_reactive_list])
    if getattr(args, "unsafe_skip_cross_reactive", False):
        cmd.append("--unsafe-skip-cross-reactive")
    if args.disable_cross_reactive:
        cmd.append("--unsafe-skip-cross-reactive")
    if getattr(args, "sex_mismatch_action", None):
        cmd.extend(["--sex-mismatch-action", args.sex_mismatch_action])
    elif args.sex_check_action:
        cmd.extend(["--sex-mismatch-action", args.sex_check_action])
    if args.sex_check_column:
        cmd.extend(["--sex_check_column", args.sex_check_column])
    if args.disable_sex_check:
        cmd.append("--disable_sex_check")
    if args.batch_column:
        cmd.extend(["--batch_column", args.batch_column])
    if args.batch_method:
        cmd.extend(["--batch_method", args.batch_method])
    if getattr(args, "tier3_min_total_n", None) is not None:
        cmd.extend(["--tier3-min-total-n", str(args.tier3_min_total_n)])
    if getattr(args, "tier3_min_per_group_per_stratum", None) is not None:
        cmd.extend(["--tier3-min-per-group-per-stratum", str(args.tier3_min_per_group_per_stratum)])
    if getattr(args, "tier3_on_fail", None):
        cmd.extend(["--tier3-on-fail", args.tier3_on_fail])
    if getattr(args, "auto_covariates_enabled", None):
        cmd.extend(["--auto-covariates-enabled", args.auto_covariates_enabled])
    if getattr(args, "auto_covariates_exclude_group_associated", None):
        cmd.extend(["--auto-covariates-exclude-group-associated", args.auto_covariates_exclude_group_associated])
    if getattr(args, "auto_covariates_group_assoc_p_threshold", None) is not None:
        cmd.extend(["--auto-covariates-group-assoc-p-threshold", str(args.auto_covariates_group_assoc_p_threshold)])
    if getattr(args, "auto_covariates_max_cor", None) is not None:
        cmd.extend(["--auto-covariates-max-cor", str(args.auto_covariates_max_cor)])
    if getattr(args, "cell_adjustment_on_high_eta2", None):
        cmd.extend(["--cell-adjustment-on-high-eta2", args.cell_adjustment_on_high_eta2])
    idat_dir = None
    if args.idat_dir:
        idat_dir = args.idat_dir
        if not os.path.isabs(idat_dir):
            idat_dir = os.path.join(os.path.dirname(os.path.abspath(config_path)), idat_dir)
        if not os.path.isdir(idat_dir):
            log_err(f"[!] IDAT directory not found: {idat_dir}")
            write_failure_summary(output_dir, "idat", "IDAT_DIR_MISSING", f"IDAT directory not found: {idat_dir}")
            sys.exit(1)
        cmd.extend(["--idat_dir", idat_dir])

    try:
        preflight = preflight_analysis(
            config_path=config_path,
            idat_dir=idat_dir,
            group_con=args.group_con,
            group_test=args.group_test,
            min_total_size=args.min_total_size,
            id_column=args.id_column,
            drop_missing_idat=not (args.fail_on_missing_idat or args.keep_missing_idat),
        )
        if preflight.get("config_path") and preflight["config_path"] != config_path:
            config_path = preflight["config_path"]
            config_idx = cmd.index("--config") + 1
            cmd[config_idx] = config_path
            log(f"[*] Using filtered config: {config_path}")
        log(f"[*] Preflight OK: samples={preflight['sample_count']}, "
            f"{args.group_con}={preflight['group_con']}, {args.group_test}={preflight['group_test']}")
        if preflight.get("warnings"):
            log("[!] Preflight warnings:")
            for warn in preflight["warnings"]:
                log(f"    - {warn}")
        preflight_report_path = os.path.join(output_dir, "preflight_report.json")
        with open(preflight_report_path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "config_path": config_path,
                    "idat_dir": preflight.get("idat_dir"),
                    "sample_count_raw": preflight.get("sample_count_raw"),
                    "sample_count": preflight.get("sample_count"),
                    "group_counts": {
                        "control": preflight.get("group_con"),
                        "test": preflight.get("group_test"),
                    },
                    "missing_idat_pairs": preflight.get("missing_idat_pairs"),
                    "missing_idat_preview": preflight.get("missing_idat_preview"),
                    "group_labels": preflight.get("group_labels", {}),
                    "batch_candidates": preflight.get("batch_candidates", []),
                    "warnings": preflight.get("warnings", []),
                    "column_profile": preflight.get("column_profile", []),
                    "auto_group": auto_group_info,
                },
                handle,
                indent=2,
                ensure_ascii=True,
            )
    except ValueError as e:
        log_err(f"[!] Preflight check failed: {e}")
        write_failure_summary(output_dir, "preflight", "PREFLIGHT_FAILED", str(e))
        sys.exit(1)
    if args.force_idat:
        cmd.append("--force_idat")
    if args.disable_auto_covariates:
        cmd.append("--disable_auto_covariates")
    if args.disable_sva:
        cmd.append("--disable_sva")
    if args.include_covariates:
        cmd.extend(["--include_covariates", args.include_covariates])
    if args.include_clock_covariates:
        cmd.append("--include_clock_covariates")
    if args.tissue:
        cmd.extend(["--tissue", args.tissue])
    if args.cell_reference:
        cmd.extend(["--cell_reference", args.cell_reference])
    if args.cell_reference_platform:
        cmd.extend(["--cell_reference_platform", args.cell_reference_platform])
    if args.positive_controls:
        cmd.extend(["--positive_controls", args.positive_controls])
    if args.skip_sesame:
        cmd.append("--skip-sesame")
    if args.sesame_typeinorm:
        cmd.append("--sesame_typeinorm")
    if args.config_yaml:
        cmd.extend(["--config_yaml", args.config_yaml])
    if args.preset:
        cmd.extend(["--preset", args.preset])
    if args.qc_detection_p_threshold is not None:
        cmd.extend(["--qc_detection_p_threshold", str(args.qc_detection_p_threshold)])
    if args.qc_sample_fail_frac is not None:
        cmd.extend(["--qc_sample_fail_frac", str(args.qc_sample_fail_frac)])
    if args.dmr_min_cpgs is not None:
        cmd.extend(["--dmr_min_cpgs", str(args.dmr_min_cpgs)])
    if args.dmr_maxgap is not None:
        cmd.extend(["--dmr_maxgap", str(args.dmr_maxgap)])
    if args.beginner_safe_delta_beta is not None:
        cmd.extend(["--beginner_safe_delta_beta", str(args.beginner_safe_delta_beta)])
    if args.logit_offset is not None:
        cmd.extend(["--logit_offset", str(args.logit_offset)])
    if args.vp_top:
        cmd.extend(["--vp_top", str(args.vp_top)])
    if args.id_column:
        cmd.extend(["--id_column", args.id_column])
    
    # Handle custom temp directory
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)
    # Sesame thread safety: default to single-thread unless explicitly disabled
    if env.get("ILLUMETA_SESAME_SINGLE_THREAD", "1") == "1":
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[key] = "1"
        log("[*] Sesame: forcing single-thread BLAS/OMP to avoid pthread errors.")
    if args.tmp_dir:
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir)
        env["TMPDIR"] = os.path.abspath(args.tmp_dir)
        log(f"[*] Using custom temporary directory: {env['TMPDIR']}")

    try:
        subprocess.run(cmd, check=True, env=env)
        log(f"[*] Analysis complete. Results are in: {output_dir}")
        
        # Generate Dashboard
        generate_dashboard(output_dir, args.group_test, args.group_con)
        
    except subprocess.CalledProcessError as e:
        log_err(f"[!] Error during analysis step: {e}")
        write_failure_summary(
            output_dir,
            "analysis",
            "R_SCRIPT_FAILED",
            "R analysis command failed",
            details={"returncode": e.returncode, "cmd": getattr(e, "cmd", cmd)},
        )
        sys.exit(1)

def safe_int(val):
    try:
        return int(val)
    except (ValueError, TypeError):
        return 0

def safe_float(val):
    if val is None:
        return None
    if isinstance(val, str) and val.strip().upper() in {"", "NA", "NAN"}:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None

def _as_bool(val):
    if isinstance(val, bool):
        return val
    if val is None:
        return False
    if isinstance(val, (int, float)):
        return bool(val)
    if isinstance(val, str):
        val = val.strip().lower()
        if val in {"true", "t", "1", "yes", "y"}:
            return True
        if val in {"false", "f", "0", "no", "n", ""}:
            return False
    return False

def collect_dashboard_warnings(stats, analysis_params, qc_summary, cell_summary=None, cell_assoc=None):
    warnings = []
    if analysis_params:
        unsafe_skip = _as_bool(analysis_params.get("unsafe_skip_cross_reactive"))
        cross_active = _as_bool(analysis_params.get("cross_reactive_active")) or _as_bool(analysis_params.get("cross_reactive_enabled"))
        cross_count = safe_int(analysis_params.get("cross_reactive_count", 0))
        if unsafe_skip or not cross_active:
            warnings.append("Cross-reactive probe filtering was skipped (unsafe).")
        elif cross_active and cross_count == 0:
            warnings.append("Cross-reactive probe list missing or empty; filtering may not have been applied.")

        sex_action = (analysis_params.get("sex_check_action") or "").lower()
        sex_mismatch = safe_int(analysis_params.get("sex_mismatch_count", 0))
        if sex_mismatch > 0:
            if sex_action == "drop":
                warnings.append(f"Sex mismatches detected ({sex_mismatch}); mismatched samples were dropped.")
            elif sex_action == "ignore":
                warnings.append(f"Sex mismatches detected ({sex_mismatch}); mismatches were ignored.")

        auto_cov_enabled = _as_bool(analysis_params.get("auto_covariates_enabled"))
        if auto_cov_enabled:
            warnings.append("Auto-selected covariates may include mediators; verify biological plausibility before interpretation.")
        sample_tier = (analysis_params.get("crf_sample_tier") or "").lower()
        if sample_tier == "minimal":
            warnings.append("CRF tier: minimal sample size; robustness assessment is exploratory only.")
        elif sample_tier == "small":
            warnings.append("CRF tier: small sample size; validate findings independently.")

    if stats:
        primary_mode = (stats.get("primary_result_mode") or "").lower()
        if primary_mode == "tier3_low_power" or _as_bool(stats.get("primary_tier3_low_power")):
            warnings.append("Tier3 stratified/meta-analysis flagged low power; interpret with extreme caution.")
        if primary_mode == "tier3_ineligible" or _as_bool(stats.get("primary_tier3_ineligible")):
            warnings.append("Tier3 confounding detected but eligibility failed; stratified/meta-analysis was not run.")
        no_signal_flag = _as_bool(stats.get("primary_no_signal"))
        if not no_signal_flag:
            inter_native_total = safe_int(stats.get("intersect_native_up")) + safe_int(stats.get("intersect_native_down"))
            inter_total = safe_int(stats.get("intersect_up")) + safe_int(stats.get("intersect_down"))
            minfi_total = safe_int(stats.get("minfi_up")) + safe_int(stats.get("minfi_down"))
            sesame_total = safe_int(stats.get("sesame_up")) + safe_int(stats.get("sesame_down"))
            sesame_native_total = safe_int(stats.get("sesame_native_up")) + safe_int(stats.get("sesame_native_down"))
            no_signal_flag = (inter_native_total == 0 and inter_total == 0 and minfi_total == 0 and
                              sesame_total == 0 and sesame_native_total == 0)
        if no_signal_flag:
            warnings.append("No significant DMPs detected at configured thresholds; likely underpowered or overly stringent.")

    if cell_summary:
        methods = [row.get("Method", "") for row in cell_summary if row.get("Method")]
        if any(m for m in methods if "reffree" not in m.lower()):
            warnings.append("Reference-based cell-type deconvolution was used; disease states can bias estimates.")

    if cell_assoc and analysis_params:
        threshold = safe_float(analysis_params.get("cell_confound_eta2_threshold"))
        if threshold is None:
            threshold = 0.5
        try:
            for row in cell_assoc:
                eta = safe_float(row.get("Eta_Squared"))
                if eta is not None and eta > threshold:
                    warnings.append(f"Cell composition strongly associated with group (Eta^2>{threshold}); over-correction risk.")
                    break
        except Exception:
            pass

    return warnings

def _detect_delimiter(line: str) -> str:
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    if ";" in line:
        return ";"
    return ""

def _read_probe_ids(path, col_candidates=None):
    if not os.path.exists(path):
        return []
    col_candidates = col_candidates or ["TargetID", "CpG", "IlmnID", "Probe", "probe", "ID", "Name", "X", "PROBE"]
    ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
        if not first:
            return []
        delim = _detect_delimiter(first)
        handle.seek(0)
        if delim:
            reader = csv.reader(handle, delimiter=delim)
            rows = list(reader)
            if not rows:
                return []
            header = rows[0]
            header_map = {h.strip(): idx for idx, h in enumerate(header)}
            use_idx = None
            for cand in col_candidates:
                if cand in header_map:
                    use_idx = header_map[cand]
                    break
            start_row = 1
            if use_idx is None:
                use_idx = 0
                start_row = 0
            for row in rows[start_row:]:
                if not row:
                    continue
                if use_idx >= len(row):
                    continue
                val = row[use_idx].strip()
                if not val or val.startswith("#"):
                    continue
                ids.append(val)
        else:
            for line in handle:
                val = line.strip()
                if not val or val.startswith("#"):
                    continue
                ids.append(val.split()[0])
    return ids

def _filter_probe_ids(ids):
    out = []
    for val in ids:
        v = val.strip()
        if not v:
            continue
        if v.lower().startswith(("cg", "ch")):
            out.append(v)
    return sorted(set(out))

def _download_file(url, dest_path):
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    req = _urllib_request.Request(url, headers={"User-Agent": "IlluMeta/1.0"})
    with _urllib_request.urlopen(req, timeout=30) as resp:
        data = resp.read()
    with open(dest_path, "wb") as handle:
        handle.write(data)
    return dest_path

def _try_download(urls, dest_path):
    last_err = None
    for url in urls:
        try:
            return _download_file(url, dest_path), url
        except (URLError, HTTPError, TimeoutError) as err:
            last_err = err
            continue
    if last_err:
        raise last_err
    raise URLError("No URLs provided")

def _read_yaml_scalar(value):
    if value is None:
        return ""
    v = str(value).strip()
    if not v or v.lower() in ("null", "none", "~"):
        return ""
    if len(v) >= 2 and v[0] == v[-1] and v[0] in ("'", '"'):
        v = v[1:-1].strip()
    return v

def _read_cross_reactive_local_dir(yaml_path):
    if not yaml_path or not os.path.exists(yaml_path):
        return ""
    try:
        with open(yaml_path, "r", encoding="utf-8") as handle:
            lines = handle.readlines()
    except Exception:
        return ""
    block_indent = None
    for raw in lines:
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        stripped = line.strip()
        if block_indent is None:
            if stripped.startswith("cross_reactive:"):
                block_indent = indent
            continue
        if indent <= block_indent:
            block_indent = None
            continue
        if stripped.startswith("local_dir:"):
            _, value = stripped.split(":", 1)
            return _read_yaml_scalar(value)
    return ""

def _file_has_min_lines(path, min_lines=10):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for idx, _ in enumerate(handle, start=1):
                if idx >= min_lines:
                    return True
    except Exception:
        return False
    return False

def _build_cross_reactive_450k(chen_path, benton_path):
    ids = []
    ids.extend(_read_probe_ids(chen_path, col_candidates=["TargetID", "CpG", "IlmnID", "probe", "ID"]))
    ids.extend(_read_probe_ids(benton_path, col_candidates=["V1", "CpG", "IlmnID", "probe", "ID"]))
    return _filter_probe_ids(ids)

def _build_cross_reactive_epic(pidsley_path, mccartney_s2_path, mccartney_s3_path):
    ids = []
    ids.extend(_read_probe_ids(pidsley_path, col_candidates=["X", "CpG", "IlmnID", "probe", "ID"]))
    ids.extend(_read_probe_ids(mccartney_s2_path, col_candidates=["CpG", "IlmnID", "probe", "ID"]))
    ids.extend(_read_probe_ids(mccartney_s3_path, col_candidates=["CpG", "IlmnID", "probe", "ID"]))
    return _filter_probe_ids(ids)

def ensure_cross_reactive_lists(dest_dir, allow_network=True):
    if os.environ.get("ILLUMETA_SKIP_XREACT_DOWNLOAD", "").strip() == "1":
        return False, "skip_env"
    os.makedirs(dest_dir, exist_ok=True)
    out_450k = os.path.join(dest_dir, "cross_reactive_450K.tsv")
    out_epic = os.path.join(dest_dir, "cross_reactive_EPIC.tsv")
    out_epicv2 = os.path.join(dest_dir, "cross_reactive_EPICv2.tsv")
    have_main = (
        os.path.exists(out_450k)
        and _file_has_min_lines(out_450k)
        and os.path.exists(out_epic)
        and _file_has_min_lines(out_epic)
    )
    need_epicv2 = not (os.path.exists(out_epicv2) and _file_has_min_lines(out_epicv2))
    if have_main and not need_epicv2:
        return False, "exists"

    sources_dir = os.path.join(dest_dir, "_sources")
    os.makedirs(sources_dir, exist_ok=True)

    raw_base = "https://raw.githubusercontent.com/markgene/maxprobes/master/inst/extdata"
    raw_fallback = "https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master"
    urls = {
        "chen": [
            f"{raw_base}/48639-non-specific-probes-Illumina450k.csv",
            f"{raw_fallback}/48639-non-specific-probes-Illumina450k.csv",
        ],
        "benton": [
            f"{raw_base}/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt",
            f"{raw_fallback}/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt",
        ],
        "pidsley": [
            f"{raw_base}/13059_2016_1066_MOESM1_ESM.csv",
            f"{raw_fallback}/EPIC/13059_2016_1066_MOESM1_ESM.csv",
        ],
        "mccartney_s2": [
            f"{raw_base}/1-s2.0-S221359601630071X-mmc2.txt",
        ],
        "mccartney_s3": [
            f"{raw_base}/1-s2.0-S221359601630071X-mmc3.txt",
        ],
    }

    downloaded = {}
    missing = []
    downloaded_any = False
    for key, url_list in urls.items():
        dest = os.path.join(sources_dir, os.path.basename(url_list[0]))
        meta = {"path": dest, "url": "", "downloaded": False}
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            downloaded[key] = meta
            continue
        if allow_network:
            try:
                path, used_url = _try_download(url_list, dest)
                meta["path"] = path
                meta["url"] = used_url
                meta["downloaded"] = True
                downloaded_any = True
                downloaded[key] = meta
                continue
            except Exception:
                pass
        meta["path"] = ""
        downloaded[key] = meta
        missing.append(key)

    if missing:
        return False, "missing_sources"

    probes_450k = _build_cross_reactive_450k(downloaded["chen"]["path"], downloaded["benton"]["path"])
    probes_epic = _build_cross_reactive_epic(
        downloaded["pidsley"]["path"], downloaded["mccartney_s2"]["path"], downloaded["mccartney_s3"]["path"]
    )

    def _write_list(path, probes):
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("CpG\n")
            for p in probes:
                handle.write(f"{p}\n")

    _write_list(out_450k, probes_450k)
    _write_list(out_epic, probes_epic)
    if not os.path.exists(out_epicv2) or not _file_has_min_lines(out_epicv2):
        _write_list(out_epicv2, probes_epic)
    meta_path = os.path.join(dest_dir, "cross_reactive_sources.json")
    with open(meta_path, "w", encoding="utf-8") as handle:
        json.dump({"sources": downloaded, "created": datetime.utcnow().isoformat() + "Z"}, handle, indent=2)
    status = "downloaded" if downloaded_any else "built_local"
    return True, status

def resolve_cross_reactive_dir(config_path, config_yaml_path=None):
    config_dir = os.path.dirname(os.path.abspath(config_path))
    yaml_path = config_yaml_path or os.path.join(config_dir, "config.yaml")
    local_dir = _read_cross_reactive_local_dir(yaml_path)
    if not local_dir:
        local_dir = os.path.join(config_dir, "references", "probe_blacklists")
    if not os.path.isabs(local_dir):
        local_dir = os.path.join(config_dir, local_dir)
    return local_dir

def generate_dashboard(output_dir, group_test, group_con):
    """Generates a beautiful HTML dashboard to navigate results."""
    # Escape user-supplied group names to prevent XSS in generated HTML
    group_test = html.escape(str(group_test))
    group_con = html.escape(str(group_con))
    results_folder_name = os.path.basename(os.path.normpath(output_dir))
    parent_dir = os.path.dirname(os.path.normpath(output_dir))
    dashboard_filename = f"{results_folder_name}_index.html"
    dashboard_path = os.path.join(parent_dir, dashboard_filename)
    
    # Reordered: Consensus + Native/Strict pipelines
    pipeline_defs = [
        ("Intersection_Native", "Intersection (Native · High-confidence)", "intersection"),
        ("Intersection", "Intersection (Strict)", "intersection"),
        ("Minfi", "Minfi (Noob)", "pipeline"),
        ("Sesame", "Sesame (Strict)", "pipeline"),
        ("Sesame_Native", "Sesame (Native)", "pipeline"),
    ]
    
    # Load Summary Statistics
    summary_path = os.path.join(output_dir, "summary.json")
    stats = {}
    try:
        with open(summary_path, "r") as f:
            stats = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        log_err("[!] Warning: summary.json not found or invalid. Dashboard will lack stats.")
        # Default empty stats
        keys = [
            "n_con", "n_test", "minfi_up", "minfi_down",
            "sesame_up", "sesame_down", "sesame_native_up", "sesame_native_down",
            "intersect_up", "intersect_down", "intersect_native_up", "intersect_native_down"
        ]
        stats = {k: 0 for k in keys}

    intersection_sections = [
        ("Consensus Highlights", [
            ("_Consensus_DMPs.html", "Consensus DMPs Table", "CpGs significant in both pipelines with consistent direction.", "TABLE"),
            ("_Consensus_DMPs.csv", "Consensus DMPs CSV", "Consensus DMP list (CSV).", "CSV"),
        ]),
        ("Agreement Checks", [
            ("_LogFC_Concordance.html", "Pipeline Concordance", "Minfi vs Sesame logFC concordance.", "PLOT"),
            ("_Significant_Overlap.html", "Pipeline Overlap", "Significant DMP counts and overlap.", "PLOT"),
        ]),
    ]

    pipeline_sections_base = [
        ("Primary Results", [
            ("_Volcano.html", "Volcano Plot", "Visualizes significant changes (LogFC vs P-value).", "PLOT"),
            ("_Manhattan.html", "Manhattan Plot", "Genomic distribution of methylation changes.", "PLOT"),
            ("_Top100_Heatmap.html", "Heatmap (Top 100)", "Clustering of the top 100 most significant CpGs.", "PLOT"),
            ("_Top_DMPs.html", "Top DMPs Table", "Searchable table of differentially methylated probes.", "TABLE"),
        ]),
        ("Region-Level Results", [
            ("_DMR_Volcano.html", "DMR Volcano Plot", "Visualizes significant DMRs (Est. Diff vs P-value).", "PLOT"),
            ("_DMR_Manhattan.html", "DMR Manhattan Plot", "Genomic distribution of DMRs.", "PLOT"),
            ("_Top_DMRs_Heatmap.html", "DMR Heatmap (Top 50)", "Average methylation levels of top 50 DMRs.", "PLOT"),
            ("_DMRs_Table.html", "DMRs Table", "Searchable table of differentially methylated regions.", "TABLE"),
        ]),
        ("Quality & Diagnostics", [
            ("_PCA_Before.html", "PCA (Raw)", "Principal Component Analysis before correction.", "PLOT"),
            ("_PCA_After_Correction.html", "PCA (Corrected)", "PCA after removing detected batch effects.", "PLOT"),
            ("_Sample_Clustering_Distance.html", "Sample Clustering", "Euclidean distance matrix of samples.", "PLOT"),
            ("_QQPlot.html", "Q-Q Plot", "Check for genomic inflation and systematic bias.", "PLOT"),
            ("_PVCA.html", "PVCA (Raw)", "Variance explained by group/batch/covariates.", "PLOT"),
            ("_AfterCorrection_PVCA.html", "PVCA (Corrected)", "Variance explained after correction.", "PLOT"),
        ]),
        ("Batch & Covariates", [
            ("_Batch_Method_Comparison.html", "Correction Method", "Batch residual vs group signal by method.", "PLOT"),
            ("_BatchMethodComparison.csv", "Correction Method (CSV)", "Batch correction scoring table (CSV).", "CSV"),
            ("_Batch_Evaluation_Before.html", "Batch Eval (Before)", "Covariate association with PCs (Raw).", "PLOT"),
            ("_Batch_Evaluation_After.html", "Batch Eval (After)", "Covariate association with PCs (Corrected).", "PLOT"),
            ("_AutoCovariates.csv", "Auto Covariates (CSV)", "Auto-selected covariate candidates with PC association P-values.", "CSV"),
            ("_DroppedCovariates.csv", "Dropped Covariates (CSV)", "Covariates removed due to confounding/instability.", "CSV"),
            ("_configure_with_clocks.tsv", "Clock-Merged Config (TSV)", "Config with clock covariates merged (if enabled).", "TSV"),
        ]),
        ("Clocks & Age", [
            ("_Epigenetic_Age_methylclock.html", "Epigenetic Age (methylclock)", "Clock vs chronological age (if provided).", "PLOT"),
            ("_Epigenetic_Age_methylclock.csv", "Epigenetic Age Table", "DNAm clock predictions (methylclock).", "CSV"),
            ("_Placental_Age_planet_scatter.html", "Placental Age Scatter", "RPC vs CPC gestational age (Placenta).", "PLOT"),
            ("_Placental_Age_planet.csv", "Placental Age (planet)", "RPC/CPC gestational age predictions (Placenta).", "CSV"),
        ]),
    ]

    pipeline_sections = {}
    def build_tier3_section(pipe_id):
        suffix = "_Tier3_Primary"
        anchor = os.path.join(output_dir, f"{pipe_id}{suffix}_DMPs.html")
        if not os.path.exists(anchor):
            return None
        return (
            "Tier3 Primary Results",
            [
                (f"{suffix}_Volcano.html", "Tier3 Volcano Plot", "Stratified/meta-analysis volcano plot.", "PLOT"),
                (f"{suffix}_Manhattan.html", "Tier3 Manhattan Plot", "Stratified/meta-analysis Manhattan plot.", "PLOT"),
                (f"{suffix}_QQPlot.html", "Tier3 Q-Q Plot", "Genomic inflation for Tier3 primary results.", "PLOT"),
                (f"{suffix}_DMPs.html", "Tier3 DMPs Table", "Tier3 primary DMP table (CSV/HTML).", "TABLE"),
                (f"{suffix}_Notes.txt", "Tier3 Primary Notes", "Why Tier3 outputs are primary for this dataset.", "TXT"),
            ],
        )
    for pipe_id, _, kind in pipeline_defs:
        if kind == "intersection":
            pipeline_sections[pipe_id] = intersection_sections
        else:
            tier3_section = build_tier3_section(pipe_id)
            sections = list(pipeline_sections_base)
            if tier3_section:
                sections = [tier3_section] + sections
            pipeline_sections[pipe_id] = sections
    
    def load_metrics(pipe_name):
        path = os.path.join(output_dir, f"{pipe_name}_Metrics.csv")
        metrics = {}
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        metrics[row["metric"]] = row["value"]
            except Exception:
                pass
        return metrics

    def build_correction_guidance(metrics):
        if not metrics:
            return None
        lam = safe_float(metrics.get("lambda"))
        batch_after = safe_int(metrics.get("batch_sig_p_lt_0.05_after")) if metrics.get("batch_sig_p_lt_0.05_after") is not None else None
        group_after = safe_float(metrics.get("group_min_p_after"))
        issues = []
        if lam is not None:
            if lam > 1.2:
                issues.append("inflation high")
            elif lam < 0.9:
                issues.append("overcorrection risk")
        if batch_after is not None and batch_after > 0:
            issues.append("residual batch signal")
        if group_after is not None and group_after > 0.1:
            issues.append("group signal weak")
        details = []
        if lam is not None:
            details.append(f"λ={lam:.3f}")
        if batch_after is not None:
            details.append(f"batch p<0.05 after={batch_after}")
        if group_after is not None:
            details.append(f"group min p after={group_after:.3g}")
        detail_txt = ", ".join(details)
        if not issues:
            if detail_txt:
                return f"Correction looks adequate ({detail_txt})."
            return "Correction looks adequate; manual adjustment usually not needed."
        issue_txt = ", ".join(issues)
        if detail_txt:
            return f"Manual review recommended: {issue_txt} ({detail_txt})."
        return f"Manual review recommended: {issue_txt}."

    def load_analysis_params():
        path = os.path.join(output_dir, "analysis_parameters.json")
        try:
            with open(path, "r") as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def load_qc_summary():
        path = os.path.join(output_dir, "QC_Summary.csv")
        metrics = {}
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        metrics[row["metric"]] = row["value"]
            except Exception:
                pass
        return metrics

    def load_cell_deconv_summary():
        path = os.path.join(output_dir, "Cell_Deconvolution_Summary.csv")
        rows = []
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    rows = [row for row in reader]
            except Exception:
                rows = []
        return rows

    def load_cell_assoc():
        for fname in ("Cell_vs_Group_Association.csv", "Cell_Group_Association.csv"):
            path = os.path.join(output_dir, fname)
            if os.path.exists(path):
                try:
                    with open(path, newline="") as f:
                        reader = csv.DictReader(f)
                        return [row for row in reader]
                except Exception:
                    return []
        return []

    def load_drop_reasons(pipe_name):
        path = os.path.join(output_dir, f"{pipe_name}_DroppedCovariates.csv")
        reasons = {}
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        reason = (row.get("Reason") or "").strip()
                        if not reason or reason == "NA":
                            continue
                        reasons[reason] = reasons.get(reason, 0) + 1
            except Exception:
                pass
        return reasons

    def load_csv_rows(path):
        if not os.path.exists(path):
            return []
        try:
            with open(path, newline="") as f:
                reader = csv.DictReader(f)
                return [row for row in reader]
        except Exception:
            return []

    def load_crf_sample_tier():
        rows = load_csv_rows(os.path.join(output_dir, "CRF_Sample_Tier.csv"))
        return rows[0] if rows else {}

    def load_crf_summary_file(name, primary_branch):
        path = os.path.join(output_dir, name)
        if os.path.exists(path):
            return load_csv_rows(path)
        if primary_branch:
            prefixed = os.path.join(output_dir, f"{primary_branch}_{name}")
            return load_csv_rows(prefixed)
        return []

    def summarize_mmc(rows):
        if not rows:
            return {}
        def _top_k(row):
            return safe_int(row.get("top_k")) if row.get("top_k") is not None else 0
        row = max(rows, key=_top_k)
        core = safe_int(row.get("core"))
        total = safe_int(row.get("total_unique"))
        core_pct = (core / total) if total else None
        spearman = safe_float(row.get("spearman_mean")) if row.get("spearman_mean") is not None else safe_float(row.get("spearman_min"))
        return {"core_pct": core_pct, "spearman": spearman, "top_k": safe_int(row.get("top_k")), "methods": safe_int(row.get("methods"))}

    def summarize_sss(rows):
        if not rows:
            return {}
        def _top_k(row):
            return safe_int(row.get("top_k")) if row.get("top_k") is not None else 0
        row = max(rows, key=_top_k)
        overlap = safe_float(row.get("overlap_mean_corr"))
        if overlap is None:
            overlap = safe_float(row.get("overlap_mean_raw"))
        sign = safe_float(row.get("sign_mean_corr"))
        if sign is None:
            sign = safe_float(row.get("sign_mean_raw"))
        return {"overlap": overlap, "sign": sign, "top_k": safe_int(row.get("top_k"))}

    def summarize_ncs(rows):
        if not rows:
            return {}
        preferred = [r for r in rows if (r.get("stage") or "").lower() == "corrected"] or rows
        preferred = [r for r in preferred if (r.get("type") or "").lower() == "snp"] or preferred
        def _n(row):
            nval = safe_int(row.get("n"))
            return nval if nval is not None else 0
        row = max(preferred, key=_n)
        return {
            "lambda": safe_float(row.get("lambda")),
            "sig_rate": safe_float(row.get("sig_rate")),
            "lambda_ci_low": safe_float(row.get("lambda_ci_low")),
            "lambda_ci_high": safe_float(row.get("lambda_ci_high")),
            "stage": row.get("stage"),
            "type": row.get("type"),
        }

    def parse_gene_list(val):
        if val is None:
            return []
        if not isinstance(val, str):
            val = str(val)
        val = val.strip()
        if not val or val.lower() in ("nan", "na", "none"):
            return []
        parts = re.split(r"[;,|/]", val)
        out = []
        for p in parts:
            p = p.strip()
            if not p or p.lower() in ("nan", "na", "none"):
                continue
            out.append(p)
        return out

    def summarize_consensus(path):
        if not os.path.exists(path):
            return {}
        gene_counts = {}
        gene_sign = {}
        abs_deltas = []
        max_abs = None
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene_val = row.get("Gene") or row.get("UCSC_RefGene_Name") or row.get("gene") or ""
                genes = parse_gene_list(gene_val)
                logfc = safe_float(row.get("logFC_mean"))
                if logfc is None:
                    logfc = safe_float(row.get("logFC.Minfi"))
                for g in genes:
                    gene_counts[g] = gene_counts.get(g, 0) + 1
                    if logfc is not None:
                        gene_sign[g] = gene_sign.get(g, 0.0) + (1.0 if logfc >= 0 else -1.0)

                deltas = []
                for key in ("Delta_Beta", "Delta_Beta.Minfi", "Delta_Beta.Sesame", "delta_beta", "delta_beta_mean"):
                    val = safe_float(row.get(key))
                    if val is not None:
                        deltas.append(val)
                if deltas:
                    mean_delta = sum(deltas) / len(deltas)
                    abs_val = abs(mean_delta)
                    abs_deltas.append(abs_val)
                    if max_abs is None or abs_val > max_abs:
                        max_abs = abs_val
        return {"gene_counts": gene_counts, "gene_sign": gene_sign, "abs_deltas": abs_deltas, "max_abs": max_abs}

    def effect_label(median_abs):
        if median_abs is None:
            return ""
        if median_abs < 0.05:
            return "small"
        if median_abs < 0.15:
            return "moderate"
        return "large"

    def classify_warnings(items):
        buckets = {"critical": [], "important": [], "info": []}
        for msg in items:
            lower = msg.lower()
            if "unsafe" in lower or "skipped (unsafe)" in lower or "critical" in lower:
                buckets["critical"].append(msg)
            elif "exploratory" in lower or "minimal" in lower or "low power" in lower or "inflation" in lower or "over-correction" in lower or "overcorrection" in lower:
                buckets["important"].append(msg)
            elif "mismatch" in lower or "confounding" in lower:
                buckets["important"].append(msg)
            else:
                buckets["info"].append(msg)
        return buckets

    def compute_verdict(tier, min_per_group, mmc_core_pct, sss_overlap, warn_count, critical_count):
        tier = (tier or "").lower()
        if critical_count > 0:
            return "EXPLORATORY"
        if tier == "minimal" and (min_per_group is not None and min_per_group < 3):
            return "EXPLORATORY"
        can_high = tier in ("moderate", "large", "ample", "high")
        if can_high and mmc_core_pct is not None and sss_overlap is not None and mmc_core_pct >= 0.60 and sss_overlap >= 0.40 and warn_count == 0:
            return "HIGH"
        if tier == "small" or (mmc_core_pct is not None and 0.40 <= mmc_core_pct < 0.60) or (sss_overlap is not None and 0.25 <= sss_overlap < 0.40) or warn_count in (1, 2):
            return "MODERATE"
        if tier == "minimal" or (mmc_core_pct is not None and mmc_core_pct < 0.40) or (sss_overlap is not None and sss_overlap < 0.25) or warn_count >= 3:
            return "LOW"
        return "MODERATE"

    def pill_link(filename, label):
        fpath = os.path.join(output_dir, filename)
        if os.path.exists(fpath):
            rel_path = f"{results_folder_name}/{filename}"
            return f'<a class="pill" href="{rel_path}" target="_blank">{label}</a>'
        return f'<span class="pill muted">{label}</span>'

    # Calculated values for Summary
    n_con = safe_int(stats.get('n_con'))
    n_test = safe_int(stats.get('n_test'))
    total_samples = n_con + n_test
    
    minfi_up = safe_int(stats.get('minfi_up'))
    minfi_down = safe_int(stats.get('minfi_down'))
    minfi_total = minfi_up + minfi_down
    
    sesame_up = safe_int(stats.get('sesame_up'))
    sesame_down = safe_int(stats.get('sesame_down'))
    sesame_total = sesame_up + sesame_down

    sesame_native_up = safe_int(stats.get('sesame_native_up'))
    sesame_native_down = safe_int(stats.get('sesame_native_down'))
    sesame_native_total = sesame_native_up + sesame_native_down
    
    intersect_up = safe_int(stats.get('intersect_up'))
    intersect_down = safe_int(stats.get('intersect_down'))
    intersect_total = intersect_up + intersect_down

    intersect_native_up = safe_int(stats.get('intersect_native_up'))
    intersect_native_down = safe_int(stats.get('intersect_native_down'))
    intersect_native_total = intersect_native_up + intersect_native_down
    analysis_params = load_analysis_params()
    qc_summary = load_qc_summary()
    cell_summary = load_cell_deconv_summary()
    cell_assoc = load_cell_assoc()
    warnings = collect_dashboard_warnings(stats, analysis_params, qc_summary, cell_summary, cell_assoc)

    primary_branch = stats.get("primary_branch") or ""
    primary_metrics = load_metrics(primary_branch) if primary_branch else {}
    lam_primary = safe_float(primary_metrics.get("lambda")) if primary_metrics else None
    if lam_primary is not None:
        if lam_primary > 1.2:
            warnings.append(f"Lambda={lam_primary:.2f} suggests p-value inflation; prioritize stronger effect sizes.")
        elif lam_primary < 0.9:
            warnings.append(f"Lambda={lam_primary:.2f} suggests potential over-correction; review batch settings.")

    lg_status = (stats.get("primary_lambda_guard_status") or "").lower()
    if lg_status in ("failed", "warn"):
        warnings.append(f"Lambda guard status: {lg_status}. Consider stricter thresholds or manual review.")

    tier_row = load_crf_sample_tier()
    min_per_group = safe_int(tier_row.get("min_per_group")) if tier_row else None
    tier_label = (analysis_params.get("crf_sample_tier") or stats.get("crf_sample_tier") or (tier_row.get("tier") if tier_row else "") or "").lower()

    mmc_rows = load_crf_summary_file("CRF_MMC_Summary.csv", primary_branch)
    ncs_rows = load_crf_summary_file("CRF_NCS_Summary.csv", primary_branch)
    sss_rows = load_crf_summary_file("CRF_SSS_Summary.csv", primary_branch)
    mmc_summary = summarize_mmc(mmc_rows)
    ncs_summary = summarize_ncs(ncs_rows)
    sss_summary = summarize_sss(sss_rows)

    warnings_by_level = classify_warnings(warnings)
    warn_critical = warnings_by_level["critical"]
    warn_important = warnings_by_level["important"]
    warn_info = warnings_by_level["info"]
    warn_count = len(warn_critical) + len(warn_important)
    verdict = compute_verdict(tier_label, min_per_group, mmc_summary.get("core_pct"), sss_summary.get("overlap"), warn_count, len(warn_critical))

    verdict_class_map = {
        "HIGH": "verdict-high",
        "MODERATE": "verdict-moderate",
        "LOW": "verdict-low",
        "EXPLORATORY": "verdict-exploratory",
    }
    verdict_class = verdict_class_map.get(verdict, "verdict-moderate")

    consensus_path = ""
    for fname in ("Intersection_Native_Consensus_DMPs.csv", "Intersection_Consensus_DMPs.csv"):
        path = os.path.join(output_dir, fname)
        if os.path.exists(path):
            consensus_path = path
            break
    consensus_summary = summarize_consensus(consensus_path) if consensus_path else {}
    abs_deltas = consensus_summary.get("abs_deltas") or []
    median_abs = None
    max_abs = consensus_summary.get("max_abs")
    if abs_deltas:
        try:
            median_abs = statistics.median(abs_deltas)
        except Exception:
            median_abs = None
    effect_tag = effect_label(median_abs)

    top_genes = []
    gene_counts = consensus_summary.get("gene_counts") or {}
    if gene_counts:
        sorted_genes = sorted(gene_counts.items(), key=lambda kv: kv[1], reverse=True)
        for gene, count in sorted_genes[:3]:
            direction = ""
            if consensus_summary.get("gene_sign") and gene in consensus_summary["gene_sign"]:
                direction = "hyper" if consensus_summary["gene_sign"][gene] >= 0 else "hypo"
            if direction:
                top_genes.append(f"{gene} ({count} CpGs, {direction})")
            else:
                top_genes.append(f"{gene} ({count} CpGs)")

    key_findings = []
    if intersect_native_total:
        key_findings.append(
            f"{intersect_native_total} high-confidence DMPs (Intersection Native: ↑{intersect_native_up} ↓{intersect_native_down})."
        )
    elif intersect_total:
        key_findings.append(
            f"{intersect_total} high-confidence DMPs (Intersection Strict: ↑{intersect_up} ↓{intersect_down})."
        )
    else:
        key_findings.append("No consensus DMPs passed the significance thresholds.")

    if top_genes:
        key_findings.append("Top genes: " + ", ".join(top_genes) + ".")
    else:
        key_findings.append("Top genes: not available (no gene annotations in consensus table).")

    if median_abs is not None:
        effect_bits = f"Median |Δβ|={median_abs:.3f}"
        if max_abs is not None:
            effect_bits += f", max |Δβ|={max_abs:.3f}"
        if effect_tag:
            effect_bits += f" ({effect_tag} effect sizes)"
        key_findings.append(effect_bits + ".")
    else:
        key_findings.append("Effect size summary: not available.")

    mmc_core_pct = mmc_summary.get("core_pct")
    mmc_pct_display = f"{mmc_core_pct * 100:.0f}%" if mmc_core_pct is not None else "N/A"
    mmc_spearman = mmc_summary.get("spearman")
    mmc_spearman_display = f"{mmc_spearman:.2f}" if mmc_spearman is not None else "N/A"
    mmc_bar = min(max(mmc_core_pct * 100, 0), 100) if mmc_core_pct is not None else 0

    sss_overlap = sss_summary.get("overlap")
    sss_overlap_display = f"{sss_overlap:.2f}" if sss_overlap is not None else "N/A"
    sss_sign = sss_summary.get("sign")
    sss_sign_display = f"{sss_sign * 100:.0f}%" if sss_sign is not None else "N/A"
    sss_bar = min(max(sss_overlap * 100, 0), 100) if sss_overlap is not None else 0

    ncs_lambda = ncs_summary.get("lambda") if ncs_summary else None
    if ncs_lambda is None:
        ncs_lambda = lam_primary
    ncs_ci_low = ncs_summary.get("lambda_ci_low") if ncs_summary else None
    ncs_ci_high = ncs_summary.get("lambda_ci_high") if ncs_summary else None
    ncs_label = "N/A"
    if ncs_lambda is not None:
        ncs_label = f"λ={ncs_lambda:.2f}"
        if ncs_ci_low is not None and ncs_ci_high is not None:
            ncs_label += f" [{ncs_ci_low:.2f}, {ncs_ci_high:.2f}]"

    def _progress_class(value):
        if value is None:
            return ""
        if value >= 60:
            return ""
        if value >= 40:
            return "warn"
        return "danger"

    mmc_progress_class = _progress_class(mmc_core_pct * 100 if mmc_core_pct is not None else None)
    sss_progress_class = _progress_class(sss_overlap * 100 if sss_overlap is not None else None)

    warn_detail_blocks = []
    if warn_critical:
        warn_detail_blocks.append("<strong>Critical</strong><ul>" + "".join([f"<li>{html.escape(w)}</li>" for w in warn_critical]) + "</ul>")
    if warn_important:
        warn_detail_blocks.append("<strong>Important</strong><ul>" + "".join([f"<li>{html.escape(w)}</li>" for w in warn_important]) + "</ul>")
    if warn_info:
        warn_detail_blocks.append("<strong>Info</strong><ul>" + "".join([f"<li>{html.escape(w)}</li>" for w in warn_info]) + "</ul>")
    warn_details_html = "".join(warn_detail_blocks) if warn_detail_blocks else "<p>No warnings detected.</p>"

    # CSS Style Block (Using format to avoid curly brace hell)
    style_block = """
        @import url('https://fonts.googleapis.com/css2?family=Manrope:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
        :root {
            --ink: #1d2628;
            --primary: #2f6b64;
            --accent: #2f7a7b;
            --accent-2: #f2b45d;
            --rose: #e07a5f;
            --bg: #f6f1e8;
            --card: #ffffff;
            --line: #e4ddd2;
            --muted: #6f7a76;
            --shadow: 0 12px 30px rgba(27, 38, 40, 0.12);
            --warning-bg: #fff3d6;
            --warning-border: #f0c36a;
            --warning-text: #7a4c14;
        }
        html { scroll-behavior: smooth; }
        body {
            font-family: "Manrope", "IBM Plex Sans", "Helvetica Neue", sans-serif;
            margin: 0;
            color: var(--ink);
            background-color: var(--bg);
            background-image:
                radial-gradient(circle at 10% 10%, rgba(242, 180, 93, 0.12), transparent 45%),
                radial-gradient(circle at 90% 20%, rgba(47, 122, 123, 0.12), transparent 40%),
                radial-gradient(circle at 50% 85%, rgba(224, 122, 95, 0.08), transparent 45%);
        }
        h1, h2, h3 { font-family: "Space Grotesk", "Manrope", sans-serif; margin: 0; }
        header {
            background: linear-gradient(120deg, #1f3f3c, #2f6b64 55%, #3d8b7d);
            color: #f7f5f0;
            padding: 3.2rem 1.5rem 3.5rem;
        }
        .hero {
            max-width: 1200px;
            margin: 0 auto;
            display: grid;
            grid-template-columns: minmax(260px, 1.2fr) minmax(240px, 0.8fr);
            gap: 28px;
            align-items: center;
        }
        .hero-eyebrow {
            text-transform: uppercase;
            font-size: 0.75rem;
            letter-spacing: 0.2em;
            color: rgba(247, 245, 240, 0.65);
        }
        .hero-title { font-size: clamp(2rem, 2.8vw, 2.8rem); margin-top: 0.5rem; }
        .hero-sub { margin-top: 0.6rem; font-size: 1rem; opacity: 0.85; }
        .hero-tags { display: flex; flex-wrap: wrap; gap: 0.5rem; margin-top: 1.2rem; }
        .tag {
            background: rgba(255, 255, 255, 0.12);
            padding: 0.35rem 0.75rem;
            border-radius: 999px;
            font-size: 0.85rem;
        }
        .hero-grid { display: grid; gap: 12px; }
        .hero-card {
            background: rgba(255, 255, 255, 0.12);
            border: 1px solid rgba(255, 255, 255, 0.2);
            border-radius: 14px;
            padding: 1rem 1.2rem;
            backdrop-filter: blur(8px);
        }
        .hero-card-title { font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.12em; opacity: 0.7; }
        .hero-card-value { font-size: 1.8rem; font-weight: 700; margin-top: 0.3rem; }
        .hero-card-sub { font-size: 0.9rem; opacity: 0.8; margin-top: 0.4rem; }

        .container { max-width: 1200px; margin: -2.2rem auto 3rem; padding: 0 20px 2rem; }
        .jump-bar {
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            padding: 0.75rem 1rem;
            border-radius: 14px;
            background: var(--card);
            box-shadow: var(--shadow);
            position: sticky;
            top: 0;
            z-index: 100;
        }
        .jump-bar a {
            text-decoration: none;
            color: var(--ink);
            font-weight: 600;
            font-size: 0.9rem;
            padding: 0.35rem 0.75rem;
            border-radius: 999px;
            background: #f2eee6;
        }
        .jump-bar a:hover { background: #e7dfd2; }

        .section-title {
            font-size: 1.3rem;
            font-weight: 700;
            color: var(--primary);
            margin: 2.2rem 0 0.8rem;
        }
        .tab-content .section-title {
            font-size: 1.05rem;
            color: var(--ink);
            margin-top: 1.6rem;
        }
        .section-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 18px; margin-bottom: 20px; }
        .callout {
            background: #fdf9f2;
            border: 1px solid var(--line);
            border-radius: 14px;
            padding: 1rem 1.2rem;
            color: var(--muted);
            box-shadow: 0 8px 22px rgba(27, 38, 40, 0.06);
        }
        .callout strong { color: var(--ink); }
        .warning-panel {
            background: var(--warning-bg);
            border: 1px solid var(--warning-border);
            border-radius: 16px;
            padding: 1rem 1.2rem;
            color: var(--warning-text);
            box-shadow: var(--shadow);
            margin: 1.6rem 0 2rem;
        }
        .warning-panel h3 {
            margin: 0 0 0.6rem;
            font-size: 1rem;
        }
        .warning-panel ul {
            margin: 0;
            padding-left: 1.2rem;
        }

        .summary-card {
            background: var(--card);
            border: 1px solid var(--line);
            border-radius: 22px;
            padding: 1.6rem;
            box-shadow: var(--shadow);
            margin-bottom: 1.6rem;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: minmax(240px, 0.9fr) minmax(260px, 1.1fr);
            gap: 18px;
            align-items: start;
        }
        .verdict-badge {
            display: inline-flex;
            align-items: center;
            gap: 8px;
            padding: 0.45rem 0.9rem;
            border-radius: 12px;
            font-weight: 700;
            letter-spacing: 0.08em;
            font-size: 0.85rem;
            text-transform: uppercase;
        }
        .verdict-high { background: #10B981; color: #fff; }
        .verdict-moderate { background: #F59E0B; color: #fff; }
        .verdict-low { background: #F97316; color: #fff; }
        .verdict-exploratory { background: #EF4444; color: #fff; }
        .summary-meta { margin-top: 0.6rem; color: var(--muted); font-size: 0.95rem; }
        .summary-kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 10px; margin-top: 1rem; }
        .kpi-card { background: #f7f2ea; border-radius: 12px; padding: 0.7rem 0.9rem; }
        .kpi-label { font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.12em; color: var(--muted); }
        .kpi-value { font-weight: 700; margin-top: 0.2rem; }
        .summary-list { margin: 0.4rem 0 0; padding-left: 1.2rem; color: var(--ink); }
        .summary-list li { margin-bottom: 0.4rem; }
        .summary-stats { margin-top: 1.2rem; display: grid; gap: 12px; }
        .stat-row { display: grid; grid-template-columns: minmax(140px, 0.6fr) minmax(140px, 1.4fr); gap: 10px; align-items: center; }
        .stat-label { font-weight: 600; color: var(--muted); font-size: 0.9rem; }
        .stat-value { font-weight: 700; }
        .progress-bar { width: 100%; height: 8px; border-radius: 999px; background: #e6dfd4; overflow: hidden; }
        .progress-bar span { display: block; height: 100%; background: var(--accent); border-radius: 999px; }
        .progress-bar.warn span { background: var(--accent-2); }
        .progress-bar.danger span { background: var(--rose); }

        .warn-summary { margin-top: 1.2rem; }
        .warn-chips { display: flex; flex-wrap: wrap; gap: 8px; margin: 0.6rem 0; }
        .warn-chip { padding: 0.35rem 0.7rem; border-radius: 999px; font-size: 0.85rem; font-weight: 700; }
        .warn-critical { background: #fde2e2; color: #b91c1c; }
        .warn-important { background: #fef3c7; color: #92400e; }
        .warn-info { background: #dbeafe; color: #1d4ed8; }
        .warn-details { background: #fdf9f2; border: 1px dashed var(--line); border-radius: 12px; padding: 0.8rem 1rem; }
        .warn-details summary { cursor: pointer; font-weight: 700; }

        .step-card {
            background: var(--card);
            border: 1px solid var(--line);
            border-radius: 16px;
            padding: 1.2rem 1.4rem;
            box-shadow: var(--shadow);
        }
        .step-card h3 { margin-bottom: 0.4rem; }
        .step-label {
            font-size: 0.75rem;
            letter-spacing: 0.18em;
            text-transform: uppercase;
            color: var(--accent);
            margin-bottom: 0.4rem;
        }
        .pill-row { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 0.8rem; }
        .pill {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 0.35rem 0.7rem;
            border-radius: 999px;
            background: #eef4f2;
            color: var(--primary);
            text-decoration: none;
            font-size: 0.85rem;
            font-weight: 600;
        }
        .pill.muted { background: #f1f0ed; color: #9aa4a1; }

        .nav-tabs { display: flex; gap: 10px; flex-wrap: wrap; margin-top: 1.5rem; }
        .nav-tab {
            padding: 0.6rem 1.4rem;
            cursor: pointer;
            font-weight: 700;
            color: var(--primary);
            border-radius: 999px;
            background: #eef4f2;
            transition: all 0.2s ease;
        }
        .nav-tab.active { background: var(--accent); color: #fff; }
        .nav-tab:hover { transform: translateY(-2px); }

        .tab-shell { margin-top: 1.2rem; }
        .tab-content { display: none; animation: riseIn 0.5s ease; }
        .tab-content.active { display: block; }

        .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 18px; }
        .card {
            background: var(--card);
            border-radius: 16px;
            border: 1px solid var(--line);
            box-shadow: 0 10px 22px rgba(27, 38, 40, 0.08);
            overflow: hidden;
            display: flex;
            flex-direction: column;
            animation: riseIn 0.55s ease both;
            animation-delay: calc(var(--i, 1) * 40ms);
        }
        .card-header {
            padding: 0.9rem 1rem;
            font-weight: 700;
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: #f7f2ea;
        }
        .card-badge {
            font-size: 0.7rem;
            padding: 2px 8px;
            border-radius: 999px;
            font-weight: 700;
            letter-spacing: 0.08em;
        }
        .badge-plot { background: rgba(47, 122, 123, 0.15); color: var(--accent); }
        .badge-table { background: rgba(242, 180, 93, 0.2); color: #a86814; }
        .badge-csv, .badge-tsv { background: rgba(224, 122, 95, 0.2); color: #b05643; }
        .badge-doc { background: rgba(31, 63, 60, 0.15); color: #1f3f3c; }
        .card-body { padding: 1.2rem 1.2rem 1.4rem; flex-grow: 1; }
        .card-desc { font-size: 0.9rem; color: var(--muted); line-height: 1.5; margin-bottom: 1.2rem; }
        .btn {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 100%;
            padding: 0.65rem 0;
            background: var(--accent);
            color: #fff;
            text-decoration: none;
            border-radius: 10px;
            font-weight: 700;
            transition: all 0.2s ease;
        }
        .btn:hover { background: #256a6b; transform: translateY(-1px); }

        .empty-msg { text-align: center; padding: 2rem; color: var(--muted); font-style: italic; }
        .metrics-card { background: var(--card); border: 1px solid var(--line); border-radius: 14px; padding: 1rem 1.2rem; box-shadow: 0 8px 18px rgba(27, 38, 40, 0.06); }
        .metrics-title { font-weight: 700; margin-bottom: 0.6rem; color: var(--primary); }
        .metrics-row { display: flex; justify-content: space-between; gap: 12px; font-size: 0.95rem; padding: 6px 0; border-bottom: 1px dashed #eee4d8; }
        .metrics-row:last-child { border-bottom: none; }

        @keyframes riseIn { from { opacity: 0; transform: translateY(12px); } to { opacity: 1; transform: translateY(0); } }
        @media (max-width: 900px) {
            header { padding: 2.6rem 1.4rem 3rem; }
            .hero { grid-template-columns: 1fr; }
            .container { margin-top: -1.6rem; }
            .summary-grid { grid-template-columns: 1fr; }
            .stat-row { grid-template-columns: 1fr; }
        }
    """

    # Construct HTML using lists and joins to be cleaner
    html_parts = []
    
    # Header
    html_parts.append(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IlluMeta Analysis: {group_test} vs {group_con}</title>
    <style>{style_block}</style>
</head>
<body>

<header>
    <div class="hero">
        <div>
            <div class="hero-eyebrow">IlluMeta Analysis</div>
            <h1 class="hero-title">Analysis Dashboard</h1>
            <p class="hero-sub">Comparison: <strong>{group_test}</strong> (Test) vs <strong>{group_con}</strong> (Control)</p>
            <div class="hero-tags">
                <span class="tag">Samples: {total_samples}</span>
                <span class="tag">Intersection (Native): {intersect_native_total} · Strict: {intersect_total}</span>
                <span class="tag">Minfi: {minfi_total} · Sesame strict: {sesame_total} · Sesame native: {sesame_native_total}</span>
            </div>
        </div>
        <div class="hero-grid">
            <div class="hero-card">
                <div class="hero-card-title">Samples</div>
                <div class="hero-card-value">{total_samples}</div>
                <div class="hero-card-sub">Control {n_con} / Test {n_test}</div>
            </div>
            <div class="hero-card">
                <div class="hero-card-title">Intersection (Native · High-confidence)</div>
                <div class="hero-card-value">{intersect_native_total}</div>
                <div class="hero-card-sub">Native ▲ {intersect_native_up} ▼ {intersect_native_down} · Strict ▲ {intersect_up} ▼ {intersect_down}</div>
            </div>
            <div class="hero-card">
                <div class="hero-card-title">Pipeline DMPs</div>
                <div class="hero-card-value">{minfi_total} | {sesame_total} | {sesame_native_total}</div>
                <div class="hero-card-sub">Minfi ▲ {minfi_up} ▼ {minfi_down} · Strict ▲ {sesame_up} ▼ {sesame_down} · Native ▲ {sesame_native_up} ▼ {sesame_native_down}</div>
            </div>
        </div>
    </div>
</header>

<div class="container">
    <div class="jump-bar">
        <a href="#summary">Executive Summary</a>
        <a href="#start">Start Here</a>
        <a href="#controls">Run Controls & QC</a>
        <a href="#docs">Run Documentation</a>
        <a href="#pipelines">Results</a>
    </div>
""")

    # -- "What should I look at first?" callout --
    html_parts.append("""
    <div class="callout" style="background:#eaf6ed; border-color:#b6d7c0; margin-bottom:1.2rem;">
        <strong>New to IlluMeta? Start here:</strong>
        (1) Check the <b>verdict badge</b> below for an overall confidence rating (hover for details).
        (2) Scroll to <a href="#start" style="color:#1a6b3c;font-weight:600;">Beginner Path</a> for a 3-step guide.
        (3) Open the <b>Intersection (Native)</b> tab for your primary, high-confidence results.
        <br><small style="color:#555;">DMPs = differentially methylated positions (individual CpG sites). DMRs = differentially methylated regions (clusters of CpGs). FDR = false discovery rate (adjusted p-value controlling for multiple testing).</small>
    </div>
""")

    key_findings_html = "".join([f"<li>{html.escape(item)}</li>" for item in key_findings])
    html_parts.append(f"""
    <div id="summary" class="section-title">Executive Summary</div>
    <div class="summary-card">
        <div class="summary-grid">
            <div>
                <div class="verdict-badge {verdict_class}" title="HIGH = Results are robust and suitable for publication. MODERATE = Usable results; review warnings carefully. LOW = Interpret with caution; significant limitations detected. EXPLORATORY = Hypothesis-generating only; not suitable for definitive conclusions.">{verdict} Confidence</div>
                <div style="font-size:0.78rem;color:#555;margin:0.25rem 0 0.4rem;">{"Results are robust and suitable for publication." if verdict == "HIGH" else "Usable results; review warnings and consider validation." if verdict == "MODERATE" else "Interpret with caution; significant limitations detected." if verdict == "LOW" else "Hypothesis-generating only; not for definitive conclusions."}</div>
                <div class="summary-meta">Samples: {total_samples} (Control {n_con} / Test {n_test}) · <span title="CRF tier reflects your study's statistical power based on sample size: Minimal (&lt;12), Small (12-23), Moderate (24-49), Large (&ge;50). Higher tiers enable more robust statistical checks.">Tier: {tier_label.upper() if tier_label else "N/A"}</span></div>
                <div class="summary-kpis">
                    <div class="kpi-card" title="The normalization pipeline used as the primary analysis branch. Results from this pipeline drive the executive summary metrics.">
                        <div class="kpi-label">Primary Branch</div>
                        <div class="kpi-value">{primary_branch or "N/A"}</div>
                    </div>
                    <div class="kpi-card" title="Number of consensus DMPs (differentially methylated positions) found significant in BOTH pipelines with concordant direction. These are your most reliable results.">
                        <div class="kpi-label">Intersection (Native)</div>
                        <div class="kpi-value">{intersect_native_total}</div>
                    </div>
                    <div class="kpi-card" title="DMP counts for each pipeline: Minfi (Noob normalization) / Sesame (strict, Minfi-aligned) / Sesame Native (pOOBAH-preserved). Higher counts in individual pipelines with low intersection may indicate pipeline-specific artifacts.">
                        <div class="kpi-label">Pipelines</div>
                        <div class="kpi-value">{minfi_total} / {sesame_total} / {sesame_native_total}</div>
                    </div>
                </div>
            </div>
            <div>
                <div class="metrics-title">Key Findings</div>
                <ul class="summary-list">{key_findings_html}</ul>
            </div>
        </div>
        <div class="summary-stats">
            <div class="stat-row" title="Method concordance measures how well Minfi and SeSAMe agree. Core % = percentage of DMPs found by both methods. Spearman rho = rank correlation of effect sizes. Higher values indicate more robust results.">
                <div class="stat-label">Method concordance</div>
                <div>
                    <div class="stat-value">{mmc_pct_display} core · ρ={mmc_spearman_display}</div>
                    <div class="progress-bar {mmc_progress_class}"><span style="width:{mmc_bar:.0f}%"></span></div>
                </div>
            </div>
            <div class="stat-row" title="Signal Stability Score (SSS) assesses reproducibility via leave-one-out analysis. Overlap = how stable the top hits are when removing one sample. Sign = how consistent the effect direction is. Values near 1.0 are ideal.">
                <div class="stat-label">Stability (SSS)</div>
                <div>
                    <div class="stat-value">Overlap {sss_overlap_display} · Sign {sss_sign_display}</div>
                    <div class="progress-bar {sss_progress_class}"><span style="width:{sss_bar:.0f}%"></span></div>
                </div>
            </div>
            <div class="stat-row" title="Negative Control Score (NCS) uses permutation testing to check if your results contain more signal than expected by chance. 'Passed' means real signal exceeds the null distribution.">
                <div class="stat-label">Negative control</div>
                <div class="stat-value">{ncs_label}</div>
            </div>
        </div>
        <div class="warn-summary">
            <div class="metrics-title">Warnings & Recommendations</div>
            <div class="warn-chips">
                <span class="warn-chip warn-critical">Critical {len(warn_critical)}</span>
                <span class="warn-chip warn-important">Important {len(warn_important)}</span>
                <span class="warn-chip warn-info">Info {len(warn_info)}</span>
            </div>
            <details class="warn-details">
                <summary>View details</summary>
                {warn_details_html}
            </details>
        </div>
    </div>
""")

    guide_steps = [
        ("Step 1", "Check sample quality",
         "Before interpreting any results, confirm that samples passed QC. Look for: >90% samples passing, "
         "detection P failure fraction < 0.20, and consistent signal intensity across samples. If many samples "
         "fail, your downstream results may be unreliable.",
         [("QC_Summary.csv", "QC Summary"), ("Sample_QC_DetectionP_FailFraction.html", "Detection P"), ("Sample_QC_Intensity_Medians.html", "Intensity")]),
        ("Step 2", "Review the High-confidence Intersection",
         "The Intersection (Native) tab shows CpGs significant in BOTH pipelines (Minfi and SeSAMe) with the "
         "same direction of change. These are your most reliable hits. Check the concordance plot to confirm "
         "pipelines agree, and use the consensus DMP table for your primary results.",
         [("Intersection_Native_Consensus_DMPs.html", "Intersection DMPs (Native)"), ("Intersection_Consensus_DMPs.html", "Intersection DMPs (Strict)"),
          ("Intersection_Native_LogFC_Concordance.html", "Concordance (Native)"), ("Intersection_LogFC_Concordance.html", "Concordance (Strict)")]),
        ("Step 3", "Explore individual pipelines",
         "For deeper analysis, check the volcano plots (effect size vs. significance), Manhattan plots (genomic "
         "distribution), and DMR tables (region-level changes). Compare Minfi and Sesame results to understand "
         "where methods agree and diverge.",
         [("Minfi_Volcano.html", "Minfi Volcano"), ("Sesame_Volcano.html", "Sesame Volcano (Strict)"),
          ("Sesame_Native_Volcano.html", "Sesame Volcano (Native)"), ("Minfi_DMRs_Table.html", "Minfi DMRs")]),
    ]
    html_parts.append('<div id="start" class="section-title">Beginner Path</div>')
    html_parts.append('<div class="section-grid">')
    for idx, (label, title, desc, links) in enumerate(guide_steps, start=1):
        link_html = " ".join([pill_link(fname, lbl) for fname, lbl in links])
        html_parts.append(f'''<div class="step-card" style="--i:{idx};">
            <div class="step-label">{label}</div>
            <h3>{title}</h3>
            <p class="card-desc">{desc}</p>
            <div class="pill-row">{link_html}</div>
        </div>''')
    html_parts.append('</div>')
    html_parts.append('<div class="callout"><strong>Beginner tip:</strong> If results look inconsistent, check QC and batch diagnostics before focusing on DMPs/DMRs. '
                      'Lambda (&lambda;) near 1.0 is ideal; values &gt;1.5 suggest bias. '
                      'Zero consensus DMPs with many pipeline-specific DMPs means the two methods disagree&mdash;investigate QC and batch correction.</div>')

    if analysis_params or qc_summary:
        html_parts.append('<div id="controls" class="section-title">Run Controls & QC</div>')
        html_parts.append('<div class="section-grid">')
        if analysis_params:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">Analysis Parameters</div>\n')
            html_parts.append(f'            <div class="metrics-row" title="FDR (False Discovery Rate): adjusted p-value threshold controlling for multiple testing. |logFC|: minimum log2 fold-change in methylation M-values. |DeltaBeta|: minimum absolute difference in beta-values (0-1 scale). Stricter thresholds yield fewer but more confident results."><span>FDR / |logFC| / |DeltaBeta|</span><span>{analysis_params.get("pval_threshold", "N/A")} / {analysis_params.get("lfc_threshold", "N/A")} / {analysis_params.get("delta_beta_threshold", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Tissue</span><span>{analysis_params.get("tissue", "N/A")} ({analysis_params.get("tissue_source", "N/A")})</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Array type</span><span>{analysis_params.get("array_type", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row" title="Correction Robustness Framework tier. Minimal (<12): exploratory only, SVA disabled. Small (12-23): limited power. Moderate (24-49): full pipeline. Large (>=50): all checks powered."><span>CRF sample tier</span><span>{analysis_params.get("crf_sample_tier", "N/A")}</span></div>\n')
            cell_ref = analysis_params.get("cell_reference", "")
            cell_ref_platform = analysis_params.get("cell_reference_platform", "")
            cell_ref_label = "default" if not cell_ref else f"{cell_ref} ({cell_ref_platform or 'auto'})"
            html_parts.append(f'            <div class="metrics-row"><span>Cell reference</span><span>{cell_ref_label}</span></div>\n')
            auto_cov_enabled = analysis_params.get("auto_covariates_enabled", "N/A")
            auto_cov_exclude = analysis_params.get("auto_covariates_exclude_group_associated", "N/A")
            auto_cov_p = analysis_params.get("auto_covariates_group_assoc_p_threshold", "N/A")
            auto_cov_cor = analysis_params.get("auto_covariates_max_cor", "N/A")
            html_parts.append(f'            <div class="metrics-row"><span>Auto covariates</span><span>enabled={auto_cov_enabled}, alpha={analysis_params.get("auto_covariate_alpha", "N/A")}, max_pcs={analysis_params.get("auto_covariate_max_pcs", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Auto covariate guards</span><span>exclude_group={auto_cov_exclude}, p&lt;{auto_cov_p}, max_cor={auto_cov_cor}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Sesame native NA max</span><span>{analysis_params.get("sesame_native_na_max_frac", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Over/Under guard</span><span>{analysis_params.get("overcorrection_guard_ratio", "N/A")} / {analysis_params.get("undercorrection_guard_min_sv", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>SVA</span><span>{analysis_params.get("sva_enabled", "N/A")} ({analysis_params.get("sva_inclusion_rule", "N/A")})</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Clock covariates</span><span>{analysis_params.get("clock_covariates_enabled", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>SV group filter</span><span>P<{analysis_params.get("sv_group_p_threshold", "N/A")}, Eta^2>{analysis_params.get("sv_group_eta2_threshold", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>PVCA min n</span><span>{analysis_params.get("pvca_min_samples", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>DMR (maxgap/p)</span><span>{analysis_params.get("dmr_maxgap", "N/A")} / {analysis_params.get("dmr_p_cutoff", "N/A")}</span></div>\n')
            html_parts.append('        </div>\n')
        if qc_summary:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">QC Summary</div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Samples (pass/fail)</span><span>{qc_summary.get("Samples_passed_QC", "N/A")} / {qc_summary.get("Samples_failed_QC", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Probes final</span><span>{qc_summary.get("Probes_final", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Probes dropped (det/cross/SNP/sex)</span><span>{qc_summary.get("Probes_failed_detection", "N/A")} / {qc_summary.get("Probes_cross_reactive", "N/A")} / {qc_summary.get("Probes_with_SNPs", "N/A")} / {qc_summary.get("Probes_sex_chromosomes", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Sex mismatches (flagged/dropped)</span><span>{qc_summary.get("Sex_mismatch_samples", "N/A")} / {qc_summary.get("Samples_failed_sex_mismatch", "N/A")}</span></div>\n')
            html_parts.append('        </div>\n')
        html_parts.append('</div>')

    run_docs = [
        ("methods.md", "Methods Summary", "Auto-generated methods text for the run.", "DOC"),
        ("analysis_parameters.json", "Analysis Parameters", "Run settings and thresholds (JSON).", "DOC"),
        ("sessionInfo.txt", "R Session Info", "Full R session and package versions.", "DOC"),
        ("code_version.txt", "Code Version", "Git commit hash when available.", "DOC"),
        ("decision_ledger.tsv", "Decision Ledger", "Automated decision log with reasons (TSV).", "DOC"),
        ("Correction_Adequacy_Report.txt", "Correction Adequacy Report", "CAF summary for the primary branch.", "DOC"),
        ("Correction_Adequacy_Summary.csv", "Correction Adequacy Summary", "CAF metrics for the primary branch.", "CSV"),
        ("Correction_Robustness_Report.txt", "Correction Robustness Report", "CRF sample-size adaptive report.", "DOC"),
        ("CRF_MMC_Summary.csv", "CRF MMC Summary (CSV)", "Multi-method concordance summary for the primary branch.", "CSV"),
        ("CRF_NCS_Summary.csv", "CRF NCS Summary (CSV)", "Negative control stability summary for the primary branch.", "CSV"),
        ("CRF_SSS_Summary.csv", "CRF SSS Summary (CSV)", "Split-sample stability summary for the primary branch.", "CSV"),
        ("CRF_Sample_Tier.csv", "CRF Sample Tier (CSV)", "Sample size tiering and warnings.", "CSV"),
        ("Preflight_Summary.csv", "Preflight Summary (CSV)", "Preflight checks and group counts.", "CSV"),
        ("Preflight_IDAT_Pairs.csv", "Preflight IDAT Pairs (CSV)", "IDAT pair existence per sample.", "CSV"),
        ("QC_Summary.csv", "QC Summary (CSV)", "Sample/probe QC counts.", "CSV"),
        ("Probe_Filter_Summary.csv", "Probe Filter Summary (CSV)", "Probe filtering counts by step.", "CSV"),
        ("CrossReactive_Removed_Probes.csv", "Cross-reactive Probes (CSV)", "Removed cross-reactive probe IDs.", "CSV"),
        ("Sample_QC_Metrics.csv", "Sample QC Metrics (CSV)", "Per-sample QC metrics.", "CSV"),
        ("Sex_Check_Summary.csv", "Sex Check Summary (CSV)", "Predicted vs metadata sex (minfi::getSex).", "CSV"),
        ("Sex_Check_Plot.html", "Sex Check Plot (HTML)", "X vs Y median intensity scatter.", "HTML"),
        ("Sex_Check_Plot.png", "Sex Check Plot (PNG)", "X vs Y median intensity scatter.", "PNG"),
        ("Sex_Check.csv", "Sex Check (Legacy CSV)", "Legacy sex check output.", "CSV"),
        ("Tier3_Eligibility.csv", "Tier3 Eligibility (CSV)", "Tier3 confounding eligibility counts.", "CSV"),
        ("Cell_Deconvolution_Summary.csv", "Cell Deconvolution Summary (CSV)", "Cell deconvolution methods and K.", "CSV"),
        ("cell_counts_merged.csv", "Cell Composition (CSV)", "Merged cell composition covariates (if available).", "CSV"),
        ("Cell_Group_Association.csv", "Cell vs Group (CSV)", "Cell composition vs group association.", "CSV"),
        ("Cell_vs_Group_Association.csv", "Cell vs Group (Alt CSV)", "Cell composition vs group association.", "CSV"),
    ]
    extra_cell_files = []
    try:
        for fname in os.listdir(output_dir):
            if fname.startswith("cell_counts_") and fname.endswith(".csv") and fname not in {"cell_counts_merged.csv"}:
                extra_cell_files.append(fname)
    except Exception:
        extra_cell_files = []
    for fname in sorted(extra_cell_files):
        run_docs.append((fname, f"{fname}", "Cell composition reference output (CSV).", "CSV"))
    doc_cards = []
    for fname, title, desc, badge in run_docs:
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            rel_path = f"{results_folder_name}/{fname}"
            doc_cards.append((title, desc, rel_path, badge))
    if doc_cards:
        html_parts.append('<div id="docs" class="section-title">Run Documentation</div>')
        html_parts.append('<div class="grid">')
        card_index = 0
        for title, desc, rel_path, badge in doc_cards:
            card_index += 1
            badge_class = f"badge-{badge.lower()}"
            card_html = f"""            <div class="card" style="--i:{card_index};">
                <div class="card-header">
                    {title}
                    <span class="card-badge {badge_class}">{badge}</span>
                </div>
                <div class="card-body">
                    <p class="card-desc">{desc}</p>
                    <a href="{rel_path}" target="_blank" class="btn">Open File</a>
                </div>
            </div>
"""
            html_parts.append(card_html)
        html_parts.append('</div>')

    html_parts.append('<div id="pipelines" class="section-title">Results by Pipeline</div>')
    html_parts.append('<div class="nav-tabs" id="navTabs">')
    for idx, (pipe_id, pipe_label, _) in enumerate(pipeline_defs):
        active_class = " active" if idx == 0 else ""
        html_parts.append(
            f'    <div class="nav-tab{active_class}" data-tab="{pipe_id}" onclick="switchTab(\'{pipe_id}\')">{pipe_label}</div>'
        )
    html_parts.append('</div>')
    html_parts.append('<div class="tab-shell">')

    # Loop for Pipelines
    for idx, (pipe_id, _, _) in enumerate(pipeline_defs):
        active_class = "active" if idx == 0 else ""
        html_parts.append(f'    <div id="{pipe_id}" class="tab-content {active_class}">\n')
        
        # Metrics block (if available)
        metrics = load_metrics(pipe_id)
        if metrics:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">Model & Batch Summary</div>\n')
            html_parts.append(f'            <div class="metrics-row" title="Genomic inflation factor. Ideal: 0.9-1.1. Values >1.2 suggest p-value inflation (possible batch effects or widespread signal). Values <0.9 suggest over-correction."><span>λ (inflation)</span><span>{metrics.get("lambda", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Batch method</span><span>{metrics.get("batch_method_applied", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Samples / CpGs</span><span>{metrics.get("n_samples", "N/A")} / {metrics.get("n_cpgs", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Covariates used (n)</span><span>{metrics.get("n_covariates_used", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Covariates used</span><span>{metrics.get("covariates_used", "None") or "None"}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>SVs used</span><span>{metrics.get("n_sv_used", "0")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Batch p<0.05 (Before→After)</span><span>{metrics.get("batch_sig_p_lt_0.05_before", "N/A")} → {metrics.get("batch_sig_p_lt_0.05_after", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Min batch p (Before→After)</span><span>{metrics.get("batch_min_p_before", "N/A")} → {metrics.get("batch_min_p_after", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Group min p (Before→After)</span><span>{metrics.get("group_min_p_before", "N/A")} → {metrics.get("group_min_p_after", "N/A")}</span></div>\n')
            guidance = build_correction_guidance(metrics)
            if guidance:
                html_parts.append(f'            <div class="metrics-row"><span>Correction guide</span><span>{guidance}</span></div>\n')
            if metrics.get("perm_mean_sig") is not None:
                html_parts.append(f'            <div class="metrics-row"><span>Perm mean/max sig (null)</span><span>{metrics.get("perm_mean_sig", "N/A")} / {metrics.get("perm_max_sig", "N/A")}</span></div>\n')
            if metrics.get("vp_primary_group_mean") is not None:
                html_parts.append(f'            <div class="metrics-row"><span>VarPart primary_group</span><span>{metrics.get("vp_primary_group_mean", "N/A")}</span></div>\n')
            if metrics.get("dropped_covariates"):
                html_parts.append(f'            <div class="metrics-row"><span>Dropped covariates</span><span>{metrics.get("dropped_covariates")}</span></div>\n')
            drop_reasons = load_drop_reasons(pipe_id)
            if drop_reasons:
                reasons_txt = ", ".join([f"{k}:{v}" for k, v in drop_reasons.items()])
                html_parts.append(f'            <div class="metrics-row"><span>Drop reasons</span><span>{reasons_txt}</span></div>\n')
            html_parts.append('        </div>\n')
        files_found = 0
        card_index = 0
        sections = pipeline_sections.get(pipe_id, [])

        for section_title, section_items in sections:
            section_cards = []
            for suffix, title, desc, badge in section_items:
                filename = f"{pipe_id}{suffix}"
                file_path = os.path.join(output_dir, filename)
                rel_path = f"{results_folder_name}/{filename}"
                if os.path.exists(file_path):
                    card_index += 1
                    badge_class = f"badge-{badge.lower()}"
                    card_html = f"""            <div class="card" style="--i:{card_index};">
                <div class="card-header">
                    {title}
                    <span class="card-badge {badge_class}">{badge}</span>
                </div>
                <div class="card-body">
                    <p class="card-desc">{desc}</p>
                    <a href="{rel_path}" target="_blank" class="btn">Open Result</a>
                </div>
            </div>
"""
                    section_cards.append(card_html)

            if section_cards:
                files_found += len(section_cards)
                html_parts.append(f'        <div class="section-title">{section_title}</div>\n')
                html_parts.append('        <div class="grid">\n')
                html_parts.extend(section_cards)
                html_parts.append('        </div>\n')

        if files_found == 0:
            html_parts.append('        <div class="empty-msg">No results found for this pipeline.</div>')
             
        html_parts.append('    </div>\n')

    html_parts.append('</div>\n')

    # Footer and Script
    html_parts.append("""
</div>

<script>
    function switchTab(tabName) {
        document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
        document.querySelectorAll('.nav-tab').forEach(el => el.classList.remove('active'));
        document.getElementById(tabName).classList.add('active');
        const match = Array.from(document.querySelectorAll('.nav-tab')).find(el => el.dataset.tab === tabName);
        if (match) match.classList.add('active');
    }
</script>

</body>
</html>
""")

    full_html = "".join(html_parts)

    with open(dashboard_path, "w") as f:
        f.write(full_html)
    
    log(f"[*] Dashboard generated: {dashboard_path}")

def main():
    parser = argparse.ArgumentParser(description="IlluMeta: Automated Illumina DNA Methylation Analysis Pipeline")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Download Command
    parser_download = subparsers.add_parser("download", help="Download GEO data and prepare configuration")
    parser_download.add_argument("gse_id", nargs="?", type=str, help="GEO Series ID (e.g., GSE12345)")
    parser_download.add_argument("-o", "--out-dir", type=str, default=None, help="Output directory (default: ./<GSE_ID>)")
    parser_download.add_argument("--platform", type=str, default="",
                                 help="Optional platform ID (e.g., GPL21145) to select when multiple platforms exist")

    # Search Command
    parser_search = subparsers.add_parser("search", help="Search GEO for Illumina IDAT-backed methylation datasets")
    parser_search.add_argument("-k", "--keywords", required=True,
                               help="Required keywords to AND with IDAT/platform filter (e.g., 'breast cancer')")
    parser_search.add_argument("-o", "--output", default="geo_idat_methylation.tsv", help="Output TSV path")
    parser_search.add_argument("--email", default=None, help="Your email (NCBI etiquette; optional)")
    parser_search.add_argument("--retmax", type=int, default=500, help="Maximum records to fetch (default: 500)")
    parser_search.add_argument("--check-suppl", dest="check_suppl", action="store_true", default=True,
                               help="Check GEO supplementary listings for .idat/.idat.gz/RAW.tar (default: on)")
    parser_search.add_argument("--no-check-suppl", dest="check_suppl", action="store_false",
                               help="Disable supplementary listing check to speed up")
    parser_search.add_argument("--sleep", type=float, default=0.5,
                               help="Sleep seconds between E-utility requests (default: 0.5)")
    
    # Analysis Command
    parser_analysis = subparsers.add_parser("analysis", help="Run analysis based on configuration")
    parser_analysis.add_argument("-i", "--input-dir", type=str, help="Project input directory (containing configure.tsv)")
    parser_analysis.add_argument("-c", "--config", type=str, help="Specific path to configure.tsv (overrides --input-dir)")
    parser_analysis.add_argument("-o", "--output", type=str, help="Output directory for results (default: [input-dir]/[test]_vs_[con]_results)")
    parser_analysis.add_argument("--idat-dir", type=str, help="Path to IDAT directory (default: [config dir]/idat)")
    parser_analysis.add_argument("--group_con", type=str, required=True, help="Control group label (case-insensitive)")
    parser_analysis.add_argument("--group_test", type=str, required=True, help="Test group label (case-insensitive)")
    parser_analysis.add_argument("--auto-group", action="store_true",
                                 help="Auto-fill primary_group from metadata (use --group-column/--group-key or auto-detect)")
    parser_analysis.add_argument("--group-column", type=str,
                                 help="Column to use for auto-group (copied into primary_group)")
    parser_analysis.add_argument("--group-key", type=str,
                                 help="Characteristics key to use for auto-group (e.g., disease, condition)")
    parser_analysis.add_argument("--group-map", type=str,
                                 help="Value mapping for auto-group (e.g., 'normal=Control,tumor=Case')")
    parser_analysis.add_argument("--auto-group-output", type=str,
                                 help="Output path for auto-grouped configure.tsv (default: configure_autogroup.tsv)")
    parser_analysis.add_argument("--auto-group-overwrite", action="store_true",
                                 help="Overwrite existing primary_group values when auto-grouping")
    parser_analysis.add_argument("--auto-group-allow-technical", action="store_true",
                                 help="Allow technical/batch columns if no biological candidate exists (not recommended)")
    parser_analysis.add_argument("--fail-on-missing-idat", action="store_true",
                                 help="Fail if any samples have missing IDAT pairs instead of auto-filtering them (default: auto-filter)")
    parser_analysis.add_argument("--keep-missing-idat", action="store_true",
                                 help="[Deprecated: use --fail-on-missing-idat] Alias for --fail-on-missing-idat")
    parser_analysis.add_argument("--max_plots", type=int, default=10000, help="Max points for interactive plots (default: 10000)")
    parser_analysis.add_argument("--pval", type=float, default=0.05, help="Adjusted P-value threshold (default: 0.05)")
    parser_analysis.add_argument("--lfc", type=float, default=0.5, help="Log2 Fold Change threshold (default: 0.5)")
    parser_analysis.add_argument("--delta-beta", type=float, default=0.0, help="Absolute Delta Beta threshold (default: 0; disabled)")
    parser_analysis.add_argument("--tmp-dir", type=str, help="Custom temporary directory (e.g., external drive)")
    parser_analysis.add_argument("--disable-auto-covariates", action="store_true", help="Disable automatic covariate selection via PCs")
    parser_analysis.add_argument("--disable-sva", action="store_true", help="Disable surrogate variable analysis")
    parser_analysis.add_argument("--include-covariates", type=str, help="Comma-separated covariate names to always try to include (if present in configure.tsv)")
    parser_analysis.add_argument("--include-clock-covariates", action="store_true", help="Compute epigenetic clocks and include them as candidate covariates")
    parser_analysis.add_argument("--beginner-safe", action="store_true",
                                 help="Enable beginner-safe mode (stricter group checks and conservative thresholds)")
    parser_analysis.add_argument("--tissue", type=str, default="Auto", help="Tissue type for cell deconvolution (default Auto = reference-free; options: Auto, Blood, CordBlood, DLPFC, Placenta)")
    parser_analysis.add_argument("--cell-reference", type=str, help="Custom cell reference (package name, package::object, or .rds/.rda path)")
    parser_analysis.add_argument("--cell-reference-platform", type=str, help="Reference platform for cell reference (e.g. IlluminaHumanMethylationEPIC or IlluminaHumanMethylation450k)")
    parser_analysis.add_argument("--positive_controls", type=str, help="Comma-separated list of known marker genes to verify (e.g. 'AHRR,CYP1A1')")
    parser_analysis.add_argument("--marker-list", type=str,
                                 help="Optional CpG marker list (TSV/CSV with CpG column or one CpG per line) for signal preservation checks")
    parser_analysis.add_argument("--cross-reactive-list", type=str,
                                 help="Optional cross-reactive probe list (TSV/CSV with CpG column or one CpG per line)")
    parser_analysis.add_argument("--unsafe-skip-cross-reactive", action="store_true",
                                 help="UNSAFE: skip mandatory cross-reactive probe filtering")
    parser_analysis.add_argument("--disable-cross-reactive", action="store_true",
                                 help="DEPRECATED: use --unsafe-skip-cross-reactive")
    parser_analysis.add_argument("--sex-mismatch-action", type=str, choices=["stop", "drop", "ignore"],
                                 help="Sex mismatch check action (stop|drop|ignore)")
    parser_analysis.add_argument("--sex-check-action", type=str, choices=["stop", "drop", "ignore"],
                                 help="DEPRECATED: use --sex-mismatch-action")
    parser_analysis.add_argument("--sex-check-column", type=str,
                                 help="Metadata column to use for sex mismatch check")
    parser_analysis.add_argument("--disable-sex-check", action="store_true",
                                 help="DEPRECATED: use --sex-mismatch-action ignore (unsafe)")
    parser_analysis.add_argument("--batch-column", type=str,
                                 help="Override batch column name for correction")
    parser_analysis.add_argument("--batch-method", type=str, choices=["none", "combat", "limma", "sva"],
                                 help="Override batch correction method (none|combat|limma|sva)")
    parser_analysis.add_argument("--tier3-min-total-n", type=int,
                                 help="Tier3 minimum total N required to run stratified/meta (default: 20)")
    parser_analysis.add_argument("--tier3-min-per-group-per-stratum", type=int,
                                 help="Tier3 minimum samples per group per stratum (default: 5)")
    parser_analysis.add_argument("--tier3-on-fail", type=str, choices=["stop", "skip"],
                                 help="Tier3 action when ineligible (stop|skip; default: stop)")
    parser_analysis.add_argument("--auto-covariates-enabled", type=str,
                                 help="Enable/disable auto covariates (true/false)")
    parser_analysis.add_argument("--auto-covariates-exclude-group-associated", type=str,
                                 help="Exclude covariates strongly associated with group (true/false)")
    parser_analysis.add_argument("--auto-covariates-group-assoc-p-threshold", type=float,
                                 help="P-value threshold for group association exclusion")
    parser_analysis.add_argument("--auto-covariates-max-cor", type=float,
                                 help="Max allowed absolute correlation among auto covariates")
    parser_analysis.add_argument("--cell-adjustment-on-high-eta2", type=str, choices=["warn", "stop"],
                                 help="Action when cell vs group Eta^2 exceeds threshold (warn|stop)")
    parser_analysis.add_argument("--skip-sesame", action="store_true", help="Skip Sesame pipeline (Minfi only)")
    parser_analysis.add_argument("--sesame-typeinorm", dest="sesame_typeinorm", action="store_true",
                                 help="Enable sesame dyeBiasCorrTypeINorm (default: disabled for stability)")
    parser_analysis.add_argument("--permutations", type=int, default=20,
                                 help="Number of label permutations for null DMP counts (set 0 to use config minimum; set calibration.permutations_min=0 to disable)")
    parser_analysis.add_argument("--config-yaml", type=str, help="Optional YAML config (default: config.yaml next to configure.tsv)")
    parser_analysis.add_argument("--preset", type=str, default="conservative", choices=["conservative", "aggressive"],
                                 help="Optimization preset (conservative or aggressive)")
    parser_analysis.add_argument("--vp-top", type=int, default=5000, help="Top-variable CpGs for variancePartition (default: 5000)")
    parser_analysis.add_argument("--id-column", type=str, help="Column in configure.tsv to treat as sample ID (useful for non-GEO datasets)")
    parser_analysis.add_argument("--min-total-size", type=int, default=6, help="Minimum total sample size required to proceed (default: 6)")
    parser_analysis.add_argument("--qc-intensity-threshold", type=float, default=9.0, help="Median M/U intensity threshold (log2). Set <=0 to disable intensity-based sample drop (default: 9.0)")
    parser_analysis.add_argument("--qc-detection-p-threshold", type=float, default=None,
                                 help="Detection P-value threshold for sample/probe QC (default: 0.05 in config).")
    parser_analysis.add_argument("--qc-sample-fail-frac", type=float, default=None,
                                 help="Sample-level fail fraction (probes with P>threshold) to exclude samples (default: 0.20).")
    parser_analysis.add_argument("--dmr-min-cpgs", type=int, default=None,
                                 help="Minimum CpGs per DMR to retain (default: 2).")
    parser_analysis.add_argument("--dmr-maxgap", type=int, default=None,
                                 help="Maximum gap between CpGs for DMRs (default: 500 or config).")
    parser_analysis.add_argument("--beginner-safe-delta-beta", type=float, default=None,
                                 help="Minimum |DeltaBeta| enforced in beginner-safe mode (default: 0.05 or config).")
    parser_analysis.add_argument("--logit-offset", type=float, default=None,
                                 help="Logit offset for M-value transform (default: 1e-4 or config).")
    parser_analysis.add_argument("--force-idat", action="store_true", help="Force reading IDATs if array sizes differ but types are similar")

    # Doctor Command
    parser_doctor = subparsers.add_parser("doctor", help="Check system and R dependencies (does not install)")
    parser_doctor.add_argument("--skip-pandoc", action="store_true", help="Skip pandoc check")
    
    args = parser.parse_args()

    # Early handling for missing subcommand or required args
    if not args.command:
        parser.print_help()
        return
    if args.command == "download" and not args.gse_id:
        parser_download.print_help()
        log("Example: python illumeta.py download GSE12345 -o /path/to/project")
        return

    if args.command == "doctor":
        run_doctor(args)
        return
    if args.command == "search":
        run_search(args)
        return
    
    check_r_installation()
    if args.command in ("download", "analysis"):
        ensure_r_dependencies()
    if args.command == "analysis":
        check_pandoc_installation()
    
    if args.command == "download":
        run_download(args)
    elif args.command == "analysis":
        run_analysis(args)

if __name__ == "__main__":
    main()
