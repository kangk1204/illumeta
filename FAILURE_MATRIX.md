# Failure-Mode Matrix — `illumeta meta` pipeline (`illumeta_meta.py`)

Tested: 2026-05-30 | Stages: 5 | Tests: 29 | Pass: 29 | Fail: 0 (1 fixed by patch)

Scope: the cross-cohort meta-analysis Python pipeline (feeds manuscript Table 2).
The per-cohort analysis (`r_scripts/analyze.R`) is R and is covered by
`tests/test_r_design_invariants.py`; `illumeta.py` CLI orchestration is a
separate follow-up scope. Tests live in `tests/test_pipeline_robustness.py`
and use only tiny synthetic inputs (no R, no `illumeta.py analysis`).

| Stage \ Failure mode | A2 empty | A4 type | A5 schema | A6 range | A7 encoding | A8 dup | A9 extreme | A10 inf/NaN | B4 header | C2 perm | D1 div0 | F1 contract |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **S1 input resolution** | ✅ | — | ✅ | — | — | — | — | — | ✅ | — | — | ✅ |
| **S2 branch loading**   | ✅ | ⚠️ | ✅ | ✅ | ✅ | ✅ | — | — | — | — | — | — |
| **S3 meta combination** | — | — | — | — | — | — | ✅¹ | ✅ | — | — | ✅ | — |
| **S4 concordance**      | ✅ | — | — | — | — | — | — | — | — | — | — | — |
| **S5 output write**     | ✅ | — | — | — | — | — | — | — | — | ✅ | — | — |

¹ `S3 × A9` was **❌ (crash)** on the first run; fixed by patch and re-verified ✅.

## Cell provenance (test → cell)

- S1·A2 `test_s2_A2_header_only_file_is_empty_with_warning` (manifest all-excluded → empty handled), S1·A5 `test_s1_A5_manifest_missing_result_dir_column`, S1·B4 `test_s1_B4_manifest_no_header`, S1·F1 `test_s1_F1_threshold_out_of_range` (×6) + `test_s1_A3_too_few_cohorts` + `test_s1_F1_pc_r_exceeds_cohorts`
- S2·A5 `test_s2_A5_missing_required_column` + `test_s2_A5_missing_se_and_t`, S2·A8 `test_s2_A8_duplicate_cpg_rejected`, S2·A2 `test_s2_A2_header_only_file_is_empty_with_warning`, S2·A4 `test_s2_A4_nonnumeric_logfc_coerced_to_nan`, S2·A6 `test_s2_A6_out_of_range_pvalue_is_nan` (×4), S2·A7 `test_s2_A7_non_utf8_bytes_do_not_crash`
- S3·A9 `test_s3_A9_tiny_se_underflow_no_crash`, S3·D1 `test_s3_D1_all_invalid_returns_nan_no_crash`, S3·A10 `test_s3_A10_all_nan_pvalues_fisher_pc`, S3·A3 (single cohort) `test_s3_A3_single_valid_cohort`
- S4·A2 `test_s4_A2_zero_candidates_no_crash`
- S5·A2 `test_s5_A2_empty_rows_writes_valid_tsv`, S5·C2 `test_s5_C2_readonly_outdir_raises_typed_error`

## Fixed cell

- **S3 × A9** (`test_s3_A9_tiny_se_underflow_no_crash`): a positive but extreme
  standard error (e.g. `1e-200`) underflows `se*se` to `0.0`, so the `se > 0`
  guard passed and `1.0/(se*se)` raised `ZeroDivisionError` (two sites:
  `_random_effect_meta_one` fixed-effect and random-effect weights). **Genuine
  missing guard.** Patch: introduced `_inv_var_weight(se, ok, extra)` which
  returns weight `0.0` for any non-finite/non-positive denominator (treating the
  cohort as a zero-weight contribution, consistent with existing `se<=0`
  handling). Numerically identical for all real (finite) SE inputs.

## Degraded cell (⚠️ — surfaced, not patched)

- **S2 × A4** (`test_s2_A4_nonnumeric_logfc_coerced_to_nan`): a non-numeric
  `logFC` in one cohort is silently coerced to `NaN` and that cohort is dropped
  for that CpG, with **no per-CpG warning**. This is safe (no crash, no
  corruption) and matches how missing values are handled, but a malformed input
  column would pass unnoticed. Recommendation (optional): count and warn when a
  branch table has a high fraction of unparseable numeric cells.

## Not-applicable rationale

- S1 × A4/A6/A7/A9/A10/D1: S1 resolves paths/manifest only — numeric-content
  pathologies are caught at S2/S3.
- S2 × A9/A10/D1: extreme-value math is exercised at S3 where weights/statistics
  are computed; S2 only parses cells (already `NaN`-coerced).
- S3 × A2/A4/A5/A7/A8/B4/C2: S3 operates on already-loaded in-memory records.
- S4/S5 × A4–A10/D1: receive already-validated rows; only emptiness (A2) and I/O
  (C2) are meaningful.
- C3 (disk full), E* (resource), B1/B2 (file format) not exercised this session.

## Coverage notes

- Categories tested: A (input pathologies), B4 (header), C2 (permission), D1, F1.
- NOT tested (recommend future): B1/B2 (truncated/wrong file types), C3 (disk
  full), E (OOM/timeout), and the `illumeta.py` orchestration pipeline
  (download-URL building, config merge, group auto-detection, results discovery).
- No slow/GPU tests. Run: `python3 -m pytest -q tests/test_pipeline_robustness.py`
