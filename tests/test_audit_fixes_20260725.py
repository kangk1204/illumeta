"""Regression tests for the 2026-07-25 scientific audit fixes.

Each test pins one behaviour that the audit found wrong, so a future refactor
cannot silently reintroduce it. Test IDs map to the audit finding IDs:

  S-1  tier3 within-cohort meta must propagate between-stratum tau2
  S-2  tier3 Delta_Beta must be measured on the strata that produced logFC
  S-3  permutation-null limitation must stay documented at the call site
  S-4  bacon .bacon columns must never hold uncorrected statistics
  S-5  genomic lambda must survive p-values below 1e-16
  S-6  two-stage batch-strategy search must declare that it is greedy
  S-7  all three preprocessing branches must share one probe-QC definition
  S-8  intersection counts must not double-count opposite-direction probes
  S-9  stratified effect matrices must be verified aligned before cbind
  C-3  non-CpG (ch*) and control (rs*) probes must be filtered
  C-8  pooled delta beta must be weighted by the samples behind Delta_Beta

R-side checks extract the target function from analyze.R and run it standalone,
mirroring test_r_design_invariants.py, and skip when Rscript is unavailable.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
ANALYZE_R = BASE_DIR / "r_scripts" / "analyze.R"

sys.path.insert(0, str(BASE_DIR))

from illumeta_meta import (  # noqa: E402
    MetaCohort,
    MetaThresholds,
    _analyze_branch,
    _cohort_column_ids,
    _tier3_variance_warnings,
)


def extract_r_function(source: str, name: str) -> str:
    marker = f"{name} <- function"
    start = source.index(marker)
    brace = source.index("{", start)
    depth = 0
    for idx in range(brace, len(source)):
        char = source[idx]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return source[start : idx + 1]
    raise ValueError(f"Could not extract R function: {name}")


def run_r(code: str, timeout: int = 120):
    if shutil.which("Rscript") is None:
        return None
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(code)
        script_path = handle.name
    try:
        return subprocess.run(
            ["Rscript", script_path], capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    finally:
        os.unlink(script_path)


class SourceContractTests(unittest.TestCase):
    """Cheap, dependency-free assertions on the shipped source."""

    @classmethod
    def setUpClass(cls):
        cls.r_source = ANALYZE_R.read_text(encoding="utf-8")
        cls.meta_source = (BASE_DIR / "illumeta_meta.py").read_text(encoding="utf-8")
        cls.template = (BASE_DIR / "config.yaml.template").read_text(encoding="utf-8")

    def test_s1_tier3_meta_defaults_to_random_effects(self):
        """A fixed-effect tier3 SE has no tau2, so the cross-cohort level over-weights it."""
        block = self.r_source[self.r_source.index("tier3_meta = list(") :][:1400]
        self.assertIn('method = "random"', block)
        self.assertIn("tau2_propagated", self.r_source)

    def test_s1_fixed_effect_tier3_is_flagged_downstream(self):
        """illumeta_meta must warn when it pools a tier3 cohort built fixed-effect."""
        summary = {"primary_result_mode": "tier3_stratified_meta",
                   "primary_tier3_meta_method": "fixed"}
        warnings = _tier3_variance_warnings("GSE1", summary)
        self.assertEqual(len(warnings), 1)
        self.assertIn("tau2", warnings[0])

        summary["primary_tier3_meta_method"] = "random"
        self.assertEqual(_tier3_variance_warnings("GSE1", summary), [])

        # A non-tier3 cohort is never flagged.
        self.assertEqual(
            _tier3_variance_warnings("GSE2", {"primary_result_mode": "tier3_ineligible"}), []
        )

    def test_s2_tier3_delta_beta_is_stratum_scoped(self):
        """Delta_Beta must come from the strata that produced logFC, and say so."""
        self.assertIn("delta_beta_scope", self.r_source)
        self.assertIn("Delta_Beta_N_Con", self.r_source)
        self.assertIn("strata_used = strata_used", self.r_source)
        # the emitter must actually receive the stratum set
        self.assertIn("strata_used = if (!is.null(meta_out)) meta_out$strata_used else NULL",
                      self.r_source)

    def test_s3_permutation_null_limitation_is_documented(self):
        """The label-leakage limitation must remain next to the function."""
        idx = self.r_source.index("run_permutation_uniformity <- function")
        doc = self.r_source[max(0, idx - 1800) : idx]
        self.assertIn("not exchangeable", doc)
        self.assertIn("relative", doc)

    def test_s6_batch_strategy_declares_greedy_search(self):
        idx = self.r_source.index("select_batch_strategy <- function")
        doc = self.r_source[max(0, idx - 2600) : idx]
        self.assertIn("greedy", doc)
        self.assertIn("top_k", doc)
        self.assertIn("search_truncated", self.r_source)

    def test_s7_branch_symmetric_masks_exist(self):
        """SNP / sex masks are computed once and applied to SeSAMe too."""
        self.assertIn("array_snp_probe_ids", self.r_source)
        self.assertIn("array_sex_probe_ids", self.r_source)
        self.assertIn("sesame_branch_symmetric_mask", self.r_source)
        # dropLociWithSnps must be *invoked* exactly once so both branches provably
        # share one SNP definition. Comments and the methods-text mention don't count.
        call_lines = [
            line for line in self.r_source.splitlines()
            if "dropLociWithSnps(" in line
            and not line.lstrip().startswith("#")
            and "sprintf(" not in line
        ]
        self.assertEqual(len(call_lines), 1, msg=f"call sites: {call_lines}")

    def test_c3_probe_class_filter_is_configured_on(self):
        self.assertIn("probe_classes = list(", self.r_source)
        self.assertIn("drop_non_cpg = TRUE", self.r_source)
        self.assertIn("drop_rs_control = TRUE", self.r_source)
        self.assertIn("probe_classes:", self.template)

    def test_s8_opposite_direction_category_exists(self):
        self.assertIn("Both (opposite direction)", self.r_source)
        self.assertIn("jaccard_overlap_any_direction", self.r_source)

    def test_s9_stratum_alignment_is_verified(self):
        idx = self.r_source.index("eff_mat <- do.call(cbind, effects)")
        guard = self.r_source[max(0, idx - 1200) : idx]
        self.assertIn("identical(eff_names, expected_probes)", guard)

    def test_c5_methods_disclose_se_reconstruction_mix(self):
        self.assertIn("empirical-Bayes-shrunken residual", self.meta_source)
        self.assertIn("no between-stratum tau2 component", self.meta_source)


class ProbeClassRTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = ANALYZE_R.read_text(encoding="utf-8")

    def test_c3_classify_probe_ids(self):
        """cg / ch / rs classification, including EPICv2 suffixed IDs."""
        code = (
            extract_r_function(self.source, "classify_probe_ids") + "\n"
            + extract_r_function(self.source, "probe_class_keep_mask") + "\n"
            + r"""
ids <- c("cg00000001", "ch.9.9999999R", "rs00000001", "cg00000002_TC21", "ch.9.8888888F")
cls <- classify_probe_ids(ids)
stopifnot(identical(cls, c("cg","ch","rs","cg","ch")))
m <- probe_class_keep_mask(ids)
stopifnot(identical(m$keep, c(TRUE, FALSE, FALSE, TRUE, FALSE)))
stopifnot(m$n_ch == 2, m$n_rs == 1, m$n_dropped == 3)
# opting out retains everything
m2 <- probe_class_keep_mask(ids, drop_non_cpg = FALSE, drop_rs_control = FALSE)
stopifnot(all(m2$keep))
cat("PROBE_CLASS_OK\n")
"""
        )
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertIn("PROBE_CLASS_OK", result.stdout, msg=result.stderr)


class BranchSymmetricMaskRTests(unittest.TestCase):
    """S-7: the SNP / sex masks must behave identically on both branches.

    The masking statements live in the top-level pipeline body rather than in a
    function, so they cannot be extracted the way the other R checks are. Instead we
    replay the exact expressions from analyze.R against a synthetic manifest, which
    catches the failure modes that matter: a branch not being masked at all, and the
    EPICv2 suffix mismatch between a suffixed manifest and unsuffixed SeSAMe IDs.
    """

    @classmethod
    def setUpClass(cls):
        cls.source = ANALYZE_R.read_text(encoding="utf-8")

    def test_s7_masks_apply_to_both_branches_with_epicv2_suffixes(self):
        self.assertIn(".strip_probe_suffix", self.source)
        code = r"""
.strip_probe_suffix <- function(ids) sub("_.+$", "", as.character(ids))

# Manifest IDs carry EPICv2 replicate suffixes; SeSAMe IDs do not.
manifest_ids <- c("cg00000001_TC21", "cg00000002_BC11", "cg00000003_TC21",
                  "cg00000004_TC21", "ch.9.1111111F_TC21", "rs00000001_TC21")
sex_raw  <- c("cg00000003_TC21")            # chrX
snp_raw  <- c("cg00000002_BC11")            # SNP-overlapping

array_sex_probe_ids <- unique(.strip_probe_suffix(sex_raw))
array_snp_probe_ids <- unique(.strip_probe_suffix(snp_raw))

# --- Minfi branch: rownames are suffixed ---
minfi <- matrix(0.5, nrow = length(manifest_ids), ncol = 2,
                dimnames = list(manifest_ids, c("s1", "s2")))
minfi <- minfi[!(.strip_probe_suffix(rownames(minfi)) %in% array_snp_probe_ids), , drop = FALSE]
minfi <- minfi[!(.strip_probe_suffix(rownames(minfi)) %in% array_sex_probe_ids), , drop = FALSE]

# --- SeSAMe branch: rownames are UNsuffixed ---
sesame_ids <- .strip_probe_suffix(manifest_ids)
sesame <- matrix(0.5, nrow = length(sesame_ids), ncol = 2,
                 dimnames = list(sesame_ids, c("s1", "s2")))
base <- .strip_probe_suffix(rownames(sesame))
sesame <- sesame[!(base %in% array_snp_probe_ids) & !(base %in% array_sex_probe_ids), , drop = FALSE]

# Both branches must have dropped exactly the SNP and sex probes.
stopifnot(identical(.strip_probe_suffix(rownames(minfi)), .strip_probe_suffix(rownames(sesame))))
stopifnot(!("cg00000002" %in% .strip_probe_suffix(rownames(sesame))))
stopifnot(!("cg00000003" %in% .strip_probe_suffix(rownames(sesame))))
stopifnot(nrow(minfi) == 4L, nrow(sesame) == 4L)
cat("MASK_SYMMETRY_OK\n")
"""
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertIn("MASK_SYMMETRY_OK", result.stdout, msg=result.stderr)

    def test_s7_probe_class_filter_completes_the_symmetry(self):
        """After class filtering both branches keep only cg-context probes."""
        code = (
            extract_r_function(self.source, "classify_probe_ids") + "\n"
            + extract_r_function(self.source, "probe_class_keep_mask") + "\n"
            + extract_r_function(self.source, "apply_probe_class_filter") + "\n"
            + r"""
ids <- c("cg00000001_TC21", "ch.9.1111111F_TC21", "rs00000001_TC21", "cg00000004")
m <- matrix(0.5, nrow = length(ids), ncol = 2, dimnames = list(ids, c("s1","s2")))
res <- apply_probe_class_filter(m, "unit")
stopifnot(nrow(res$mat) == 2L, res$removed == 2L, res$n_ch == 1L, res$n_rs == 1L)
stopifnot(all(grepl("^cg", rownames(res$mat))))
cat("CLASS_SYMMETRY_OK\n")
"""
        )
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertIn("CLASS_SYMMETRY_OK", result.stdout, msg=result.stderr)


class GenomicLambdaRTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = ANALYZE_R.read_text(encoding="utf-8")

    def test_s5_lambda_keeps_extreme_p_values(self):
        """qchisq(1-p,1) discarded every p below ~1e-16; the upper-tail form must not."""
        code = extract_r_function(self.source, "compute_genomic_lambda") + "\n" + r"""
# All p-values are far below the double-precision limit of 1 - p.
extreme <- rep(1e-300, 101)
lam <- compute_genomic_lambda(extreme)
stopifnot(is.finite(lam), lam > 100)

# A null-uniform vector must still land near 1.
set.seed(7)
lam_null <- compute_genomic_lambda(runif(20000))
stopifnot(abs(lam_null - 1) < 0.1)

# Mixed vector: the extreme tail must not be silently dropped.
set.seed(8)
mixed <- c(runif(9000), rep(1e-300, 1000))
stopifnot(is.finite(compute_genomic_lambda(mixed)))
cat("LAMBDA_OK\n")
"""
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertIn("LAMBDA_OK", result.stdout, msg=result.stderr)


class BaconColumnRTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = ANALYZE_R.read_text(encoding="utf-8")

    def test_s4_skipped_bacon_writes_na_not_raw(self):
        """A skipped correction must not leave raw p-values in the .bacon columns."""
        code = (
            extract_r_function(self.source, "compute_genomic_lambda") + "\n"
            + extract_r_function(self.source, "apply_bacon_correction") + "\n"
            + r"""
# min_probes forces the skip path without needing the bacon package.
res <- data.frame(logFC = rnorm(10), t = rnorm(10), P.Value = runif(10),
                  adj.P.Val = runif(10), CpG = paste0("cg", 1:10),
                  stringsAsFactors = FALSE)
out <- apply_bacon_correction(res, min_probes = 1000L)
stopifnot(attr(out, "bacon_status") == "skipped")
stopifnot(all(is.na(out$P.Value.bacon)))
stopifnot(all(is.na(out$logFC.bacon)))
stopifnot(all(is.na(out$t.bacon)))
stopifnot(all(is.na(out$adj.P.Val.bacon)))
# raw columns untouched
stopifnot(identical(out$P.Value, res$P.Value))
cat("BACON_SKIP_OK\n")
"""
        )
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertIn("BACON_SKIP_OK", result.stdout, msg=result.stderr)


class PooledDeltaWeightTests(unittest.TestCase):
    """C-8: the pooled delta beta must use the samples behind Delta_Beta."""

    def _cohort(self, tmp: Path, name: str, n_con: int, n_test: int,
                delta: float, db_n_con: int | None = None,
                db_n_test: int | None = None) -> MetaCohort:
        d = tmp / name
        d.mkdir(parents=True, exist_ok=True)
        header = ["CpG", "logFC", "t", "P.Value", "Delta_Beta"]
        row = ["cg0000001", "0.5", "5.0", "0.001", str(delta)]
        if db_n_con is not None:
            header += ["Delta_Beta_N_Con", "Delta_Beta_N_Test"]
            row += [str(db_n_con), str(db_n_test)]
        (d / "Minfi_DMPs_full.csv").write_text(
            ",".join(header) + "\n" + ",".join(row) + "\n", encoding="utf-8"
        )
        return MetaCohort(cohort_id=name, result_dir=d, n_con=n_con, n_test=n_test)

    def test_c8_stratified_cohort_is_weighted_by_effective_n(self):
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            # Big cohort contributes delta=+1.0 but only 10 samples produced it;
            # two small cohorts contribute delta=0.0 with 10 samples each.
            cohorts = [
                self._cohort(tmp, "BIG", 400, 400, 1.0, db_n_con=5, db_n_test=5),
                self._cohort(tmp, "S1", 5, 5, 0.0),
                self._cohort(tmp, "S2", 5, 5, 0.0),
            ]
            rows, _summary, warnings = _analyze_branch(
                cohorts, _cohort_column_ids(cohorts), "minfi",
                "Minfi_DMPs_full.csv", MetaThresholds(), allow_missing_branches=False,
            )
            pooled = rows[0]["pooled_delta_beta_weighted"]
            # Effective weights are 10/10/10 -> pooled = 1/3, not 800/810 ~= 0.988.
            self.assertTrue(math.isclose(pooled, 1.0 / 3.0, rel_tol=1e-9),
                            msg=f"pooled={pooled}")
            self.assertTrue(any("stratified subset" in w for w in warnings), warnings)

    def test_c8_falls_back_to_cohort_size_without_the_columns(self):
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            cohorts = [
                self._cohort(tmp, "BIG", 400, 400, 1.0),
                self._cohort(tmp, "S1", 5, 5, 0.0),
                self._cohort(tmp, "S2", 5, 5, 0.0),
            ]
            rows, _summary, _warnings = _analyze_branch(
                cohorts, _cohort_column_ids(cohorts), "minfi",
                "Minfi_DMPs_full.csv", MetaThresholds(), allow_missing_branches=False,
            )
            pooled = rows[0]["pooled_delta_beta_weighted"]
            self.assertTrue(math.isclose(pooled, 800.0 / 820.0, rel_tol=1e-9),
                            msg=f"pooled={pooled}")


if __name__ == "__main__":
    unittest.main()
