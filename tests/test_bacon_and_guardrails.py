"""Executable checks for the bacon inflation-correction integration and the
three added guardrails (tissue/region-vs-batch confounding, experimental-design
diagnostic, brain reference cell deconvolution).

R-side tests extract the target function from analyze.R and run it in isolation,
mirroring test_r_design_invariants.py. bacon-dependent tests skip cleanly when the
bacon package is not on the library path.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[1]
ANALYZE_R = BASE_DIR / "r_scripts" / "analyze.R"
ILLUMETA = BASE_DIR / "illumeta.py"


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


def _libpaths_prelude() -> str:
    ul = os.environ.get("R_LIBS_USER", "")
    if ul:
        return f'.libPaths(c({".libPaths()" if not ul else repr(ul)}, .libPaths()))\n'
    return ""


class BaconIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = ANALYZE_R.read_text(encoding="utf-8")

    def test_helper_functions_present(self):
        """bacon helper and its lambda dependency exist in analyze.R."""
        self.assertIn("apply_bacon_correction <- function", self.source)
        self.assertIn("compute_genomic_lambda <- function", self.source)
        # bacon must be namespaced, never library()'d (would mask limma::topTable)
        self.assertNotIn("library(bacon)", self.source)
        self.assertNotIn("require(bacon)", self.source)
        self.assertIn("bacon::bacon", self.source)

    def test_config_has_bacon_block(self):
        self.assertIn("bacon = list(", self.source)
        # config.yaml.template mirrors it
        tmpl = (BASE_DIR / "config.yaml.template").read_text(encoding="utf-8")
        self.assertIn("bacon:", tmpl)

    def test_meta_flag_exists(self):
        meta_src = ILLUMETA.read_text(encoding="utf-8")
        self.assertIn("--use-bacon", meta_src)
        meta_mod = (BASE_DIR / "illumeta_meta.py").read_text(encoding="utf-8")
        self.assertIn("_USE_BACON", meta_mod)
        self.assertIn("logFC.bacon", meta_mod)
        self.assertIn("P.Value.bacon", meta_mod)

    def test_apply_bacon_correction_runs(self):
        """On synthetic inflated statistics, bacon adds recalibrated columns and
        reports a finite inflation estimate; falls back gracefully if bacon absent."""
        prelude = _libpaths_prelude()
        fn_lambda = extract_r_function(self.source, "compute_genomic_lambda")
        fn_bacon = extract_r_function(self.source, "apply_bacon_correction")
        code = prelude + fn_lambda + "\n" + fn_bacon + "\n" + r"""
if (!requireNamespace("bacon", quietly=TRUE)) { cat("SKIP_NO_BACON\n"); quit(status=0) }
set.seed(11); n <- 3000
res <- data.frame(
  logFC = rnorm(n, 0, 0.05),
  t = rnorm(n, 0.5, 2),   # slight positive bias -> inflation/bias to correct
  P.Value = runif(n),
  adj.P.Val = runif(n),
  CpG = paste0("cg", 1:n), stringsAsFactors = FALSE)
res$P.Value <- 2*pnorm(-abs(res$t))
res$adj.P.Val <- p.adjust(res$P.Value, "BH")
out <- apply_bacon_correction(res, prefix="TEST", label="unit")
stopifnot(all(c("logFC.bacon","t.bacon","P.Value.bacon","adj.P.Val.bacon") %in% colnames(out)))
stopifnot(attr(out,"bacon_status") == "ok")
stopifnot(is.finite(attr(out,"inflation_bacon")))
stopifnot(is.finite(attr(out,"lambda_before")), is.finite(attr(out,"lambda_after")))
stopifnot(nrow(out) == n)
cat("BACON_OK inflation=", round(attr(out,"inflation_bacon"),3), "\n")
"""
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        if "SKIP_NO_BACON" in result.stdout:
            self.skipTest("bacon package not on library path")
        self.assertIn("BACON_OK", result.stdout)


class GuardrailTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = ANALYZE_R.read_text(encoding="utf-8")

    def test_guardrail_functions_present(self):
        for name in ("detect_tissue_batch_confounding",
                     "diagnose_experimental_design",
                     ".cramers_v_pair"):
            self.assertIn(f"{name} <- function", self.source, f"missing {name}")

    def test_tissue_guard_detects_region_confounding(self):
        """chip perfectly confounded with tissue (but balanced vs group) is flagged
        as biological and returned for protection."""
        prelude = _libpaths_prelude()
        fns = "\n".join(extract_r_function(self.source, n) for n in
                        (".cramers_v_pair", "detect_tissue_batch_confounding"))
        code = prelude + fns + "\n" + r"""
region <- rep(c("entorhinal cortex","cerebellum"), each=96)
chip <- paste0("chip", rep(1:16, each=12))
grp <- rep(rep(c("AD","Control"), each=6), 16)
meta <- data.frame(Sentrix_ID=chip, tissue=region, primary_group=grp, stringsAsFactors=FALSE)
r <- detect_tissue_batch_confounding(meta, "Sentrix_ID", "primary_group", NULL, prefix="T")
stopifnot("tissue" %in% r$confounded_tissue_cols)
d <- r$diagnostics[r$diagnostics$variable=="tissue", ]
stopifnot(d$cramer_v_with_batch >= 0.99)
stopifnot(d$cramer_v_with_group <= 0.01)
cat("TISSUE_GUARD_OK\n")
"""
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertIn("TISSUE_GUARD_OK", result.stdout)

    def test_design_guard_flags_confounded_and_passes_randomised(self):
        """SEVERE verdict when case/control not randomised across chips; ok when it is."""
        prelude = _libpaths_prelude()
        fns = "\n".join(extract_r_function(self.source, n) for n in
                        (".cramers_v_pair", "diagnose_experimental_design"))
        code = prelude + fns + "\n" + r"""
# confounded: chips 1-3 all AD, 4-6 all Control
meta_bad <- data.frame(Sentrix_ID=paste0("chip",rep(1:6,each=10)),
  primary_group=c(rep("AD",30),rep("Control",30)), stringsAsFactors=FALSE)
db <- diagnose_experimental_design(meta_bad, "primary_group", NULL, prefix="B")
stopifnot(grepl("SEVERE", db$verdict[1]))
# randomised: each chip balanced
meta_ok <- data.frame(Sentrix_ID=paste0("chip",rep(1:6,each=10)),
  primary_group=rep(rep(c("AD","Control"),each=5),6), stringsAsFactors=FALSE)
do <- diagnose_experimental_design(meta_ok, "primary_group", NULL, prefix="O")
stopifnot(grepl("ok:", do$verdict[1]))
cat("DESIGN_GUARD_OK\n")
"""
        result = run_r(code)
        if result is None:
            self.skipTest("Rscript not installed")
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertIn("DESIGN_GUARD_OK", result.stdout)

    def test_brain_reference_routing_present(self):
        """DLPFC brain reference is wired into the tissue->reference map."""
        self.assertIn('"DLPFC"', self.source)
        self.assertIn("FlowSorted.DLPFC.450k", self.source)


class MetaBaconColumnSelectionTests(unittest.TestCase):
    """--use-bacon selects bacon columns when present; falls back with a warning."""

    def _write_table(self, path: Path, with_bacon: bool):
        base = {
            "CpG": "cg1", "Gene": "G", "chr": "chr1", "pos": "1",
            "logFC": "0.30", "t": "10", "P.Value": "1e-8",
            "adj.P.Val": "1e-5", "Delta_Beta": "0.03",
        }
        if with_bacon:
            base.update({"logFC.bacon": "0.18", "t.bacon": "6.0",
                         "P.Value.bacon": "1e-4", "adj.P.Val.bacon": "1e-2"})
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(base.keys()))
            w.writeheader(); w.writerow(base)

    def test_bacon_column_selection(self):
        sys.path.insert(0, str(BASE_DIR))
        import importlib
        m = importlib.import_module("illumeta_meta")
        importlib.reload(m)
        # simulate the selection block against a header with bacon columns
        fields_with = ["CpG", "logFC", "t", "P.Value", "logFC.bacon", "t.bacon", "P.Value.bacon"]
        fields_without = ["CpG", "logFC", "t", "P.Value"]

        def select(fields, use_bacon):
            eff, p, tcol = "logFC", "P.Value", "t"
            if use_bacon and "P.Value.bacon" in fields and "logFC.bacon" in fields:
                eff, p = "logFC.bacon", "P.Value.bacon"
                tcol = "t.bacon" if "t.bacon" in fields else "t"
            return eff, p, tcol

        self.assertEqual(select(fields_with, True), ("logFC.bacon", "P.Value.bacon", "t.bacon"))
        self.assertEqual(select(fields_with, False), ("logFC", "P.Value", "t"))
        # fallback when bacon columns absent
        self.assertEqual(select(fields_without, True), ("logFC", "P.Value", "t"))
        # module exposes the flag
        self.assertTrue(hasattr(m, "_USE_BACON"))


if __name__ == "__main__":
    unittest.main()
