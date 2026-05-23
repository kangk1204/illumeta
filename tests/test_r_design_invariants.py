"""Executable checks for R design-matrix invariants."""

from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import tempfile
import textwrap
import unittest


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYZE_R = os.path.join(BASE_DIR, "r_scripts", "analyze.R")
DOWNLOAD_R = os.path.join(BASE_DIR, "r_scripts", "download.R")


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


class RDesignInvariantTests(unittest.TestCase):
    def run_r(self, code: str) -> None:
        if shutil.which("Rscript") is None:
            self.skipTest("Rscript is not installed")
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
            handle.write(code)
            script_path = handle.name
        try:
            result = subprocess.run(
                ["Rscript", script_path],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        finally:
            os.unlink(script_path)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)

    def test_rank_pruning_preserves_group_contrast_columns(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        drop_fn = extract_r_function(source, "drop_linear_dependencies")

        self.run_r(
            textwrap.dedent(
                f"""
                {drop_fn}
                mat <- cbind(
                  Control = c(1, 1, 0, 0),
                  Case = c(0, 0, 1, 1),
                  ConfoundedControl = c(1, 1, 0, 0),
                  NumericCovariate = c(1, 2, 3, 4)
                )
                res <- drop_linear_dependencies(mat, c("Control", "Case"))
                if (!all(c("Control", "Case") %in% colnames(res$mat))) {{
                  stop("group contrast columns were not preserved")
                }}
                if ("ConfoundedControl" %in% colnames(res$mat)) {{
                  stop("confounded covariate was retained instead of the group column")
                }}
                if (!("ConfoundedControl" %in% res$dropped)) {{
                  stop("confounded covariate was not reported as dropped")
                }}
                """
            )
        )

    def test_covariate_collinearity_tolerates_all_na_correlations(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        map_fn = extract_r_function(source, "map_mm_columns_to_vars")
        collinearity_fn = extract_r_function(source, "compute_covariate_collinearity")

        self.run_r(
            textwrap.dedent(
                f"""
                {map_fn}
                {collinearity_fn}
                targets <- data.frame(
                  Constant = c(1, 1, 1, 1),
                  Alternating = c(0, 1, 0, 1),
                  NumericCovariate = c(1, 2, 3, 4)
                )
                res <- compute_covariate_collinearity(
                  targets,
                  c("Constant", "Alternating", "NumericCovariate")
                )
                if (nrow(res) != 3) stop("unexpected row count")
                constant_row <- res[res$Variable == "Constant", ]
                if (nrow(constant_row) != 1 || !is.na(constant_row$Max_Correlation)) {{
                  stop("all-NA correlations should be represented as NA")
                }}
                """
            )
        )

    def test_tier3_eligibility_allows_two_valid_overlap_strata(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        tier3_fn = extract_r_function(source, "check_tier3_eligibility")

        self.run_r(
            textwrap.dedent(
                f"""
                {tier3_fn}
                targets <- data.frame(
                  batch = c(rep("A", 10), rep("B", 10), rep("C", 6)),
                  primary_group = c(
                    rep(c("Control", "Case"), each = 5),
                    rep(c("Control", "Case"), each = 5),
                    "Control", rep("Case", 5)
                  )
                )
                out_dir <- tempdir()
                res <- check_tier3_eligibility(
                  targets,
                  "batch",
                  "primary_group",
                  min_total_n = 20,
                  min_per_group_per_stratum = 5,
                  out_dir = out_dir
                )
                if (!isTRUE(res$eligible)) stop("two valid overlap strata should be eligible")
                if (res$n_eligible_strata != 2) stop("expected exactly two eligible strata")
                if (res$min_overlap_group_n != 1) stop("expected weak-stratum minimum to remain auditable")
                if (res$min_eligible_stratum_n != 5) stop("eligible-stratum minimum should exclude failed strata")
                if (res$min_stratum_n != 5) stop("legacy min_stratum_n should report eligible-stratum minimum")
                if (!identical(sort(res$failed_strata), "C")) stop("weak stratum C should be reported as failed")
                csv <- read.csv(file.path(out_dir, "Tier3_Eligibility.csv"))
                if (!all(c("Batch_Min_Group_N", "Batch_Eligible", "Min_Overlap_Group_N", "Min_Eligible_Stratum_N", "Eligible_Strata_N", "Failed_Strata") %in% colnames(csv))) {{
                  stop("eligibility CSV missing per-stratum audit columns")
                }}
                c_rows <- csv[csv$Batch == "C", ]
                if (!all(c_rows$Batch_Min_Group_N == 1)) stop("stratum C local minimum was not reported")
                if (any(c_rows$Batch_Eligible)) stop("stratum C should not be batch-eligible")
                """
            )
        )

    def test_tier3_eligibility_ignores_unused_factor_group_levels(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        tier3_fn = extract_r_function(source, "check_tier3_eligibility")

        self.run_r(
            textwrap.dedent(
                f"""
                {tier3_fn}
                targets <- data.frame(
                  batch = c(rep("A", 10), rep("B", 10)),
                  primary_group = factor(
                    c(rep(c("Control", "Case"), each = 5), rep(c("Control", "Case"), each = 5)),
                    levels = c("Control", "Case", "Unused")
                  )
                )
                res <- check_tier3_eligibility(
                  targets,
                  "batch",
                  "primary_group",
                  min_total_n = 20,
                  min_per_group_per_stratum = 5,
                  out_dir = tempdir()
                )
                if (!isTRUE(res$eligible)) stop("unused factor level should not make Tier3 ineligible")
                if (res$min_overlap_group_n != 5) stop("unused factor level leaked into minimum group count")
                if ("Unused" %in% res$counts$Group) stop("unused factor level should not be audited as an observed group")
                """
            )
        )

    def test_batch_overlap_helpers_ignore_unused_factor_group_levels(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        safe_chisq_fn = extract_r_function(source, "safe_chisq_p")
        identify_fn = extract_r_function(source, "identify_overlap_batches")
        confounded_fn = extract_r_function(source, "is_batch_confounded")

        self.run_r(
            textwrap.dedent(
                f"""
                {safe_chisq_fn}
                {identify_fn}
                {confounded_fn}
                meta <- data.frame(
                  batch = c(rep("A", 10), rep("B", 10)),
                  primary_group = factor(
                    c(rep(c("Control", "Case"), each = 5), rep(c("Control", "Case"), each = 5)),
                    levels = c("Control", "Case", "Unused")
                  )
                )
                overlap <- identify_overlap_batches(meta, "batch", "primary_group")
                if (!identical(sort(overlap), c("A", "B"))) stop("unused factor level changed overlap batches")
                if (is_batch_confounded(meta, "batch", "primary_group")) {{
                  stop("unused factor level created a false batch-confounding zero cell")
                }}
                """
            )
        )

    def test_r_summary_exposes_sesame_native_imputation_contract(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        summary_block = source.split("summary_payload <- list(", 1)[1].split("write_json(summary_payload", 1)[0]

        required = [
            "sesame_native_before",
            "sesame_native_retained",
            "sesame_native_dropped_all_na",
            "sesame_native_dropped_na",
            "sesame_native_imputed",
            "sesame_native_impute_method",
        ]
        missing = [field for field in required if field not in summary_block]
        self.assertEqual(missing, [])

    def test_skip_sesame_initializes_absent_branch_summaries(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        start_marker = "sesame_out <- NULL"
        end_marker = "# --- Consensus (Intersection) between Minfi and Sesame ---"
        start = source.find(start_marker)
        end = source.find(end_marker)
        self.assertNotEqual(start, -1, f"Missing marker: {start_marker}")
        self.assertNotEqual(end, -1, f"Missing marker: {end_marker}")
        self.assertLess(start, end)
        sesame_block = source[start:end]
        self.assertIn("sesame_summary <- NULL", sesame_block)
        self.assertLess(sesame_block.index("sesame_summary <- NULL"), sesame_block.index("if (!is.null(beta_sesame_strict))"))
        self.assertIn("sesame_native_summary <- NULL", sesame_block)
        self.assertLess(sesame_block.index("sesame_native_summary <- NULL"), sesame_block.index("if (!is.null(beta_sesame_native))"))

    def test_config_loading_fails_fast_on_invalid_existing_config(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        load_fn = extract_r_function(source, "load_yaml_config")
        preset_fn = extract_r_function(source, "normalize_preset")

        self.run_r(
            textwrap.dedent(
                f"""
                {load_fn}
                {preset_fn}
                absent <- load_yaml_config(file.path(tempdir(), "definitely_missing_config.yaml"))
                if (!is.null(absent)) stop("missing config should use defaults")
                bad_config <- tempfile(fileext = ".yaml")
                writeLines("preset: [", bad_config)
                bad_config_stopped <- tryCatch({{
                  load_yaml_config(bad_config)
                  FALSE
                }}, error = function(e) TRUE)
                if (!bad_config_stopped) stop("existing malformed config should fail fast")
                bad_preset_stopped <- tryCatch({{
                  normalize_preset("fast-and-loose", source = "test")
                  FALSE
                }}, error = function(e) TRUE)
                if (!bad_preset_stopped) stop("invalid preset should fail fast")
                """
            )
        )

    def test_download_rejects_unsafe_raw_tar_member_paths(self):
        with open(DOWNLOAD_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        safe_extract_fn = extract_r_function(source, "safe_extract_idats_from_tar")

        with tempfile.TemporaryDirectory() as tmpdir:
            tar_path = os.path.join(tmpdir, "GSE999_RAW.tar")
            payload_path = os.path.join(tmpdir, "payload.idat")
            with open(payload_path, "wb") as handle:
                handle.write(b"mock")
            with tarfile.open(tar_path, "w") as tar:
                tar.add(payload_path, arcname="../escape_Grn.idat")

            self.run_r(
                textwrap.dedent(
                    f"""
                    {safe_extract_fn}
                    stopped <- tryCatch({{
                      safe_extract_idats_from_tar("{tar_path}", exdir = tempfile("extract_"))
                      FALSE
                    }}, error = function(e) grepl("Unsafe RAW tar member", conditionMessage(e)))
                    if (!isTRUE(stopped)) stop("unsafe tar member was not rejected")
                    """
                )
            )
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "escape_Grn.idat")))

    def test_download_extracts_only_idats_from_safe_raw_tar(self):
        with open(DOWNLOAD_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        safe_extract_fn = extract_r_function(source, "safe_extract_idats_from_tar")

        with tempfile.TemporaryDirectory() as tmpdir:
            tar_path = os.path.join(tmpdir, "GSE999_RAW.tar")
            idat_path = os.path.join(tmpdir, "sample_Grn.idat")
            note_path = os.path.join(tmpdir, "README.txt")
            with open(idat_path, "wb") as handle:
                handle.write(b"mock")
            with open(note_path, "w", encoding="utf-8") as handle:
                handle.write("metadata")
            with tarfile.open(tar_path, "w") as tar:
                tar.add(idat_path, arcname="nested/sample_Grn.idat")
                tar.add(note_path, arcname="README.txt")

            self.run_r(
                textwrap.dedent(
                    f"""
                    {safe_extract_fn}
                    out_dir <- tempfile("extract_")
                    extracted <- safe_extract_idats_from_tar("{tar_path}", exdir = out_dir)
                    if (length(extracted) != 1) stop("expected exactly one extracted IDAT")
                    if (!grepl("sample_Grn.idat$", extracted[[1]])) stop("wrong extracted file")
                    if (file.exists(file.path(out_dir, "README.txt"))) stop("non-IDAT file should not be extracted")
                    """
                )
            )

    def test_download_rejects_idat_symlink_members(self):
        with open(DOWNLOAD_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        safe_extract_fn = extract_r_function(source, "safe_extract_idats_from_tar")

        with tempfile.TemporaryDirectory() as tmpdir:
            tar_path = os.path.join(tmpdir, "GSE999_RAW.tar")
            info = tarfile.TarInfo("linked_Grn.idat")
            info.type = tarfile.SYMTYPE
            info.linkname = "/etc/passwd"
            with tarfile.open(tar_path, "w") as tar:
                tar.addfile(info)

            self.run_r(
                textwrap.dedent(
                    f"""
                    {safe_extract_fn}
                    stopped <- tryCatch({{
                      safe_extract_idats_from_tar("{tar_path}", exdir = tempfile("extract_"))
                      FALSE
                    }}, error = function(e) grepl("Unsafe RAW tar symlink", conditionMessage(e)))
                    if (!isTRUE(stopped)) stop("IDAT symlink member was not rejected")
                    """
                )
            )


if __name__ == "__main__":
    unittest.main()
