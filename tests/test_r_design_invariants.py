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

    def test_manuscript_facing_effect_labels_use_m_value_scale(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertNotIn("log2FC", source)
        self.assertIn("limma logFC on M-value scale", source)

    def test_gene_annotation_exports_keep_full_gene_names(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertNotIn("res$Gene <- sapply(res$Gene, clean_gene_names)", source)
        self.assertIn("res$Gene_Display <- sapply(res$Gene, clean_gene_names)", source)

    def test_sesame_reference_cell_counts_are_tissue_gated(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertIn("sesame_reference_cell_counts_supported <- function", source)
        sesame_block = source.split("estimate_sesame_cell_counts(sdf_list", 1)[0]
        self.assertIn("sesame_reference_cell_counts_supported(tissue_use)", sesame_block)

    def test_crf_rss_summary_filename_is_used_by_figure_script(self):
        script = os.path.join(BASE_DIR, "scripts", "generate_application_note_figures.py")
        with open(script, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertIn("CRF_RSS_Summary.csv", source)
        self.assertNotIn("CRF_SSS_Summary.csv", source)

    def test_lambda_guard_simplify_is_labeled_diagnostic(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertNotIn('lambda_guard_status <- "simplified"', source)
        self.assertIn('lambda_guard_status <- "diagnostic_simplified"', source)

    def test_lambda_guard_preserves_existing_group_level_order(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertNotIn("targets_use[[group_col]] <- factor(as.character(targets_use[[group_col]]))", source)
        self.assertIn("levels(targets[[group_col]])", source)
        self.assertIn("targets_use[[group_col]] <- factor(group_values, levels = group_levels)", source)

    def test_consensus_failure_is_fatal_by_default(self):
        """[contract] A consensus computation error must NOT silently collapse to a
        biological 0/0 result: run_intersection must stop() by default and only
        downgrade to a recorded consensus_status='failed' when allow_consensus_failure
        is explicitly enabled. Currently this contract had zero test coverage."""
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        body = extract_r_function(source, "run_intersection")
        # Fail-loud-by-default gate on the explicit opt-in flag.
        self.assertIn('get0("allow_consensus_failure", ifnotfound = FALSE)', body)
        self.assertIn("stop(sprintf(\"Consensus (%s) computation failed", body)
        # The silent-zero path is reachable ONLY inside the allowed branch.
        gate = body.index('get0("allow_consensus_failure"')
        stop_idx = body.index("stop(sprintf(\"Consensus (%s) computation failed", gate)
        zero_idx = body.index('consensus_counts$up <<- 0', gate)
        # The stop() (default path) precedes the silent up<<-0 (allowed-only path).
        self.assertLess(stop_idx, zero_idx)
        self.assertIn('consensus_counts$status <<- "failed"', body)

    def test_genomic_lambda_and_guard_boundary_are_functionally_correct(self):
        """[contract] Functionally exercise the statistic the guard consumes:
        compute_genomic_lambda must be ~1 for uniform p-values and >>1.5 for an
        inflated set, and the guard decision (lambda<=threshold -> ok, else triggered)
        must label those two regimes correctly. Previously only the label string was
        grep-checked, never the numeric threshold behavior."""
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        lambda_fn = extract_r_function(source, "compute_genomic_lambda")
        code = textwrap.dedent(
            """
            {fn}
            thr <- 1.5
            decide <- function(l) if (!is.finite(l)) "missing_lambda" else if (l <= thr) "ok" else "triggered"
            lam_unif <- compute_genomic_lambda(ppoints(5000))          # well-calibrated
            lam_infl <- compute_genomic_lambda(rep(1e-10, 2000))       # grossly inflated
            stopifnot(abs(lam_unif - 1) < 0.1)
            stopifnot(lam_infl > 1.5)
            stopifnot(decide(lam_unif) == "ok")
            stopifnot(decide(lam_infl) == "triggered")
            # Boundary: exactly at threshold is 'ok' (<= is inclusive), just above is 'triggered'.
            stopifnot(decide(1.5) == "ok")
            stopifnot(decide(1.5000001) == "triggered")
            cat("OK\\n")
            """
        ).format(fn=lambda_fn)
        self.run_r(code)
        # Pin the source decision operator so the functional expectation stays bound to code.
        self.assertIn("lambda_val <= lambda_guard_threshold", source)
        self.assertIn('reason = "lambda_above_threshold"', source)

    def test_sex_check_maps_to_genome_before_getsex(self):
        """[contract] getSex has no RGChannelSet method; it must be called on a mapped
        object or the mandatory sex-mismatch QC silently no-ops via the error handler."""
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        self.assertIn("minfi::getSex(minfi::mapToGenome(minfi::preprocessRaw(rgSet)))", source)
        # The raw-rgSet call must not return: getSex(rgSet) directly is the bug.
        self.assertNotIn("minfi::getSex(rgSet)", source)

    def test_consensus_row_selection_is_na_safe(self):
        """[contract] Consensus rows must be selected NA-safely so the emitted table
        matches the na.rm-based up/down counts (no phantom all-NA rows)."""
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        self.assertIn("(is_up %in% TRUE) | (is_down %in% TRUE)", source)
        self.assertNotIn("concord[is_up | is_down, , drop = FALSE]", source)

    def test_consensus_reports_selection_p_values_as_primary(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertNotIn("consensus_df$P.Value <- consensus_df$P.Fisher", source)
        self.assertNotIn("consensus_df$adj.P.Val <- consensus_df$adj.P.Fisher", source)
        self.assertIn("consensus_df$P.Value.fisher_ranking <- consensus_df$P.Fisher", source)
        self.assertIn("consensus_df$P.Value <- consensus_df$P.Value.selection", source)
        self.assertIn("consensus_df$adj.P.Val <- consensus_df$adj.P.Val.selection", source)

    def test_concordance_plot_keeps_all_consensus_points(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertIn("consensus_plot_df <- concord[concord$Consensus %in% TRUE", source)
        self.assertIn("background_df <- concord[!(concord$Consensus %in% TRUE)", source)
        self.assertIn("all consensus points retained", source)

    def test_batch_eval_heatmaps_share_fill_scale(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertIn("make_batch_eval_heatmap <- function", source)
        self.assertIn("batch_fill_limits <- c(0, max(1, max(batch_eval_logp", source)
        self.assertIn("fill_limits = batch_fill_limits", source)

    def test_r_start_removes_stale_success_artifacts(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertIn("stale_success_files <- c(", source)
        self.assertIn('"summary.json"', source)
        self.assertIn('"methods.md"', source)
        self.assertIn("unlink(stale_success_paths, force = TRUE)", source)

    def test_cdx_public_metric_label_is_cdx_score(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        cdx_block = source.split('write.csv(cdx_df, file.path(out_dir, paste0(prefix, "_CDx_Summary.csv"))', 1)[0]
        cdx_block = cdx_block.rsplit("cdx_df <- data.frame", 1)[1]
        self.assertIn('"cdx_score"', cdx_block)
        self.assertNotIn('"cai"', cdx_block)

    def test_calibration_is_unified_two_sided_lambda(self):
        """Batch-method selection AND the CDx report must use ONE calibration
        definition: the shared two-sided null-lambda helper. No KS-ratio score and
        no saturating null-FPR-vs-alpha comparison drive selection/reporting."""
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        # Single shared helper implementing the two-sided log2 lambda closeness.
        self.assertIn("calibration_score_from_lambda <- function(lambda, lambda_tol = 1.5)", source)
        self.assertIn("clamp01(1 - abs(log2(lambda)) / log2(tol))", source)
        # Both the selector and the CDx report call the SAME helper.
        self.assertGreaterEqual(source.count("calibration_score_from_lambda("), 2)
        # The selector no longer derives calibration from the KS-uniformity ratio.
        self.assertNotIn("clamp01(ks_p_med / config_settings$calibration$ks_p)", source)
        # The internal composite field is unified to 'cdx' (no legacy 'cai').
        self.assertIn("cdx = cdx,", source)
        self.assertNotIn("cdx_res$cai", source)

    def test_cdx_robustness_contracts(self):
        """Lock the reliability fixes for public use: legacy caf config honored,
        lambda_tol recorded for reproducibility, stale legacy artifacts cleaned,
        NaN-safe composite, and the selector guard tied to the lambda_tol score."""
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        # Legacy `caf:` config block is honored (not silently ignored) after the rename.
        self.assertIn("config_settings$caf", source)
        self.assertIn("merge_config_defaults(cdx_cfg, config_settings$caf)", source)
        # The new calibration parameter is recorded in the reproducibility outputs.
        self.assertIn("calibration_lambda_tol = config_settings$calibration$lambda_tol", source)
        # Same-dir reruns also clean the pre-rename Correction_Adequacy_* artifacts.
        self.assertIn('"Correction_Adequacy_Report.txt"', source)
        self.assertIn('"Correction_Adequacy_Summary.csv"', source)
        # Composite weighted mean cannot become NaN from a 0-weight finite component.
        self.assertIn("keep <- is.finite(scores) & is.finite(w) & w > 0", source)
        # Selector calibration guard derives from the lambda_tol-driven score.
        self.assertIn("guard_cal <- is.finite(cal_score) && cal_score >= 0.5", source)

    def test_plot_tooltips_use_short_gene_display_labels(self):
        with open(ANALYZE_R, "r", encoding="utf-8") as handle:
            source = handle.read()

        self.assertIn("plot_res$Gene_Label <- if (\"Gene_Display\" %in% colnames(plot_res))", source)
        self.assertIn("<br>Gene:\", Gene_Label", source)
        self.assertIn('gene_map <- if ("Gene_Display" %in% colnames(top_100)) top_100$Gene_Display else top_100$Gene', source)

    def test_download_fails_on_raw_tar_listing_warnings_and_status(self):
        with open(DOWNLOAD_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        safe_extract_fn = extract_r_function(source, "safe_extract_idats_from_tar")

        self.assertIn("list_warnings <- character(0)", safe_extract_fn)
        self.assertIn("Unsafe RAW tar listing failed", safe_extract_fn)
        self.assertIn("bad_status <- is.numeric(untar_status)", safe_extract_fn)
        self.assertIn("Unsafe RAW tar extraction failed", safe_extract_fn)

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

    def test_download_treats_untar_warnings_as_failed_extraction(self):
        with open(DOWNLOAD_R, "r", encoding="utf-8") as handle:
            source = handle.read()
        safe_extract_fn = extract_r_function(source, "safe_extract_idats_from_tar")

        with tempfile.TemporaryDirectory() as tmpdir:
            tar_path = os.path.join(tmpdir, "GSE999_RAW.tar")
            payload_path = os.path.join(tmpdir, "sample_Grn.idat")
            with open(payload_path, "wb") as handle:
                handle.write(b"mock")
            info = tarfile.TarInfo("linked_Grn.idat")
            info.type = tarfile.LNKTYPE
            info.linkname = "/etc/passwd"
            with tarfile.open(tar_path, "w") as tar:
                tar.add(payload_path, arcname="sample_Grn.idat")
                tar.addfile(info)

            self.run_r(
                textwrap.dedent(
                    f"""
                    {safe_extract_fn}
                    stopped <- tryCatch({{
                      safe_extract_idats_from_tar("{tar_path}", exdir = tempfile("extract_"))
                      FALSE
                    }}, error = function(e) grepl("Unsafe RAW tar extraction failed", conditionMessage(e)))
                    if (!isTRUE(stopped)) stop("untar warning/failure was not fatal")
                    """
                )
            )


if __name__ == "__main__":
    unittest.main()
