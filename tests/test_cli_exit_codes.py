"""CLI failure paths must return non-zero status codes."""

from types import SimpleNamespace
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ILLUMETA = os.path.join(BASE_DIR, "illumeta.py")


class CliExitCodeTests(unittest.TestCase):
    def run_illumeta(self, *args):
        return subprocess.run(
            [sys.executable, ILLUMETA, *args],
            capture_output=True,
            text=True,
            timeout=30,
        )

    def make_analysis_args(self, config_path: Path, output_dir: Path, **overrides):
        values = {
            "config": str(config_path),
            "input_dir": None,
            "output": str(output_dir),
            "group_con": "Control",
            "group_test": "Case",
            "beginner_safe": False,
            "min_total_size": 2,
            "delta_beta": 0.0,
            "disable_auto_covariates": False,
            "auto_covariates_enabled": None,
            "skip_sesame": False,
            "sesame_typeinorm": False,
            "beginner_safe_delta_beta": None,
            "auto_group": False,
            "group_column": None,
            "group_key": None,
            "group_map": None,
            "auto_group_output": None,
            "auto_group_overwrite": False,
            "auto_group_allow_technical": False,
            "cross_reactive_list": str(output_dir / "cross_reactive.tsv"),
            "unsafe_skip_cross_reactive": False,
            "config_yaml": None,
            "max_plots": 10000,
            "pval": 0.05,
            "lfc": 0.5,
            "permutations": 0,
            "qc_intensity_threshold": 9.0,
            "marker_list": None,
            "sex_mismatch_action": None,
            "sex_check_column": None,
            "batch_column": None,
            "batch_method": None,
            "tier3_min_total_n": None,
            "tier3_min_per_group_per_stratum": None,
            "tier3_on_fail": None,
            "auto_covariates_exclude_group_associated": None,
            "auto_covariates_group_assoc_p_threshold": None,
            "auto_covariates_max_cor": None,
            "cell_adjustment_on_high_eta2": None,
            "idat_dir": None,
            "fail_on_missing_idat": False,
            "force_idat": False,
            "disable_sva": False,
            "include_covariates": None,
            "include_clock_covariates": False,
            "tissue": "Auto",
            "cell_reference": None,
            "cell_reference_platform": None,
            "positive_controls": None,
            "preset": None,
            "qc_detection_p_threshold": None,
            "qc_sample_fail_frac": None,
            "dmr_min_cpgs": None,
            "dmr_maxgap": None,
            "logit_offset": None,
            "vp_top": None,
            "id_column": None,
            "tmp_dir": None,
        }
        values.update(overrides)
        return SimpleNamespace(**values)

    def preflight_payload(self, config_path: Path):
        return {
            "config_path": str(config_path),
            "idat_dir": None,
            "sample_count_raw": 4,
            "sample_count": 4,
            "group_con": 2,
            "group_test": 2,
            "missing_idat_pairs": 0,
            "missing_idat_preview": [],
            "group_labels": {},
            "batch_candidates": [],
            "warnings": [],
            "column_profile": [],
        }

    def import_illumeta(self):
        sys.path.insert(0, BASE_DIR)
        try:
            import illumeta
        finally:
            if sys.path[0] == BASE_DIR:
                sys.path.pop(0)
        return illumeta

    def test_search_empty_keywords_is_nonzero(self):
        result = self.run_illumeta("search", "--keywords", "", "--no-check-suppl")

        self.assertEqual(result.returncode, 2)
        self.assertIn("--keywords is required", result.stderr)

    def test_download_missing_gse_is_nonzero(self):
        result = self.run_illumeta("download")

        self.assertEqual(result.returncode, 2)
        self.assertIn("Example: python illumeta.py download", result.stdout)

    def test_download_invalid_gse_fails_before_runtime_setup(self):
        result = self.run_illumeta("download", "GSE12345;rm")

        self.assertEqual(result.returncode, 2)
        self.assertIn("GEO Series ID must be GSE followed by digits only", result.stderr)
        self.assertNotIn("Conda env", result.stderr)

    def test_missing_subcommand_is_nonzero(self):
        result = self.run_illumeta()

        self.assertEqual(result.returncode, 2)
        self.assertIn("Available commands", result.stdout)

    def test_download_rscript_launch_failure_is_clean_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "download"
            illumeta = self.import_illumeta()
            args = SimpleNamespace(gse_id="GSE12345", out_dir=str(output_dir), platform="")

            with mock.patch.object(illumeta, "ensure_r_lib_env", side_effect=lambda env: env):
                with mock.patch.object(illumeta, "add_conda_paths", side_effect=lambda env: env):
                    with mock.patch.object(illumeta.subprocess, "run", side_effect=FileNotFoundError("Rscript")):
                        with mock.patch("sys.stderr", new_callable=io.StringIO) as stderr:
                            with self.assertRaises(SystemExit) as ctx:
                                illumeta.run_download(args)

            self.assertEqual(ctx.exception.code, 1)
            self.assertIn("Failed to prepare or launch download step", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())

    def test_download_output_directory_oserror_is_clean_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "download"
            illumeta = self.import_illumeta()
            args = SimpleNamespace(gse_id="GSE12345", out_dir=str(output_dir), platform="")

            with mock.patch.object(illumeta, "ensure_r_lib_env", side_effect=lambda env: env):
                with mock.patch.object(illumeta, "add_conda_paths", side_effect=lambda env: env):
                    with mock.patch.object(illumeta.os, "makedirs", side_effect=PermissionError("readonly")):
                        with mock.patch.object(illumeta.subprocess, "run") as run_mock:
                            with mock.patch("sys.stderr", new_callable=io.StringIO) as stderr:
                                with self.assertRaises(SystemExit) as ctx:
                                    illumeta.run_download(args)

            self.assertEqual(ctx.exception.code, 1)
            run_mock.assert_not_called()
            self.assertIn("Failed to prepare or launch download step", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())

    def test_r_dependency_setup_launch_failure_is_clean_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / ".r_setup_done"
            marker.write_text("stale marker\n", encoding="utf-8")
            illumeta = self.import_illumeta()

            with mock.patch.dict(os.environ, {"ILLUMETA_FORCE_SETUP": "1"}):
                with mock.patch.object(illumeta, "SETUP_MARKER", str(marker)):
                    with mock.patch.object(illumeta, "ensure_r_lib_env", side_effect=lambda env: env):
                        with mock.patch.object(illumeta, "add_conda_paths", side_effect=lambda env: env):
                            with mock.patch.object(illumeta, "detect_r_major_minor", return_value=None):
                                with mock.patch.object(illumeta.subprocess, "run", side_effect=FileNotFoundError("Rscript")):
                                    with mock.patch("sys.stderr", new_callable=io.StringIO) as stderr:
                                        with self.assertRaises(SystemExit) as ctx:
                                            illumeta.ensure_r_dependencies()

            self.assertEqual(ctx.exception.code, 1)
            self.assertFalse(marker.exists())
            self.assertIn("Failed to launch R dependency setup", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())

    def test_missing_config_does_not_delete_unrelated_failure_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work_dir = root / "work"
            config_dir = root / "config"
            output_dir = root / "out"
            work_dir.mkdir()
            config_dir.mkdir()
            old_summary = config_dir / "failure_summary.json"
            old_reason = config_dir / "failure_reason.txt"
            old_summary.write_text('{"code":"OLD"}\n', encoding="utf-8")
            old_reason.write_text("OLD: previous run\n", encoding="utf-8")
            sys.path.insert(0, BASE_DIR)
            import illumeta

            old_cwd = os.getcwd()
            try:
                os.chdir(work_dir)
                with self.assertRaises(SystemExit) as ctx:
                    illumeta.run_analysis(
                        SimpleNamespace(
                            config=str(config_dir / "missing_config.tsv"),
                            input_dir=None,
                            output=str(output_dir),
                            group_con="control",
                            group_test="case",
                        )
                    )
            finally:
                os.chdir(old_cwd)
                if sys.path[0] == BASE_DIR:
                    sys.path.pop(0)

            self.assertEqual(ctx.exception.code, 1)
            self.assertEqual(old_summary.read_text(encoding="utf-8"), '{"code":"OLD"}\n')
            self.assertEqual(old_reason.read_text(encoding="utf-8"), "OLD: previous run\n")

    def test_missing_config_fails_before_r_dependency_setup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work_dir = root / "work"
            work_dir.mkdir()
            missing_config = root / "missing_config.tsv"
            output_dir = root / "out"
            env = os.environ.copy()
            env["ILLUMETA_ALLOW_NON_CONDA"] = "1"
            env["ILLUMETA_FORCE_SETUP"] = "1"

            result = subprocess.run(
                [
                    sys.executable,
                    ILLUMETA,
                    "analysis",
                    "--config",
                    str(missing_config),
                    "--output",
                    str(output_dir),
                    "--group_con",
                    "control",
                    "--group_test",
                    "case",
                ],
                capture_output=True,
                cwd=work_dir,
                text=True,
                timeout=5,
                env=env,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("Configuration file", result.stderr)
            self.assertNotIn("Ensuring R dependencies", result.stdout + result.stderr)

    def test_analysis_timeout_parser_accepts_unlimited_values(self):
        sys.path.insert(0, BASE_DIR)
        try:
            import illumeta

            for value in ("0", "-1", "none", "off", "unlimited"):
                with self.subTest(value=value):
                    self.assertIsNone(illumeta.parse_analysis_timeout(value))
            self.assertEqual(illumeta.parse_analysis_timeout("123"), 123)
            self.assertEqual(illumeta.parse_analysis_timeout("100.5"), 100.5)
        finally:
            if sys.path[0] == BASE_DIR:
                sys.path.pop(0)

    def test_analysis_rscript_launch_failure_writes_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "configure.tsv"
            config_path.write_text("Basename\tprimary_group\n", encoding="utf-8")
            output_dir = root / "out"
            illumeta = self.import_illumeta()

            with mock.patch.object(illumeta, "preflight_analysis", return_value=self.preflight_payload(config_path)):
                with mock.patch.object(illumeta.subprocess, "run", side_effect=FileNotFoundError("Rscript")):
                    with self.assertRaises(SystemExit) as ctx:
                        illumeta.run_analysis(self.make_analysis_args(config_path, output_dir))

            self.assertEqual(ctx.exception.code, 1)
            payload = json.loads((output_dir / "failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["code"], "R_BINARY_MISSING")
            self.assertEqual(payload["stage"], "analysis")

    def test_cross_reactive_prep_failure_writes_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "configure.tsv"
            config_path.write_text("Basename\tprimary_group\n", encoding="utf-8")
            output_dir = root / "out"
            illumeta = self.import_illumeta()
            args = self.make_analysis_args(config_path, output_dir, cross_reactive_list=None)

            with mock.patch.object(illumeta, "ensure_cross_reactive_lists", side_effect=OSError("readonly")):
                with self.assertRaises(SystemExit) as ctx:
                    illumeta.run_analysis(args)

            self.assertEqual(ctx.exception.code, 1)
            payload = json.loads((output_dir / "failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["code"], "CROSS_REACTIVE_PREP_FAILED")
            self.assertEqual(payload["stage"], "cross_reactive")

    def test_output_path_file_writes_failure_marker_to_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "configure.tsv"
            config_path.write_text("Basename\tprimary_group\n", encoding="utf-8")
            output_file = root / "out.html"
            output_file.write_text("not a directory\n", encoding="utf-8")
            illumeta = self.import_illumeta()

            with self.assertRaises(SystemExit) as ctx:
                illumeta.run_analysis(self.make_analysis_args(config_path, output_file))

            self.assertEqual(ctx.exception.code, 1)
            payload = json.loads((root / "failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["code"], "OUTPUT_PATH_IS_FILE")
            self.assertEqual(payload["stage"], "output")

    def test_tmp_dir_file_fails_before_launch_and_writes_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "configure.tsv"
            config_path.write_text("Basename\tprimary_group\n", encoding="utf-8")
            output_dir = root / "out"
            tmp_file = root / "tmpfile"
            tmp_file.write_text("not a directory\n", encoding="utf-8")
            illumeta = self.import_illumeta()

            with mock.patch.object(illumeta, "preflight_analysis", return_value=self.preflight_payload(config_path)) as preflight_mock:
                with mock.patch.object(illumeta, "ensure_r_lib_env", side_effect=lambda env: env):
                    with mock.patch.object(illumeta, "add_conda_paths", side_effect=lambda env: env):
                        with mock.patch.object(illumeta.subprocess, "run") as run_mock:
                            with self.assertRaises(SystemExit) as ctx:
                                illumeta.run_analysis(self.make_analysis_args(config_path, output_dir, tmp_dir=str(tmp_file)))

            self.assertEqual(ctx.exception.code, 1)
            preflight_mock.assert_not_called()
            run_mock.assert_not_called()
            payload = json.loads((output_dir / "failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["code"], "TMP_DIR_NOT_DIRECTORY")
            self.assertEqual(payload["stage"], "tmp_dir")

    def test_dashboard_failure_is_nonfatal_after_successful_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "configure.tsv"
            config_path.write_text("Basename\tprimary_group\n", encoding="utf-8")
            output_dir = root / "out"
            illumeta = self.import_illumeta()

            with mock.patch.object(illumeta, "preflight_analysis", return_value=self.preflight_payload(config_path)):
                with mock.patch.object(illumeta, "ensure_r_lib_env", side_effect=lambda env: env):
                    with mock.patch.object(illumeta.subprocess, "run", return_value=subprocess.CompletedProcess([], 0)):
                        with mock.patch.object(illumeta, "generate_dashboard", side_effect=RuntimeError("bad dashboard")):
                            illumeta.run_analysis(self.make_analysis_args(config_path, output_dir))

            self.assertFalse((output_dir / "failure_summary.json").exists())
            payload = json.loads((output_dir / "dashboard_failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["code"], "DASHBOARD_FAILED")
            self.assertEqual(payload["stage"], "dashboard")


if __name__ == "__main__":
    unittest.main()
