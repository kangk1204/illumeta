"""CLI failure paths must return non-zero status codes."""

from types import SimpleNamespace
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


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

    def test_search_empty_keywords_is_nonzero(self):
        result = self.run_illumeta("search", "--keywords", "", "--no-check-suppl")

        self.assertEqual(result.returncode, 2)
        self.assertIn("--keywords is required", result.stderr)

    def test_download_missing_gse_is_nonzero(self):
        result = self.run_illumeta("download")

        self.assertEqual(result.returncode, 2)
        self.assertIn("Example: python illumeta.py download", result.stdout)

    def test_missing_subcommand_is_nonzero(self):
        result = self.run_illumeta()

        self.assertEqual(result.returncode, 2)
        self.assertIn("Available commands", result.stdout)

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


if __name__ == "__main__":
    unittest.main()
