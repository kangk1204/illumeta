"""Tests for the one-command `illumeta demo` flow.

These verify command construction, cache detection, grouping resolution, and the
offline/download branching WITHOUT running the real pipeline (which needs
minfi/sesame). The demo reuses the validated download+analysis CLIs, so testing
the wiring here is sufficient to guarantee the beginner one-liner dispatches
correctly.
"""

from __future__ import annotations

import json
import os
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import illumeta


def _args(out_dir=None, offline=False):
    return types.SimpleNamespace(out_dir=out_dir, offline=offline)


def _seed_cache(demo_dir: Path):
    """Create a minimal cached demo dataset (configure.tsv + one IDAT)."""
    (demo_dir / "idat").mkdir(parents=True, exist_ok=True)
    (demo_dir / "configure.tsv").write_text("Sample_Name\tprimary_group\n", encoding="utf-8")
    (demo_dir / "idat" / "GSM_x_Grn.idat").write_bytes(b"\x00")


def _clean_demo_env():
    """Patch out any ambient ILLUMETA_DEMO_* env so defaults are deterministic."""
    return mock.patch.dict(
        os.environ,
        {k: "" for k in ("ILLUMETA_DEMO_GSE", "ILLUMETA_DEMO_GROUP_COLUMN",
                         "ILLUMETA_DEMO_GROUP_CON", "ILLUMETA_DEMO_GROUP_TEST")},
        clear=False,
    )


class DemoCommandTests(unittest.TestCase):
    def test_demo_commands_download_then_analysis(self):
        dl, an = illumeta._demo_commands(
            "GSE125605", "/p/demo", need_download=True,
            group_con="Control", group_test="Case", group_column="description",
        )
        self.assertEqual(dl[-4:], ["download", "GSE125605", "-o", "/p/demo"])
        self.assertEqual(
            an[2:],
            ["analysis", "-i", "/p/demo", "--group_con", "Control", "--group_test", "Case",
             "--auto-group", "--group-column", "description", "--tier3-on-fail", "skip"],
        )

    def test_demo_commands_cached_skips_download(self):
        dl, an = illumeta._demo_commands(
            "GSE125605", "/p/demo", need_download=False,
            group_con="Control", group_test="Case", group_column="description",
        )
        self.assertIsNone(dl)
        self.assertIn("analysis", an)

    def test_demo_commands_use_provided_grouping(self):
        _dl, an = illumeta._demo_commands(
            "GSE999", "/p/demo", need_download=True,
            group_con="ctrl", group_test="tumor", group_column="disease state",
        )
        self.assertIn("ctrl", an)
        self.assertIn("tumor", an)
        self.assertIn("disease state", an)

    def test_idats_present_detects_idat_and_gz(self):
        with mock.patch("os.path.isdir", return_value=True), \
             mock.patch("os.listdir", return_value=["a_Grn.idat.gz"]):
            self.assertTrue(illumeta._demo_idats_present("/x"))
        with mock.patch("os.path.isdir", return_value=True), \
             mock.patch("os.listdir", return_value=["readme.txt"]):
            self.assertFalse(illumeta._demo_idats_present("/x"))

    def test_offline_without_cache_exits_nonzero(self):
        import tempfile
        with _clean_demo_env(), mock.patch("illumeta.subprocess.run") as run:
            with tempfile.TemporaryDirectory() as tmp:
                with self.assertRaises(SystemExit) as cm:
                    illumeta.run_demo(_args(out_dir=tmp, offline=True))
                self.assertEqual(cm.exception.code, 1)
            run.assert_not_called()

    def test_cached_offline_runs_analysis_only(self):
        import tempfile
        with _clean_demo_env(), tempfile.TemporaryDirectory() as tmp:
            _seed_cache(Path(tmp))
            calls = []

            def fake_run(cmd, *a, **k):
                calls.append(cmd)
                return types.SimpleNamespace(returncode=0)

            with mock.patch("illumeta.subprocess.run", side_effect=fake_run):
                with self.assertRaises(SystemExit) as cm:
                    illumeta.run_demo(_args(out_dir=tmp, offline=True))
                self.assertEqual(cm.exception.code, 0)
            self.assertEqual(len(calls), 1)
            self.assertIn("analysis", calls[0])
            self.assertNotIn("download", calls[0])

    def test_uncached_downloads_then_analyzes(self):
        import tempfile
        with _clean_demo_env(), tempfile.TemporaryDirectory() as tmp:
            demo_dir = os.path.join(tmp, "projects", "demo")
            calls = []

            def fake_run(cmd, *a, **k):
                calls.append(cmd)
                return types.SimpleNamespace(returncode=0)

            with mock.patch("illumeta.subprocess.run", side_effect=fake_run):
                with self.assertRaises(SystemExit) as cm:
                    illumeta.run_demo(_args(out_dir=demo_dir, offline=False))
                self.assertEqual(cm.exception.code, 0)
            self.assertEqual(len(calls), 2)
            self.assertIn("download", calls[0])
            self.assertIn("analysis", calls[1])

    def test_download_failure_aborts_before_analysis(self):
        import tempfile
        with _clean_demo_env(), tempfile.TemporaryDirectory() as tmp:
            demo_dir = os.path.join(tmp, "projects", "demo")
            calls = []

            def fake_run(cmd, *a, **k):
                calls.append(cmd)
                return types.SimpleNamespace(returncode=2)

            with mock.patch("illumeta.subprocess.run", side_effect=fake_run):
                with self.assertRaises(SystemExit) as cm:
                    illumeta.run_demo(_args(out_dir=demo_dir, offline=False))
                self.assertEqual(cm.exception.code, 2)
            self.assertEqual(len(calls), 1)
            self.assertIn("download", calls[0])

    def test_custom_gse_with_grouping_env_runs(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            demo_dir = os.path.join(tmp, "projects", "demo")
            calls = []

            def fake_run(cmd, *a, **k):
                calls.append(cmd)
                return types.SimpleNamespace(returncode=0)

            env = {
                "ILLUMETA_DEMO_GSE": "GSE99999",
                "ILLUMETA_DEMO_GROUP_COLUMN": "disease state",
                "ILLUMETA_DEMO_GROUP_CON": "control",
                "ILLUMETA_DEMO_GROUP_TEST": "case",
            }
            with mock.patch.dict(os.environ, env, clear=False), \
                 mock.patch("illumeta.subprocess.run", side_effect=fake_run):
                with self.assertRaises(SystemExit):
                    illumeta.run_demo(_args(out_dir=demo_dir, offline=False))
            self.assertIn("GSE99999", calls[0])
            self.assertIn("disease state", calls[1])
            self.assertIn("case", calls[1])

    def test_custom_gse_without_grouping_exits_2(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            demo_dir = os.path.join(tmp, "projects", "demo")
            with mock.patch.dict(os.environ, {"ILLUMETA_DEMO_GSE": "GSE99999"}, clear=False), \
                 mock.patch.dict(os.environ, {k: "" for k in
                                 ("ILLUMETA_DEMO_GROUP_COLUMN", "ILLUMETA_DEMO_GROUP_CON",
                                  "ILLUMETA_DEMO_GROUP_TEST")}, clear=False), \
                 mock.patch("illumeta.subprocess.run") as run:
                with self.assertRaises(SystemExit) as cm:
                    illumeta.run_demo(_args(out_dir=demo_dir, offline=False))
                self.assertEqual(cm.exception.code, 2)
            run.assert_not_called()

    def test_invalid_demo_gse_exits_2(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.dict(os.environ, {"ILLUMETA_DEMO_GSE": "not-a-gse"}, clear=False), \
                 mock.patch("illumeta.subprocess.run") as run:
                with self.assertRaises(SystemExit) as cm:
                    illumeta.run_demo(_args(out_dir=tmp, offline=False))
                self.assertEqual(cm.exception.code, 2)
            run.assert_not_called()

    def test_failure_detail_reads_failure_summary(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "failure_summary.json").write_text(
                json.dumps({"code": "AUTO_GROUP_FAILED", "message": "column not found"}),
                encoding="utf-8",
            )
            detail = illumeta._demo_failure_detail(tmp)
            self.assertIn("AUTO_GROUP_FAILED", detail)
            self.assertIn("column not found", detail)
            self.assertEqual(illumeta._demo_failure_detail("/nonexistent_dir_xyz"), "")


if __name__ == "__main__":
    unittest.main()
