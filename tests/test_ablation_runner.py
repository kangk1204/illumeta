"""Tests for ablation runner state reuse contracts."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import ablation_runner


class AblationRunnerReuseTests(unittest.TestCase):
    def write_metric_csv(self, path: str) -> None:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
            writer.writeheader()
            writer.writerow({"metric": "pipeline", "value": "Minfi"})
            writer.writerow({"metric": "lambda", "value": "1.02"})

    def test_summary_only_variant_is_not_reusable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            variant_dir = os.path.join(tmpdir, "baseline")
            os.makedirs(variant_dir)
            with open(os.path.join(variant_dir, "summary.json"), "w", encoding="utf-8") as handle:
                json.dump({"n_con": 2, "n_test": 2, "primary_branch": "Minfi"}, handle)

            reusable, reason = ablation_runner.classify_reusable_variant(variant_dir)
            self.assertFalse(reusable)
            self.assertEqual(reason, "missing_or_invalid_analysis_parameters")

    def test_complete_variant_is_reusable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            variant_dir = os.path.join(tmpdir, "baseline")
            os.makedirs(variant_dir)
            with open(os.path.join(variant_dir, "summary.json"), "w", encoding="utf-8") as handle:
                json.dump({"n_con": 2, "n_test": 2, "primary_branch": "Minfi"}, handle)
            with open(os.path.join(variant_dir, "analysis_parameters.json"), "w", encoding="utf-8") as handle:
                json.dump({"group_con": "Control", "group_test": "Case"}, handle)
            self.write_metric_csv(os.path.join(variant_dir, "Minfi_Metrics.csv"))

            reusable, reason = ablation_runner.classify_reusable_variant(variant_dir)
            self.assertTrue(reusable)
            self.assertEqual(reason, "complete")

    def test_reuse_requires_metrics_for_declared_primary_branch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            variant_dir = os.path.join(tmpdir, "baseline")
            os.makedirs(variant_dir)
            with open(os.path.join(variant_dir, "summary.json"), "w", encoding="utf-8") as handle:
                json.dump({"n_con": 2, "n_test": 2, "primary_branch": "Sesame"}, handle)
            with open(os.path.join(variant_dir, "analysis_parameters.json"), "w", encoding="utf-8") as handle:
                json.dump({"group_con": "Control", "group_test": "Case"}, handle)
            self.write_metric_csv(os.path.join(variant_dir, "Minfi_Metrics.csv"))

            reusable, reason = ablation_runner.classify_reusable_variant(variant_dir)
            self.assertFalse(reusable)
            self.assertEqual(reason, "missing_or_invalid_primary_metrics")

    def test_reuse_existing_dry_run_does_not_collect_stale_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            variant_dir = os.path.join(tmpdir, "baseline")
            os.makedirs(variant_dir)
            with open(os.path.join(variant_dir, "summary.json"), "w", encoding="utf-8") as handle:
                json.dump({"n_con": 2, "n_test": 2, "primary_branch": "Minfi"}, handle)

            result = subprocess.run(
                [
                    sys.executable,
                    "scripts/ablation_runner.py",
                    "--config",
                    os.path.join(tmpdir, "configure.tsv"),
                    "--group-con",
                    "Control",
                    "--group-test",
                    "Case",
                    "--out-root",
                    tmpdir,
                    "--variants",
                    "baseline",
                    "--reuse-existing",
                    "--dry-run",
                ],
                cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            with open(os.path.join(tmpdir, "ablation_manifest.tsv"), encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertTrue(rows[0]["status"].startswith("would_rerun_stale_or_incomplete("))
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "ablation_counts.tsv")))


if __name__ == "__main__":
    unittest.main()
