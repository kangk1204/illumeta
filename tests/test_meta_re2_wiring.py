"""The RE2 flag has to reach the output table, not just exist.

`illumeta_re2.py` is covered by its own unit tests, but a correct module wired up wrongly
produces exactly the same green suite and an empty column in the file a user reads. These
tests exercise the path from the argparse flag through the meta pipeline to the header and
the rows, and pin the two properties that make the column safe to read: it is absent unless
asked for, and it never displaces an existing column.
"""
from __future__ import annotations

import csv
import gzip
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from illumeta_meta import _base_fieldnames, run_meta_cli


def write_branch_table(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["CpG", "logFC", "SE", "t", "P.Value", "adj.P.Val", "Delta_Beta",
                    "chr", "pos", "Gene", "Region", "Island_Context"])
        for cpg, lfc, se in rows:
            w.writerow([cpg, lfc, se, lfc / se if se else 0, 0.01, 0.05, lfc / 10.0,
                        "chr1", 1000 + 100 * int(cpg[2:]), "GENE", "Body", "Island"])


def make_cohort(root: Path, name: str, effects, n_con=30, n_test=30):
    res = root / name / "AD_vs_Control_results"
    res.mkdir(parents=True, exist_ok=True)
    (res / "summary.json").write_text(
        '{"primary_result_mode": "standard", "n_control": %d, "n_test": %d}' % (n_con, n_test)
    )
    write_branch_table(res / "Minfi_DMPs_full.csv", effects)
    return res


class FieldnameTests(unittest.TestCase):
    def test_columns_absent_unless_requested(self):
        default = _base_fieldnames(["a", "b"])
        for column in ("re2_stat", "re2_p", "re2_fdr"):
            self.assertNotIn(column, default)

    def test_columns_are_added_not_reordered(self):
        default = _base_fieldnames(["a", "b"])
        requested = _base_fieldnames(["a", "b"], report_re2=True)
        re2 = {"re2_stat", "re2_p", "re2_fdr"}
        for column in re2:
            self.assertIn(column, requested)
        self.assertEqual([c for c in requested if c not in re2], default)
        self.assertEqual(len(requested), len(default) + 3)

    def test_re2_and_knapp_hartung_compose(self):
        """Both flags at once must give both blocks, and still only add."""
        default = _base_fieldnames(["a", "b"])
        both = _base_fieldnames(["a", "b"], report_knapp_hartung=True, report_re2=True)
        extra = {"random_se_hk", "random_p_hk", "random_fdr_hk",
                 "re2_stat", "re2_p", "re2_fdr"}
        self.assertEqual([c for c in both if c not in extra], default)
        self.assertEqual(len(both), len(default) + 6)


class EndToEndTests(unittest.TestCase):
    """Run the real CLI entry point and read the file it writes."""

    def _run(self, tmp: Path, report_re2: bool):
        cohorts = [
            make_cohort(tmp, "C1", [("cg01", 0.30, 0.05), ("cg02", 0.02, 0.05),
                                    ("cg03", 0.40, 0.05)]),
            make_cohort(tmp, "C2", [("cg01", 0.28, 0.05), ("cg02", -0.02, 0.05),
                                    ("cg03", -0.40, 0.05)]),
            make_cohort(tmp, "C3", [("cg01", 0.32, 0.05), ("cg02", 0.01, 0.05),
                                    ("cg03", 0.40, 0.05)]),
        ]
        out = tmp / ("out_re2" if report_re2 else "out_plain")
        args = SimpleNamespace(
            result_dirs=[str(c) for c in cohorts], manifest=None, output=str(out),
            project_root=str(tmp), branches="minfi", branch_file=[],
            allow_missing_branches=False, allow_missing_summary=True,
            tier3_primary=True, allow_missing_tier3_primary=True,
            report_knapp_hartung=False, report_re2=report_re2, use_bacon=False,
            min_cohorts=3, meta_fdr=0.05, min_direction_fraction=0.70,
            min_loo_direction_fraction=0.80, max_i2=60.0, min_abs_delta_beta=0.0,
            partial_conjunction_r=3, top_n=100,
        )
        self.assertEqual(run_meta_cli(args), 0)
        with gzip.open(out / "minfi_meta_full.tsv.gz", "rt") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_flag_off_leaves_no_re2_columns_in_the_written_file(self):
        with tempfile.TemporaryDirectory() as td:
            rows = self._run(Path(td), report_re2=False)
            self.assertTrue(rows)
            for column in ("re2_stat", "re2_p", "re2_fdr"):
                self.assertNotIn(column, rows[0])

    def test_flag_on_writes_populated_re2_columns(self):
        with tempfile.TemporaryDirectory() as td:
            rows = self._run(Path(td), report_re2=True)
            self.assertTrue(rows)
            for column in ("re2_stat", "re2_p", "re2_fdr"):
                self.assertIn(column, rows[0])
            # The column must carry numbers. An empty or all-NaN column is the exact
            # failure a module-only test cannot see.
            values = [float(r["re2_p"]) for r in rows if r["re2_p"] not in ("", "nan")]
            self.assertTrue(values)
            self.assertTrue(all(0.0 <= v <= 1.0 for v in values))
            stats = [float(r["re2_stat"]) for r in rows if r["re2_stat"] not in ("", "nan")]
            self.assertTrue(all(s >= 0.0 for s in stats))

    def test_the_flag_changes_nothing_else(self):
        """RE2 is a reported diagnostic; it must not move the pooled estimate."""
        with tempfile.TemporaryDirectory() as td:
            plain = self._run(Path(td), report_re2=False)
        with tempfile.TemporaryDirectory() as td:
            withre2 = self._run(Path(td), report_re2=True)
        self.assertEqual(len(plain), len(withre2))
        for a, b in zip(plain, withre2):
            self.assertEqual(a["CpG"], b["CpG"])
            for column in ("random_effect_logFC", "random_se", "random_p",
                           "random_fdr", "core_candidate", "I2", "tau2"):
                self.assertEqual(a[column], b[column], column)

    def test_a_discordant_site_separates_the_two_references(self):
        """cg03 has opposite signs across cohorts: no mean effect, real heterogeneity.

        This is the case RE2 exists for, and it is the one that proves the column is
        computed from the data rather than copied from random_p.
        """
        with tempfile.TemporaryDirectory() as td:
            rows = self._run(Path(td), report_re2=True)
        by_cpg = {r["CpG"]: r for r in rows}
        self.assertIn("cg03", by_cpg)
        discordant = by_cpg["cg03"]
        self.assertGreater(float(discordant["random_p"]), 0.05)
        self.assertLess(float(discordant["re2_p"]), float(discordant["random_p"]))


if __name__ == "__main__":
    unittest.main()
