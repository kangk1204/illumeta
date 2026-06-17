"""Regression guard for the dementia paper's load-bearing numbers.

These pin the headline figures so an accidental edit to a committed result table,
or a meta re-run that silently changes a number, fails CI instead of sailing into
the manuscript. The canonical TSVs live under benchmarks/ (gitignored paper data);
when they are absent (clean checkout / CI) each test SKIPS rather than fails.
"""

from __future__ import annotations

import csv
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEM = os.path.join(ROOT, "benchmarks", "dementia_special_issue")
COHORT = os.path.join(DEM, "cross_cohort", "cohort_summary.tsv")
META = os.path.join(DEM, "meta_cpg_analysis_cli", "meta_branch_summary.tsv")
CONCORDANT = os.path.join(DEM, "meta_cpg_analysis_cli", "branch_concordant_core_candidates.tsv")
CELLTYPE = os.path.join(DEM, "cell_type_reframe_summary.tsv")


def _rows(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CanonicalNumberTests(unittest.TestCase):
    def test_aggregate_consensus_sums(self):
        if not os.path.isfile(COHORT):
            self.skipTest("cohort_summary.tsv not present (gitignored benchmark data)")
        rows = _rows(COHORT)
        strict = sum(int(r["consensus_strict_dmp_rows"]) for r in rows)
        native = sum(int(r["consensus_native_dmp_rows"]) for r in rows)
        n_con = sum(int(r["n_control"]) for r in rows)
        n_ad = sum(int(r["n_AD"]) for r in rows)
        self.assertEqual(strict, 2793, "aggregate strict consensus DMPs drifted")
        self.assertEqual(native, 3629, "aggregate native consensus DMPs drifted")
        self.assertEqual((n_con, n_ad), (771, 883), "cohort sample totals drifted")

    def test_meta_branch_summary(self):
        if not os.path.isfile(META):
            self.skipTest("meta_branch_summary.tsv not present")
        by = {r["branch"]: r for r in _rows(META)}
        expected_analyzable = {"minfi": 330686, "sesame_strict": 253277, "sesame_native": 292045}
        expected_core = {"minfi": 1273, "sesame_strict": 1554, "sesame_native": 1721}
        for branch, val in expected_analyzable.items():
            self.assertEqual(int(by[branch]["n_meta_analyzable_cpgs"]), val, f"{branch} analyzable drifted")
        for branch, val in expected_core.items():
            self.assertEqual(int(by[branch]["n_core_candidates"]), val, f"{branch} core candidates drifted")

    def test_branch_concordant_counts_and_lead(self):
        if not os.path.isfile(CONCORDANT):
            self.skipTest("branch_concordant_core_candidates.tsv not present")
        rows = _rows(CONCORDANT)
        self.assertEqual(len(rows), 2995, "core-candidate union drifted")
        ge2 = [r for r in rows if int(r["n_candidate_branches"]) >= 2]
        all3 = [r for r in rows if int(r["n_candidate_branches"]) == 3]
        self.assertEqual(len(ge2), 1333, ">=2-branch concordant count drifted")
        self.assertEqual(len(all3), 220, "all-3-branch count drifted")
        up = sum(1 for r in ge2 if r["direction"] == "up")
        down = sum(1 for r in ge2 if r["direction"] == "down")
        self.assertEqual((up, down), (704, 629), "concordant up/down split drifted")
        # Lead three-branch candidate by smallest max FDR.
        three = sorted(all3, key=lambda r: float(r["max_random_fdr"]))
        self.assertEqual(three[0]["CpG"], "cg11979621", "lead 3-branch CpG drifted")
        self.assertEqual(three[0]["Gene"], "RBM33", "lead 3-branch gene drifted")

    def test_celltype_reference_axis(self):
        if not os.path.isfile(CELLTYPE):
            self.skipTest("cell_type_reframe_summary.tsv not present")
        rows = _rows(CELLTYPE)
        gse306226 = next((r for r in rows if "GSE306226" in r.get("source", "") or "GSE306226" in r.get("axis", "")), None)
        self.assertIsNotNone(gse306226, "GSE306226 reference row not found")
        self.assertEqual(int(gse306226["strict_consensus_dmps"]), 6398, "GSE306226 strict drifted")
        self.assertEqual(int(gse306226["native_consensus_dmps"]), 6269, "GSE306226 native drifted")


if __name__ == "__main__":
    unittest.main()
