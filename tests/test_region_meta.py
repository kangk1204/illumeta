"""Tests for region-level cross-cohort meta-analysis.

The grouping rules are the whole method: a region is defined by where the runs break, so
each break condition gets its own test with a hand-built case whose answer is obvious by
inspection. The statistic and the empirical calibration are checked against values computed
by hand rather than against the implementation's own output.
"""
from __future__ import annotations

import math
import unittest

from illumeta_region import (
    CpGRecord,
    Region,
    calibrate_regions,
    call_regions,
    combine_region,
    summarize,
)


def rec(cpg, chrom, pos, effect, se=0.1, p=0.001, delta=0.02, gene=""):
    return CpGRecord(cpg=cpg, chrom=chrom, pos=pos, effect=effect, se=se, p=p,
                     delta_beta=delta, gene=gene)


class TestGrouping(unittest.TestCase):
    def test_consecutive_same_sign_within_gap_form_one_region(self):
        recs = [rec(f"cg{i}", "chr1", 1000 + 100 * i, 0.3) for i in range(4)]
        regions = call_regions(recs, max_gap=500)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].n_cpgs, 4)
        self.assertEqual((regions[0].start, regions[0].end), (1000, 1300))
        self.assertEqual(regions[0].direction, "up")

    def test_gap_larger_than_max_gap_splits(self):
        recs = [rec("a", "chr1", 1000, 0.3), rec("b", "chr1", 1100, 0.3),
                rec("c", "chr1", 9000, 0.3), rec("d", "chr1", 9100, 0.3)]
        regions = call_regions(recs, max_gap=500)
        self.assertEqual([r.n_cpgs for r in regions], [2, 2])

    def test_sign_change_splits(self):
        recs = [rec("a", "chr1", 1000, 0.3), rec("b", "chr1", 1100, 0.3),
                rec("c", "chr1", 1200, -0.3), rec("d", "chr1", 1300, -0.3)]
        regions = call_regions(recs, max_gap=500)
        self.assertEqual([r.direction for r in regions], ["up", "down"])

    def test_chromosome_change_splits(self):
        recs = [rec("a", "chr1", 1000, 0.3), rec("b", "chr1", 1100, 0.3),
                rec("c", "chr2", 1200, 0.3), rec("d", "chr2", 1300, 0.3)]
        regions = call_regions(recs, max_gap=500)
        self.assertEqual([r.chrom for r in regions], ["chr1", "chr2"])

    def test_cpg_failing_seed_p_breaks_the_run(self):
        recs = [rec("a", "chr1", 1000, 0.3), rec("b", "chr1", 1100, 0.3),
                rec("x", "chr1", 1150, 0.3, p=0.9),
                rec("c", "chr1", 1200, 0.3), rec("d", "chr1", 1300, 0.3)]
        regions = call_regions(recs, max_gap=500, seed_p=0.05)
        self.assertEqual(len(regions), 2)
        self.assertNotIn("x", [c for r in regions for c in r.cpgs])

    def test_singleton_runs_are_dropped(self):
        recs = [rec("a", "chr1", 1000, 0.3), rec("b", "chr1", 9000, 0.3)]
        self.assertEqual(call_regions(recs, max_gap=500), [])

    def test_min_cpgs_is_honoured(self):
        recs = [rec(f"cg{i}", "chr1", 1000 + 100 * i, 0.3) for i in range(2)]
        self.assertEqual(len(call_regions(recs, min_cpgs=2)), 1)
        self.assertEqual(call_regions(recs, min_cpgs=3), [])

    def test_unusable_records_are_ignored(self):
        recs = [rec("a", "chr1", 1000, 0.3), rec("b", "chr1", 1100, 0.3),
                CpGRecord("bad", "chr1", 1150, float("nan"), 0.1, 0.001),
                CpGRecord("zero_se", "chr1", 1160, 0.3, 0.0, 0.001)]
        regions = call_regions(recs, max_gap=500)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].n_cpgs, 2)

    def test_input_order_does_not_matter(self):
        forward = [rec(f"cg{i}", "chr1", 1000 + 100 * i, 0.3) for i in range(4)]
        shuffled = [forward[2], forward[0], forward[3], forward[1]]
        self.assertEqual(
            [r.as_row() for r in call_regions(forward)],
            [r.as_row() for r in call_regions(shuffled)],
        )


class TestStatistic(unittest.TestCase):
    def test_summed_effect_and_independent_se(self):
        members = [rec("a", "chr1", 1, 0.2, se=0.1), rec("b", "chr1", 2, 0.4, se=0.1)]
        estimate, se, z = combine_region(members)
        self.assertAlmostEqual(estimate, 0.6)
        self.assertAlmostEqual(se, math.sqrt(0.01 + 0.01))
        self.assertAlmostEqual(z, 0.6 / math.sqrt(0.02))

    def test_degenerate_variance_yields_nan_not_an_exception(self):
        members = [CpGRecord("a", "chr1", 1, 0.2, 0.0, 0.01)]
        self.assertTrue(all(math.isnan(v) for v in combine_region(members)))

    def test_more_cpgs_raise_z_when_effects_agree(self):
        two = combine_region([rec("a", "chr1", 1, 0.2), rec("b", "chr1", 2, 0.2)])
        four = combine_region([rec(f"c{i}", "chr1", i, 0.2) for i in range(4)])
        self.assertGreater(abs(four[2]), abs(two[2]))


class TestCalibration(unittest.TestCase):
    def _region(self, z):
        return Region(chrom="chr1", start=1, end=2, n_cpgs=2, cpgs=("a", "b"),
                      genes=(), direction="up", estimate=0.1, se=0.1 / max(abs(z), 1e-9),
                      z=z, p_independent=1.0, mean_delta_beta=0.02, min_cpg_p=0.001)

    def test_empirical_p_counts_null_statistics_at_least_as_extreme(self):
        obs = [self._region(3.0)]
        nulls = [[1.0, 2.0], [1.5, 4.0]]          # one null z (4.0) exceeds 3.0
        calibrate_regions(obs, nulls)
        self.assertAlmostEqual(obs[0].p_empirical, (1 + 1) / (4 + 1))
        self.assertFalse(obs[0].p_empirical_is_bound)

    def test_p_hits_the_floor_and_is_flagged_when_no_null_is_as_extreme(self):
        obs = [self._region(9.0)]
        nulls = [[1.0, 2.0], [1.5, 3.0]]
        calibrate_regions(obs, nulls)
        self.assertAlmostEqual(obs[0].p_empirical, 1 / (4 + 1))
        self.assertTrue(obs[0].p_empirical_is_bound)

    def test_sign_of_z_does_not_matter(self):
        up, down = self._region(3.0), self._region(-3.0)
        calibrate_regions([up, down], [[1.0, 4.0], [1.5, 2.0]])
        self.assertAlmostEqual(up.p_empirical, down.p_empirical)

    def test_fdr_is_expected_false_over_observed_at_the_same_threshold(self):
        obs = [self._region(5.0), self._region(4.0), self._region(3.0)]
        nulls = [[3.5], [2.0]]   # 2 replicates; one null z exceeds 3.0
        calibrate_regions(obs, nulls)
        # at |z| >= 3.0 there are 3 observed regions and 1 null over 2 replicates
        self.assertAlmostEqual(obs[2].fdr_empirical, (1 / 2) / 3)
        # at |z| >= 5.0 no null reaches it, so the estimated false count is zero
        self.assertAlmostEqual(obs[0].fdr_empirical, 0.0)

    def test_fdr_is_capped_at_one(self):
        obs = [self._region(1.0)]
        nulls = [[5.0, 5.0, 5.0]]
        calibrate_regions(obs, nulls)
        self.assertLessEqual(obs[0].fdr_empirical, 1.0)

    def test_empty_null_leaves_regions_untouched(self):
        obs = [self._region(3.0)]
        calibrate_regions(obs, [])
        self.assertTrue(math.isnan(obs[0].p_empirical))

    def test_more_extreme_region_never_gets_a_larger_p(self):
        obs = [self._region(2.0), self._region(6.0)]
        calibrate_regions(obs, [[1.0, 3.0, 5.0], [2.5, 4.0, 7.0]])
        self.assertLessEqual(obs[1].p_empirical, obs[0].p_empirical)


class TestSummary(unittest.TestCase):
    def test_counts_only_regions_below_the_fdr_threshold(self):
        regions = [
            Region("chr1", 1, 200, 3, ("a", "b", "c"), ("G1",), "up", 0.3, 0.1, 3.0,
                   0.002, 0.02, 0.001, 0.01, False, 0.01),
            Region("chr2", 1, 100, 2, ("d", "e"), (), "down", -0.2, 0.1, -2.0,
                   0.04, -0.03, 0.01, 0.4, False, 0.4),
        ]
        out = summarize(regions, fdr_threshold=0.05)
        self.assertEqual(out["n_candidate_regions"], 2)
        self.assertEqual(out["n_regions_fdr"], 1)
        self.assertEqual((out["n_up"], out["n_down"]), (1, 0))
        self.assertEqual(out["median_n_cpgs"], 3)

    def test_empty_input_reports_nan_medians_not_a_crash(self):
        out = summarize([])
        self.assertEqual(out["n_regions_fdr"], 0)
        self.assertTrue(math.isnan(out["median_n_cpgs"]))


if __name__ == "__main__":
    unittest.main()


class TestCliWiring(unittest.TestCase):
    """The module is reachable from the command line and its guards actually fire.

    A region module with passing unit tests and no call path is an analysis script, not a
    feature. These tests run the real entry point on synthetic cohorts and check the two
    refusals that keep an uncalibrated table from ever being written.
    """

    def _cohorts(self, root, n_cohorts, n_cpgs=40):
        import csv as _csv
        dirs = []
        for c in range(n_cohorts):
            res = root / f"C{c}" / "AD_vs_Control_results"
            res.mkdir(parents=True, exist_ok=True)
            (res / "summary.json").write_text(
                '{"primary_result_mode": "standard", "n_control": 30, "n_test": 30}'
            )
            with (res / "Minfi_DMPs_full.csv").open("w", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["CpG", "logFC", "SE", "t", "P.Value", "adj.P.Val",
                            "Delta_Beta", "chr", "pos", "Gene", "Region", "Island_Context"])
                for i in range(n_cpgs):
                    # A contiguous block of concordant CpGs, then unrelated noise.
                    lfc = 0.30 if i < 12 else (0.01 if i % 2 else -0.01)
                    w.writerow([f"cg{i:04d}", lfc, 0.05, lfc / 0.05, 0.001, 0.01,
                                lfc / 10, "chr1", 1000 + 100 * i, "GENE", "Body", "Island"])
            dirs.append(str(res))
        return dirs

    def _args(self, root, dirs, **over):
        from types import SimpleNamespace
        base = dict(result_dirs=dirs, manifest=None, output=str(root / "out"),
                    project_root=str(root), branches="minfi",
                    allow_missing_branches=False, allow_missing_summary=True,
                    tier3_primary=True, allow_missing_tier3_primary=True,
                    max_gap=500, seed_p=0.05, min_cpgs=2, min_cohorts=3,
                    region_fdr=0.05, min_null_patterns=7, max_null_patterns=127)
        base.update(over)
        return SimpleNamespace(**base)

    def test_runs_end_to_end_and_writes_a_summary(self):
        import csv as _csv
        import tempfile
        from pathlib import Path as P
        from illumeta_region import run_regions_cli
        with tempfile.TemporaryDirectory() as td:
            root = P(td)
            dirs = self._cohorts(root, 4)
            self.assertEqual(run_regions_cli(self._args(root, dirs)), 0)
            summary = list(_csv.DictReader((root / "out" / "region_summary.tsv").open(),
                                           delimiter="\t"))
            self.assertEqual(len(summary), 1)
            self.assertGreater(int(summary[0]["n_candidate_regions"]), 0)
            self.assertEqual(int(summary[0]["n_null_patterns"]), 2 ** (4 - 1) - 1)
            self.assertTrue((root / "out" / "region_manifest.json").is_file())
            self.assertTrue((root / "out" / "minfi_regions.tsv").is_file())

    def test_too_few_cohorts_is_refused_not_run_with_a_thin_null(self):
        import tempfile
        from pathlib import Path as P
        from illumeta_region import run_regions_cli
        with tempfile.TemporaryDirectory() as td:
            root = P(td)
            dirs = self._cohorts(root, 3)          # 3 cohorts -> 3 sign patterns
            with self.assertRaises(ValueError) as ctx:
                run_regions_cli(self._args(root, dirs))
            self.assertIn("sign patterns", str(ctx.exception))

    def test_too_many_cohorts_is_refused_until_the_limit_is_raised(self):
        import tempfile
        from pathlib import Path as P
        from illumeta_region import run_regions_cli
        with tempfile.TemporaryDirectory() as td:
            root = P(td)
            dirs = self._cohorts(root, 5)
            with self.assertRaises(ValueError) as ctx:
                run_regions_cli(self._args(root, dirs, max_null_patterns=8))
            self.assertIn("max-null-patterns", str(ctx.exception))
            # and passes once the operator raises it deliberately
            self.assertEqual(
                run_regions_cli(self._args(root, dirs, max_null_patterns=31)), 0
            )

    def test_no_cohorts_is_refused(self):
        import tempfile
        from pathlib import Path as P
        from illumeta_region import run_regions_cli
        with tempfile.TemporaryDirectory() as td:
            root = P(td)
            with self.assertRaises(ValueError):
                run_regions_cli(self._args(root, []))
