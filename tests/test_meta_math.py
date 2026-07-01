import inspect
import unittest

import math

import illumeta_meta
from illumeta_meta import (
    _chi2_sf_even,
    _directional_pc_p,
    _normal_two_sided_p,
    _random_effect_meta_one,
    bh_fdr,
)


class MetaMathTests(unittest.TestCase):
    def test_even_chi_square_survival_matches_closed_forms(self):
        self.assertAlmostEqual(_chi2_sf_even(2.0, 2), math.exp(-1.0), places=15)
        self.assertAlmostEqual(_chi2_sf_even(4.0, 4), math.exp(-2.0) * 3.0, places=15)

    def test_even_chi_square_survival_stays_stable_near_large_mean(self):
        # For df=2000, stat=df is near the distribution center and must not
        # underflow to zero while summing the finite Poisson tail identity.
        value = _chi2_sf_even(2000.0, 2000)

        self.assertTrue(math.isfinite(value))
        self.assertGreater(value, 0.45)
        self.assertLess(value, 0.55)

    def test_directional_partial_conjunction_pays_for_direction_choice(self):
        value = _directional_pc_p([1.0], [0.04], [True], 1)

        self.assertAlmostEqual(value, 0.04, places=15)

    def test_directional_partial_conjunction_caps_two_direction_correction(self):
        value = _directional_pc_p([1.0], [0.90], [True], 1)

        self.assertAlmostEqual(value, 0.90, places=15)

    def test_directionless_zero_effect_supports_neither_direction(self):
        # A cohort with exactly logFC==0 is directionless: it must NOT count as strong
        # one-sided support for both the up and down families (which would anti-conservatively
        # report ~0.04 here). With strict sign handling it contributes ~1.0 (no support).
        value = _directional_pc_p([0.0], [0.04], [True], 1)

        self.assertGreater(value, 0.95)

    # --- ground-truth regression tests (expected values derived independently of the code) ---

    def test_bh_fdr_matches_hand_computed_values(self):
        # Reordering + step-up + monotone enforcement, checked against hand calculations
        # equivalent to R p.adjust(method="BH").
        self.assertEqual([round(v, 6) for v in bh_fdr([0.9, 0.001])], [0.9, 0.002])
        # Ties get the same adjusted value.
        self.assertEqual([round(v, 6) for v in bh_fdr([0.01, 0.01, 0.03])], [0.015, 0.015, 0.03])
        # A uniform ramp collapses to the max under step-up.
        self.assertTrue(all(abs(v - 0.05) < 1e-9 for v in bh_fdr([0.01, 0.02, 0.03, 0.04, 0.05])))

    def test_bh_fdr_excludes_nan_from_m_and_passes_it_through(self):
        out = bh_fdr([0.01, float("nan"), 0.02])
        self.assertTrue(math.isnan(out[1]))
        self.assertAlmostEqual(out[0], 0.02, places=12)  # m=2, not 3
        self.assertAlmostEqual(out[2], 0.02, places=12)

    def test_bh_fdr_is_monotone_nondecreasing_in_sorted_p(self):
        ps = [i / 200.0 for i in range(1, 201)]
        adj = bh_fdr(ps)
        pairs = sorted(zip(ps, adj))
        self.assertTrue(all(pairs[i][1] <= pairs[i + 1][1] + 1e-12 for i in range(len(pairs) - 1)))

    def test_dersimonian_laird_ground_truth_and_clamps(self):
        r3 = _random_effect_meta_one([-2.0, 0.0, 2.0], [1.0, 1.0, 1.0], [True, True, True])
        self.assertAlmostEqual(r3["Q"], 8.0, places=9)
        self.assertAlmostEqual(r3["tau2"], 3.0, places=9)   # (Q-df)/C = (8-2)/2
        self.assertAlmostEqual(r3["I2"], 75.0, places=6)    # (Q-df)/Q = 6/8
        # Q <= df must clamp tau2 to 0 (no negative variance) and I2 to 0.
        rc = _random_effect_meta_one([-1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [True, True, True])
        self.assertEqual(rc["tau2"], 0.0)
        # Single cohort: Q/I2/tau2 must be NaN (heterogeneity undefined at k<2).
        r1 = _random_effect_meta_one([1.0], [0.1], [True])
        self.assertTrue(math.isnan(r1["Q"]) and math.isnan(r1["I2"]) and math.isnan(r1["tau2"]))

    def test_two_sided_normal_p_underflow_and_nan(self):
        self.assertAlmostEqual(_normal_two_sided_p(1.959964), 0.05, places=6)
        self.assertEqual(_normal_two_sided_p(40.0), 0.0)
        self.assertTrue(math.isnan(_normal_two_sided_p(float("nan"))))

    def test_two_cohort_heterogeneity_excluded_from_summary_aggregate(self):
        # F1 regression: df=1 (k=2) I^2 is degenerate and must be kept out of the summary
        # median even when min_cohorts==2. Locked at source level because the aggregate has
        # no standalone unit entry point.
        src = inspect.getsource(illumeta_meta._analyze_branch)
        self.assertIn("max(3, thresholds.min_cohorts)", src)


if __name__ == "__main__":
    unittest.main()
