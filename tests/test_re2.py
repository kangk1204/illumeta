"""Tests for the RE2 modified-null random-effects test.

The statistic is a likelihood ratio, so the tests check the two likelihoods against values
computed directly from the normal density rather than against the optimiser's own answer, and
check the properties the test is supposed to have: monotone in the mean, monotone in the
heterogeneity, and equal to zero when the data sit exactly on the null.
"""
from __future__ import annotations

import math
import unittest

from illumeta_re2 import _loglik, _null_loglik, _chi2_sf, re2_one


def normal_loglik(xs, sds, mu):
    return sum(
        -0.5 * (math.log(2 * math.pi * s * s) + ((x - mu) / s) ** 2)
        for x, s in zip(xs, sds)
    )


class TestLikelihoods(unittest.TestCase):
    def test_null_loglik_matches_the_normal_density_at_zero(self):
        xs, sds = [0.2, -0.1, 0.05], [0.1, 0.1, 0.2]
        self.assertAlmostEqual(
            _null_loglik(xs, sds, [True] * 3), normal_loglik(xs, sds, 0.0), places=12
        )

    def test_loglik_at_zero_tau_matches_the_fixed_effect_fit(self):
        xs, sds = [0.2, 0.25, 0.22], [0.1, 0.1, 0.1]
        mu = sum(x / s**2 for x, s in zip(xs, sds)) / sum(1 / s**2 for s in sds)
        self.assertAlmostEqual(
            _loglik(xs, sds, [True] * 3, 0.0), normal_loglik(xs, sds, mu), places=12
        )

    def test_invalid_entries_are_excluded_not_imputed(self):
        xs, sds = [0.2, 0.25, 99.0], [0.1, 0.1, 0.1]
        with_all = _null_loglik(xs, sds, [True, True, True])
        without = _null_loglik(xs, sds, [True, True, False])
        self.assertAlmostEqual(without, normal_loglik(xs[:2], sds[:2], 0.0), places=12)
        self.assertLess(with_all, without)


class TestChiSquaredTails(unittest.TestCase):
    def test_one_df_tail_at_known_points(self):
        self.assertAlmostEqual(_chi2_sf(3.841458820694124, 1), 0.05, places=6)
        self.assertAlmostEqual(_chi2_sf(6.634896601021214, 1), 0.01, places=6)

    def test_two_df_tail_is_the_exponential(self):
        self.assertAlmostEqual(_chi2_sf(5.991464547107979, 2), 0.05, places=9)
        self.assertAlmostEqual(_chi2_sf(0.0, 2), 1.0)

    def test_nonpositive_statistic_gives_probability_one(self):
        self.assertEqual(_chi2_sf(0.0, 1), 1.0)
        self.assertEqual(_chi2_sf(-1.0, 2), 1.0)

    def test_unsupported_df_is_refused_rather_than_approximated(self):
        with self.assertRaises(ValueError):
            _chi2_sf(1.0, 3)


class TestRE2(unittest.TestCase):
    def test_data_exactly_on_the_null_gives_a_zero_statistic(self):
        out = re2_one([0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [True] * 3)
        self.assertAlmostEqual(out["re2_stat"], 0.0, places=9)
        self.assertAlmostEqual(out["re2_p"], 1.0, places=9)
        self.assertAlmostEqual(out["re2_tau2"], 0.0, places=9)

    def test_statistic_grows_with_a_consistent_mean_effect(self):
        small = re2_one([0.05] * 4, [0.1] * 4, [True] * 4)["re2_stat"]
        large = re2_one([0.50] * 4, [0.1] * 4, [True] * 4)["re2_stat"]
        self.assertGreater(large, small)

    def test_statistic_grows_with_heterogeneity_at_zero_mean(self):
        """A signal RE2 is meant to see and a plain fixed-effect test is not."""
        homogeneous = re2_one([0.0] * 4, [0.1] * 4, [True] * 4)["re2_stat"]
        opposed = re2_one([0.5, -0.5, 0.5, -0.5], [0.1] * 4, [True] * 4)["re2_stat"]
        self.assertGreater(opposed, homogeneous)
        self.assertGreater(opposed, 0.0)

    def test_opposite_direction_effects_are_invisible_to_the_pooled_mean(self):
        out = re2_one([0.5, -0.5, 0.5, -0.5], [0.1] * 4, [True] * 4)
        self.assertAlmostEqual(out["re2_mu"], 0.0, places=6)
        self.assertGreater(out["re2_tau2"], 0.0)
        self.assertLess(out["re2_p"], 0.05)

    def test_homogeneous_data_puts_tau2_on_the_boundary(self):
        out = re2_one([0.30, 0.30, 0.30], [0.1] * 3, [True] * 3)
        self.assertEqual(out["re2_tau2"], 0.0)

    def test_p_value_is_a_probability(self):
        for effects in ([0.0] * 3, [0.3] * 3, [1.0, -1.0, 1.0], [0.01, 0.02, -0.01]):
            p = re2_one(effects, [0.1] * 3, [True] * 3)["re2_p"]
            self.assertTrue(0.0 <= p <= 1.0, effects)

    def test_fewer_than_two_usable_studies_returns_nan_not_a_guess(self):
        out = re2_one([0.3, 0.3], [0.1, 0.1], [True, False])
        self.assertTrue(math.isnan(out["re2_stat"]))
        self.assertEqual(out["re2_k"], 1)

    def test_degenerate_standard_errors_are_dropped(self):
        out = re2_one([0.3, 0.3, 0.3], [0.1, 0.0, float("nan")], [True] * 3)
        self.assertEqual(out["re2_k"], 1)
        self.assertTrue(math.isnan(out["re2_p"]))

    def test_k_counts_only_usable_studies(self):
        out = re2_one([0.3, 0.3, 0.3, 0.3], [0.1, 0.1, -1.0, 0.1], [True, True, True, False])
        self.assertEqual(out["re2_k"], 2)

    def test_statistic_is_never_negative(self):
        for effects in ([1e-12] * 3, [0.0, 1e-9, -1e-9]):
            self.assertGreaterEqual(re2_one(effects, [0.1] * 3, [True] * 3)["re2_stat"], 0.0)

    def test_scaling_effects_and_errors_together_leaves_the_statistic_unchanged(self):
        """The likelihood ratio depends on effects only through effect/SE."""
        base = re2_one([0.2, 0.3, 0.25], [0.1, 0.1, 0.1], [True] * 3)["re2_stat"]
        scaled = re2_one([2.0, 3.0, 2.5], [1.0, 1.0, 1.0], [True] * 3)["re2_stat"]
        self.assertAlmostEqual(base, scaled, places=8)


if __name__ == "__main__":
    unittest.main()
