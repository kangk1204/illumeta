import unittest

import math

from illumeta_meta import _chi2_sf_even, _directional_pc_p


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


if __name__ == "__main__":
    unittest.main()
