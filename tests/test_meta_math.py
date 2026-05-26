import math
import unittest

from illumeta_meta import _chi2_sf_even


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


if __name__ == "__main__":
    unittest.main()
