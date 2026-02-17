import unittest

from illumeta import safe_int


class NumericUtilsTests(unittest.TestCase):
    def test_safe_int_accepts_float_like_strings(self):
        self.assertEqual(safe_int("1.0"), 1)
        self.assertEqual(safe_int("2.9"), 2)

    def test_safe_int_handles_missing_and_non_finite(self):
        self.assertEqual(safe_int(""), 0)
        self.assertEqual(safe_int("NA"), 0)
        self.assertEqual(safe_int("NaN"), 0)
        self.assertEqual(safe_int(float("inf")), 0)
        self.assertEqual(safe_int(float("-inf")), 0)


if __name__ == "__main__":
    unittest.main()
