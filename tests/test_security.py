"""Security tests for input validation and HTML escaping."""

import re
import unittest


class TestGSEValidation(unittest.TestCase):
    """GSE ID must be GSE followed by digits only."""

    PATTERN = re.compile(r'^GSE\d+$')

    def test_valid_gse(self):
        self.assertIsNotNone(self.PATTERN.match("GSE12345"))
        self.assertIsNotNone(self.PATTERN.match("GSE1"))
        self.assertIsNotNone(self.PATTERN.match("GSE999999"))

    def test_rejects_injection(self):
        self.assertIsNone(self.PATTERN.match('GSE12345; rm -rf /'))
        self.assertIsNone(self.PATTERN.match('GSE12345\nGSE99999'))
        self.assertIsNone(self.PATTERN.match('../etc/passwd'))
        self.assertIsNone(self.PATTERN.match(''))
        self.assertIsNone(self.PATTERN.match('GSE'))
        self.assertIsNone(self.PATTERN.match('gse12345'))

    def test_geo_suppl_url_validates(self):
        from illumeta import geo_suppl_url

        url = geo_suppl_url("GSE12345")
        self.assertIn("GSE12nnn", url)
        self.assertIn("GSE12345", url)

        with self.assertRaises(ValueError):
            geo_suppl_url("GS")
        with self.assertRaises(ValueError):
            geo_suppl_url("")
        with self.assertRaises(ValueError):
            geo_suppl_url("notGSE123")


class TestRPackageNameValidation(unittest.TestCase):
    """Package names with shell/R metacharacters must be rejected."""

    PATTERN = re.compile(r'^[A-Za-z][A-Za-z0-9._]+$')

    def test_valid_package_names(self):
        for name in ["minfi", "sesame", "FlowSorted.Blood.EPIC", "data.table", "lme4"]:
            self.assertIsNotNone(self.PATTERN.match(name), f"Should accept: {name}")

    def test_rejects_injection(self):
        for name in ['"; system("id")', "pkg\ncat /etc/passwd", "", "1package", "-flag",
                      "a", "pkg;rm", "pkg$(cmd)"]:
            self.assertIsNone(self.PATTERN.match(name), f"Should reject: {name!r}")


class TestHTMLEscaping(unittest.TestCase):
    """_h() must properly escape HTML entities."""

    def test_escapes_script_tags(self):
        from illumeta import _h
        self.assertEqual(_h("<script>alert(1)</script>"),
                         "&lt;script&gt;alert(1)&lt;/script&gt;")

    def test_escapes_attribute_injection(self):
        from illumeta import _h
        result = _h('"><img onerror=alert(1)>')
        self.assertNotIn('">', result)
        self.assertIn("&lt;", result)

    def test_none_returns_na(self):
        from illumeta import _h
        self.assertEqual(_h(None), "N/A")

    def test_numeric(self):
        from illumeta import _h
        self.assertEqual(_h(42), "42")
        self.assertEqual(_h(3.14), "3.14")


class TestParseGroupMap(unittest.TestCase):
    """parse_group_map input validation."""

    def test_valid_comma_separated(self):
        from illumeta import parse_group_map
        m = parse_group_map("normal=Control,tumor=Case", "Control", "Case")
        self.assertEqual(m["normal"], "Control")
        self.assertEqual(m["tumor"], "Case")

    def test_valid_semicolon_separated(self):
        from illumeta import parse_group_map
        m = parse_group_map("normal=Control;tumor=Case", "Control", "Case")
        self.assertEqual(m["normal"], "Control")
        self.assertEqual(m["tumor"], "Case")

    def test_none_returns_empty(self):
        from illumeta import parse_group_map
        self.assertEqual(parse_group_map(None, "Control", "Case"), {})

    def test_empty_returns_empty(self):
        from illumeta import parse_group_map
        self.assertEqual(parse_group_map("", "Control", "Case"), {})


class TestSafeFloat(unittest.TestCase):
    """safe_float must handle edge cases without crashing."""

    def test_valid_float(self):
        from illumeta import safe_float
        self.assertEqual(safe_float("1.5"), 1.5)
        self.assertEqual(safe_float(3), 3.0)

    def test_na_variants(self):
        from illumeta import safe_float
        self.assertIsNone(safe_float("NA"))
        self.assertIsNone(safe_float(""))
        self.assertIsNone(safe_float(None))
        self.assertIsNone(safe_float("NaN"))
        self.assertIsNone(safe_float("not a number"))


class TestAsBool(unittest.TestCase):
    """_as_bool must handle various truthy/falsy inputs."""

    def test_true_values(self):
        from illumeta import _as_bool
        for val in [True, "true", "TRUE", "yes", "1", 1]:
            self.assertTrue(_as_bool(val), f"Expected True for {val!r}")

    def test_false_values(self):
        from illumeta import _as_bool
        for val in [False, "false", "FALSE", "no", "0", 0, "", None]:
            self.assertFalse(_as_bool(val), f"Expected False for {val!r}")


if __name__ == "__main__":
    unittest.main()
