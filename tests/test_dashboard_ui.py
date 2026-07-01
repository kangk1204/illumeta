"""Regression tests for dashboard accessibility/usability improvements.

Renders generate_dashboard() on an empty output dir (tolerant path: missing
summary.json -> all-zeros stats) and asserts the a11y/design markers are present.
"""

import os
import sys
import tempfile
import unittest

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class DashboardUiTests(unittest.TestCase):
    def _render(self):
        sys.path.insert(0, BASE_DIR)
        try:
            import illumeta
        finally:
            if sys.path and sys.path[0] == BASE_DIR:
                sys.path.pop(0)
        tmp = tempfile.mkdtemp()
        out = os.path.join(tmp, "Case_vs_Control_results")
        os.makedirs(out)
        illumeta.generate_dashboard(out, "Case", "Control")
        with open(os.path.join(tmp, "Case_vs_Control_results_index.html"), encoding="utf-8") as fh:
            return fh.read()

    def test_pipeline_tabs_are_a_real_tablist(self):
        h = self._render()
        self.assertIn('role="tablist"', h)
        self.assertIn('<button type="button" class="nav-tab', h)  # real buttons, not <div onclick>
        self.assertIn('role="tab"', h)
        self.assertIn('aria-selected=', h)
        self.assertIn('role="tabpanel"', h)
        self.assertIn("initTablistKeys", h)  # arrow-key keyboard navigation

    def test_icon_buttons_have_accessible_names(self):
        h = self._render()
        self.assertIn('aria-label="Toggle dark mode"', h)
        self.assertIn('aria-label="Copy results path to clipboard"', h)
        self.assertIn('aria-label="Open results directory"', h)

    def test_focus_and_print_styles_exist(self):
        h = self._render()
        self.assertIn(":focus-visible", h)   # keyboard focus indicator
        self.assertIn("@media print", h)     # all tabs printable / flattened header

    def test_dark_mode_callout_uses_variables_not_hardcoded_hex(self):
        h = self._render()
        self.assertIn("var(--callout-good-bg)", h)
        # the old hardcoded light-green callout background must be gone
        self.assertNotIn("background:#eaf6ed", h)

    def test_no_redundant_hero_tag_row_and_pipeline_numbers_are_labeled(self):
        h = self._render()
        self.assertEqual(h.count('class="hero-tags"'), 0)  # duplicate skim row removed
        self.assertIn("Minfi <b>", h)                      # labeled, not "N | N | N"

    def test_base64_thumbnails_drop_the_noop_lazy_attribute(self):
        h = self._render()
        self.assertNotIn('loading="lazy"', h)  # no-op on inline data: URIs


if __name__ == "__main__":
    unittest.main()
