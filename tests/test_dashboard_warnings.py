import unittest
from pathlib import Path

from illumeta import collect_dashboard_warnings

BASE_DIR = Path(__file__).resolve().parents[1]


class DashboardWarningsTests(unittest.TestCase):
    def test_warnings_basic(self):
        stats = {
            "primary_result_mode": "tier3_ineligible",
            "primary_tier3_low_power": False,
        }
        analysis_params = {
            "unsafe_skip_cross_reactive": True,
            "cross_reactive_active": False,
            "sex_check_action": "ignore",
            "sex_mismatch_count": 2,
            "auto_covariates_enabled": True,
            "cell_confound_eta2_threshold": 0.5,
        }
        cell_summary = [
            {"Method": "EpiDISH::RPC(centDHSbloodDMC.m)"},
        ]
        cell_assoc = [
            {"Eta_Squared": "0.8"},
        ]
        warnings = collect_dashboard_warnings(stats, analysis_params, {}, cell_summary, cell_assoc)
        joined = " | ".join(warnings)
        self.assertIn("Cross-reactive probe filtering was skipped", joined)
        self.assertIn("Sex mismatches detected", joined)
        self.assertIn("Auto-selected covariates", joined)
        self.assertIn("Tier3 confounding detected but eligibility failed", joined)
        self.assertIn("Reference-based cell-type deconvolution was used", joined)
        self.assertIn("Cell composition strongly associated with group", joined)

    def test_missing_cross_reactive_metadata_is_not_treated_as_explicit_skip(self):
        warnings = collect_dashboard_warnings(
            {"minfi_up": 1},
            {"pval_threshold": 0.05},
            {},
        )
        self.assertFalse(
            any("Cross-reactive probe" in warning for warning in warnings),
            warnings,
        )

    def test_explicit_cross_reactive_disable_is_reported(self):
        warnings = collect_dashboard_warnings(
            {"minfi_up": 1},
            {"cross_reactive_active": False},
            {},
        )
        self.assertIn(
            "Cross-reactive probe filtering was skipped (unsafe).",
            warnings,
        )

    def test_active_cross_reactive_filter_only_warns_for_explicit_empty_list(self):
        missing_count = collect_dashboard_warnings(
            {"minfi_up": 1},
            {"cross_reactive_active": True},
            {},
        )
        explicit_empty = collect_dashboard_warnings(
            {"minfi_up": 1},
            {"cross_reactive_active": True, "cross_reactive_count": 0},
            {},
        )
        self.assertFalse(
            any("probe list missing or empty" in warning for warning in missing_count),
            missing_count,
        )
        self.assertIn(
            "Cross-reactive probe list missing or empty; filtering may not have been applied.",
            explicit_empty,
        )

    def test_dashboard_effect_threshold_label_uses_m_value_scale(self):
        source = (BASE_DIR / "illumeta.py").read_text(encoding="utf-8")
        self.assertNotIn("minimum log2 fold-change in methylation M-values", source)
        self.assertIn("minimum absolute limma effect on the methylation M-value scale", source)

    def test_dashboard_results_path_escapes_script_close(self):
        source = (BASE_DIR / "illumeta.py").read_text(encoding="utf-8")
        self.assertIn('.replace("</", "<\\\\/")', source)

    def test_small_tier_does_not_mask_low_verdict_conditions(self):
        source = (BASE_DIR / "illumeta.py").read_text(encoding="utf-8")
        compute_block = source.split("def compute_verdict(", 1)[1].split("def pill_link(", 1)[0]
        self.assertLess(compute_block.index('tier == "minimal"'), compute_block.index('tier == "small"'))


if __name__ == "__main__":
    unittest.main()
