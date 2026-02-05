import unittest

from illumeta import collect_dashboard_warnings


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


if __name__ == "__main__":
    unittest.main()
