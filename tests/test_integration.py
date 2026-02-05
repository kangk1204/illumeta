"""
Integration tests for IlluMeta pipeline.

These tests verify end-to-end functionality of the pipeline components
without requiring actual IDAT files. For full integration tests with
real data, use the smoke test framework in scripts/run_smoke_pipeline.py.
"""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from illumeta import (
    auto_group_config,
    preflight_analysis,
    safe_path_component,
    load_config_rows,
    normalize_group_value,
    score_group_candidate,
    collect_group_candidates,
    collect_dashboard_warnings,
)


def write_config(path, headers, rows):
    """Helper to write a TSV config file."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


def create_mock_idat_pair(idat_dir, basename):
    """Create empty mock IDAT files for testing."""
    os.makedirs(idat_dir, exist_ok=True)
    for suffix in ("_Grn.idat", "_Red.idat"):
        path = os.path.join(idat_dir, basename + suffix)
        open(path, "wb").close()


class TestSafePathComponent(unittest.TestCase):
    """Tests for filesystem-safe path generation."""

    def test_ascii_passthrough(self):
        self.assertEqual(safe_path_component("test_name"), "test_name")

    def test_unicode_normalization(self):
        # Korean characters should be stripped, and leading _ is trimmed
        result = safe_path_component("테스트_test")
        self.assertEqual(result, "test")

    def test_space_replacement(self):
        self.assertEqual(safe_path_component("test name"), "test_name")

    def test_special_chars_replacement(self):
        result = safe_path_component("test@#$%name")
        self.assertEqual(result, "test_name")

    def test_empty_fallback(self):
        self.assertEqual(safe_path_component("", "fallback"), "fallback")
        self.assertEqual(safe_path_component("@#$%", "fallback"), "fallback")


class TestLoadConfigRows(unittest.TestCase):
    """Tests for config file loading."""

    def test_tsv_loading(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group", "age"]
            rows = [
                {"Basename": "S1", "primary_group": "control", "age": "30"},
                {"Basename": "S2", "primary_group": "case", "age": "45"},
            ]
            write_config(config_path, headers, rows)

            loaded_rows, loaded_headers, delim = load_config_rows(config_path)
            self.assertEqual(delim, "\t")
            self.assertEqual(len(loaded_rows), 2)
            self.assertIn("Basename", loaded_headers)

    def test_missing_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group", "optional"]
            rows = [
                {"Basename": "S1", "primary_group": "control", "optional": ""},
            ]
            write_config(config_path, headers, rows)

            loaded_rows, _, _ = load_config_rows(config_path)
            self.assertEqual(loaded_rows[0].get("optional"), "")


class TestNormalizeGroupValue(unittest.TestCase):
    """Tests for group value normalization."""

    def test_lowercase_strip(self):
        self.assertEqual(normalize_group_value("  Control  "), "control")

    def test_remove_special_chars(self):
        self.assertEqual(normalize_group_value("case-control"), "casecontrol")

    def test_none_handling(self):
        self.assertEqual(normalize_group_value(None), "")


class TestScoreGroupCandidate(unittest.TestCase):
    """Tests for auto-group scoring logic."""

    def test_binary_group_scores_high(self):
        values = ["control", "case", "control", "case"]
        result = score_group_candidate("disease_state", values)
        self.assertIsNotNone(result)
        self.assertGreater(result["score"], 0)

    def test_single_value_rejected(self):
        values = ["control", "control", "control"]
        result = score_group_candidate("constant_col", values)
        self.assertIsNone(result)

    def test_high_cardinality_rejected(self):
        values = [f"sample_{i}" for i in range(20)]
        result = score_group_candidate("sample_id", values)
        self.assertIsNone(result)


class TestCollectDashboardWarnings(unittest.TestCase):
    """Tests for dashboard warning generation."""

    def test_cross_reactive_skipped_warning(self):
        stats = {}
        analysis_params = {"unsafe_skip_cross_reactive": True}
        warnings = collect_dashboard_warnings(stats, analysis_params, None, None)
        self.assertTrue(any("cross-reactive" in w.lower() for w in warnings))

    def test_small_sample_tier_warning(self):
        stats = {}
        analysis_params = {"crf_sample_tier": "minimal"}
        warnings = collect_dashboard_warnings(stats, analysis_params, None, None)
        self.assertTrue(any("minimal" in w.lower() for w in warnings))

    def test_no_signal_warning(self):
        stats = {
            "minfi_up": 0, "minfi_down": 0,
            "sesame_up": 0, "sesame_down": 0,
            "sesame_native_up": 0, "sesame_native_down": 0,
            "intersect_up": 0, "intersect_down": 0,
            "intersect_native_up": 0, "intersect_native_down": 0,
        }
        analysis_params = {}
        warnings = collect_dashboard_warnings(stats, analysis_params, None, None)
        self.assertTrue(any("no significant" in w.lower() for w in warnings))


class TestPreflightIntegration(unittest.TestCase):
    """Integration tests for preflight analysis."""

    def test_preflight_with_mock_idats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")

            # Create mock IDAT pairs
            for sample in ("S1_R01C01", "S2_R01C01", "S3_R01C01", "S4_R01C01"):
                create_mock_idat_pair(idat_dir, sample)

            # Create config
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                {"Basename": "S1_R01C01", "primary_group": "control"},
                {"Basename": "S2_R01C01", "primary_group": "control"},
                {"Basename": "S3_R01C01", "primary_group": "case"},
                {"Basename": "S4_R01C01", "primary_group": "case"},
            ]
            write_config(config_path, headers, rows)

            result = preflight_analysis(
                config_path=config_path,
                idat_dir=idat_dir,
                group_con="control",
                group_test="case",
                min_total_size=4,
            )

            self.assertEqual(result["group_con"], 2)
            self.assertEqual(result["group_test"], 2)
            self.assertEqual(result["sample_count"], 4)

    def test_preflight_filters_missing_idats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")

            # Create only partial IDAT pairs (S3 missing)
            for sample in ("S1_R01C01", "S2_R01C01"):
                create_mock_idat_pair(idat_dir, sample)

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                {"Basename": "S1_R01C01", "primary_group": "control"},
                {"Basename": "S2_R01C01", "primary_group": "case"},
                {"Basename": "S3_R01C01", "primary_group": "case"},  # Missing IDAT
            ]
            write_config(config_path, headers, rows)

            result = preflight_analysis(
                config_path=config_path,
                idat_dir=idat_dir,
                group_con="control",
                group_test="case",
                min_total_size=2,
                drop_missing_idat=True,
            )

            # Should filter out S3
            self.assertEqual(result["sample_count"], 2)
            self.assertEqual(result["missing_idat_pairs"], 1)


class TestAutoGroupIntegration(unittest.TestCase):
    """Integration tests for auto-group functionality."""

    def test_auto_group_full_workflow(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["primary_group", "disease_state", "Basename", "age"]
            rows = [
                {"primary_group": "", "disease_state": "healthy", "Basename": "S1", "age": "30"},
                {"primary_group": "", "disease_state": "healthy", "Basename": "S2", "age": "35"},
                {"primary_group": "", "disease_state": "disease", "Basename": "S3", "age": "40"},
                {"primary_group": "", "disease_state": "disease", "Basename": "S4", "age": "45"},
            ]
            write_config(config_path, headers, rows)

            out_path, info = auto_group_config(
                config_path=config_path,
                group_con="Control",
                group_test="Case",
                group_column="disease_state",
            )

            self.assertTrue(info["updated"])
            self.assertEqual(info["source"], "column:disease_state")

            # Verify output file
            loaded_rows, _, _ = load_config_rows(out_path)
            groups = [row["primary_group"] for row in loaded_rows]
            self.assertEqual(groups.count("Control"), 2)
            self.assertEqual(groups.count("Case"), 2)


class TestDoctorCommand(unittest.TestCase):
    """Tests for the doctor command."""

    @unittest.skipIf(
        os.environ.get("CI") == "true" or os.environ.get("SKIP_SLOW_TESTS"),
        "Skipping slow doctor test in CI"
    )
    def test_doctor_runs_without_crash(self):
        """Verify doctor command executes without errors."""
        illumeta_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "illumeta.py"
        )

        result = subprocess.run(
            [sys.executable, illumeta_path, "doctor", "--skip-pandoc"],
            capture_output=True,
            text=True,
            timeout=300,  # Extended timeout for R package checks
        )

        # Doctor should exit 0 if dependencies are installed
        # or exit 1 with clear error message
        self.assertIn(result.returncode, [0, 1])

        # Should not crash with Python exception
        self.assertNotIn("Traceback", result.stderr)


class TestVersionConsistency(unittest.TestCase):
    """Tests for version consistency across files."""

    def test_version_matches_citation(self):
        """Verify illumeta.py version matches CITATION.cff."""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # Get version from illumeta.py
        illumeta_path = os.path.join(base_dir, "illumeta.py")
        with open(illumeta_path, "r") as f:
            content = f.read()

        import re
        match = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', content)
        self.assertIsNotNone(match)
        py_version = match.group(1)

        # Get version from CITATION.cff
        citation_path = os.path.join(base_dir, "CITATION.cff")
        with open(citation_path, "r") as f:
            content = f.read()

        match = re.search(r'^version:\s*(.+)$', content, re.MULTILINE)
        self.assertIsNotNone(match)
        cff_version = match.group(1).strip()

        self.assertEqual(py_version, cff_version)


if __name__ == "__main__":
    unittest.main()
