import os
import tempfile
import unittest

from illumeta import preflight_analysis


def write_config(path, headers, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


class PreflightTests(unittest.TestCase):
    def test_preflight_ok(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")
            os.makedirs(idat_dir, exist_ok=True)
            for sample in ("S1_R01C01", "S2_R01C01"):
                open(os.path.join(idat_dir, f"{sample}_Grn.idat"), "wb").close()
                open(os.path.join(idat_dir, f"{sample}_Red.idat"), "wb").close()

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                {"Basename": "S1_R01C01", "primary_group": "control"},
                {"Basename": "S2_R01C01", "primary_group": "all"},
            ]
            write_config(config_path, headers, rows)

            result = preflight_analysis(
                config_path=config_path,
                idat_dir=idat_dir,
                group_con="control",
                group_test="all",
                min_total_size=2,
                id_column=None,
            )
            self.assertEqual(result["group_con"], 1)
            self.assertEqual(result["group_test"], 1)

    def test_preflight_missing_pairs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")
            os.makedirs(idat_dir, exist_ok=True)
            open(os.path.join(idat_dir, "S1_R01C01_Grn.idat"), "wb").close()

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [{"Basename": "S1_R01C01", "primary_group": "control"}]
            write_config(config_path, headers, rows)

            with self.assertRaises(ValueError):
                preflight_analysis(
                    config_path=config_path,
                    idat_dir=idat_dir,
                    group_con="control",
                    group_test="all",
                    min_total_size=1,
                    id_column=None,
                )


    def test_preflight_basename_with_idat_suffix(self):
        """Basename already ending with _Grn.idat should still resolve correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")
            os.makedirs(idat_dir, exist_ok=True)
            for sample in ("S1_R01C01", "S2_R01C01"):
                open(os.path.join(idat_dir, f"{sample}_Grn.idat"), "wb").close()
                open(os.path.join(idat_dir, f"{sample}_Red.idat"), "wb").close()

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                # Basename has _Grn.idat suffix — must be stripped to avoid double-suffix
                {"Basename": os.path.join(idat_dir, "S1_R01C01_Grn.idat"), "primary_group": "control"},
                {"Basename": os.path.join(idat_dir, "S2_R01C01_Red.idat.gz"), "primary_group": "test"},
            ]
            write_config(config_path, headers, rows)

            # Should NOT raise — suffix must be stripped before pair check
            result = preflight_analysis(
                config_path=config_path,
                idat_dir=idat_dir,
                group_con="control",
                group_test="test",
                min_total_size=2,
                id_column=None,
            )
            self.assertEqual(result["group_con"], 1)
            self.assertEqual(result["group_test"], 1)

    def test_preflight_filters_non_target_groups(self):
        """Samples with labels outside control/test are excluded and config rewritten."""
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")
            os.makedirs(idat_dir, exist_ok=True)
            for sample in ("S1_R01C01", "S2_R01C01", "S3_R01C01"):
                open(os.path.join(idat_dir, f"{sample}_Grn.idat"), "wb").close()
                open(os.path.join(idat_dir, f"{sample}_Red.idat"), "wb").close()

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                {"Basename": "S1_R01C01", "primary_group": "control"},
                {"Basename": "S2_R01C01", "primary_group": "test"},
                {"Basename": "S3_R01C01", "primary_group": "other"},
            ]
            write_config(config_path, headers, rows)

            result = preflight_analysis(
                config_path=config_path,
                idat_dir=idat_dir,
                group_con="control",
                group_test="test",
                min_total_size=2,
                id_column=None,
            )
            self.assertEqual(result["group_con"], 1)
            self.assertEqual(result["group_test"], 1)
            self.assertEqual(result["sample_count"], 2)
            self.assertTrue(result["config_path"].endswith("_groupfiltered.tsv"))
            self.assertTrue(any("Excluded 1 sample" in w for w in result["warnings"]))

    def test_preflight_normalizes_group_labels(self):
        """Group matching should tolerate underscores/hyphens vs spaces."""
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")
            os.makedirs(idat_dir, exist_ok=True)
            for sample in ("S1_R01C01", "S2_R01C01"):
                open(os.path.join(idat_dir, f"{sample}_Grn.idat"), "wb").close()
                open(os.path.join(idat_dir, f"{sample}_Red.idat"), "wb").close()

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                {"Basename": "S1_R01C01", "primary_group": "tumor_adjacent"},
                {"Basename": "S2_R01C01", "primary_group": "normal-control"},
            ]
            write_config(config_path, headers, rows)

            # CLI args use spaces, configure.tsv uses underscores/hyphens
            result = preflight_analysis(
                config_path=config_path,
                idat_dir=idat_dir,
                group_con="normal control",
                group_test="tumor adjacent",
                min_total_size=2,
                id_column=None,
            )
            self.assertEqual(result["group_con"], 1)
            self.assertEqual(result["group_test"], 1)

    def test_preflight_rejects_normalized_group_collisions(self):
        """Distinct config labels must not collapse into the same normalized group."""
        with tempfile.TemporaryDirectory() as tmpdir:
            idat_dir = os.path.join(tmpdir, "idat")
            os.makedirs(idat_dir, exist_ok=True)
            for sample in ("S1_R01C01", "S2_R01C01", "S3_R01C01"):
                open(os.path.join(idat_dir, f"{sample}_Grn.idat"), "wb").close()
                open(os.path.join(idat_dir, f"{sample}_Red.idat"), "wb").close()

            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["Basename", "primary_group"]
            rows = [
                {"Basename": "S1_R01C01", "primary_group": "case-control"},
                {"Basename": "S2_R01C01", "primary_group": "casecontrol"},
                {"Basename": "S3_R01C01", "primary_group": "test"},
            ]
            write_config(config_path, headers, rows)

            with self.assertRaisesRegex(ValueError, "Ambiguous normalized group labels"):
                preflight_analysis(
                    config_path=config_path,
                    idat_dir=idat_dir,
                    group_con="case-control",
                    group_test="test",
                    min_total_size=2,
                    id_column=None,
                )

    def test_strip_idat_suffix(self):
        """Verify strip_idat_suffix removes all IDAT suffix variants."""
        from illumeta import strip_idat_suffix
        self.assertEqual(strip_idat_suffix("/path/sample_Grn.idat"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample_Red.idat"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample_Grn.idat.gz"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample_Red.idat.gz"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample"), "/path/sample")
        self.assertEqual(strip_idat_suffix(""), "")

    def test_r_analysis_uses_same_group_normalization_contract(self):
        """R analysis group matching must stay aligned with Python preflight."""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        analyze_path = os.path.join(base_dir, "r_scripts", "analyze.R")
        with open(analyze_path, "r", encoding="utf-8") as handle:
            src = handle.read()

        self.assertIn("normalize_group_label <- function", src)
        self.assertIn("target_groups_norm <- normalize_group_label(targets$primary_group)", src)
        self.assertIn("con_norm <- normalize_group_label(group_con_in)", src)
        self.assertIn("test_norm <- normalize_group_label(group_test_in)", src)
        self.assertNotIn("target_groups_lower <- tolower(targets$primary_group)", src)

    def test_r_analysis_rejects_group_normalization_collisions(self):
        """R analysis must reject labels that normalize into ambiguous groups."""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        analyze_path = os.path.join(base_dir, "r_scripts", "analyze.R")
        with open(analyze_path, "r", encoding="utf-8") as handle:
            src = handle.read()

        self.assertIn("Control/test group labels must not normalize to empty values.", src)
        self.assertIn("identical(con_norm, test_norm)", src)
        self.assertIn("Ambiguous normalized group labels detected in configuration", src)
        self.assertIn("target_groups_norm <- normalize_group_label(targets$primary_group)", src)

    def test_r_pipeline_failure_returns_standard_contract(self):
        """run_pipeline failure branches must keep callers' $res/$summary contract."""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        analyze_path = os.path.join(base_dir, "r_scripts", "analyze.R")
        with open(analyze_path, "r", encoding="utf-8") as handle:
            src = handle.read()

        self.assertNotIn("return(branch_summary)", src)
        self.assertIn("return(list(res = NULL", src)
        self.assertIn("summary = branch_summary))", src)

    def test_r_download_validates_direct_gse_argument(self):
        """Direct R entrypoint must reject malformed GSE IDs before network work."""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        download_path = os.path.join(base_dir, "r_scripts", "download.R")
        with open(download_path, "r", encoding="utf-8") as handle:
            src = handle.read()

        validation_pos = src.index('grepl("^GSE[0-9]+$", gse_id)')
        network_pos = src.index("getGEO(gse_id")
        self.assertLess(validation_pos, network_pos)


if __name__ == "__main__":
    unittest.main()
