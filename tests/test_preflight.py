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

    def test_strip_idat_suffix(self):
        """Verify strip_idat_suffix removes all IDAT suffix variants."""
        from illumeta import strip_idat_suffix
        self.assertEqual(strip_idat_suffix("/path/sample_Grn.idat"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample_Red.idat"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample_Grn.idat.gz"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample_Red.idat.gz"), "/path/sample")
        self.assertEqual(strip_idat_suffix("/path/sample"), "/path/sample")
        self.assertEqual(strip_idat_suffix(""), "")


if __name__ == "__main__":
    unittest.main()
