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


if __name__ == "__main__":
    unittest.main()
