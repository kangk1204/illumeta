import csv
import os
import tempfile
import unittest

from illumeta import auto_group_config


def write_config(path, headers, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


def read_primary_groups(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row.get("primary_group", "") for row in reader]


class AutoGroupTests(unittest.TestCase):
    def test_auto_group_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["primary_group", "disease_state", "Basename"]
            rows = [
                {"primary_group": "", "disease_state": "control", "Basename": "S1"},
                {"primary_group": "", "disease_state": "case", "Basename": "S2"},
            ]
            write_config(config_path, headers, rows)

            out_path, info = auto_group_config(
                config_path=config_path,
                group_con="Control",
                group_test="Case",
                group_column="disease_state",
            )

            self.assertTrue(info["updated"])
            self.assertTrue(os.path.exists(out_path))
            self.assertEqual(read_primary_groups(out_path), ["Control", "Case"])

    def test_auto_group_infer_keyword(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["primary_group", "condition", "Basename"]
            rows = [
                {"primary_group": "", "condition": "Control", "Basename": "S1"},
                {"primary_group": "", "condition": "Case", "Basename": "S2"},
            ]
            write_config(config_path, headers, rows)

            out_path, info = auto_group_config(
                config_path=config_path,
                group_con="Control",
                group_test="Case",
            )

            self.assertTrue(info["updated"])
            self.assertTrue(os.path.exists(out_path))
            self.assertEqual(read_primary_groups(out_path), ["Control", "Case"])

    def test_auto_group_requires_signal(self):
        """Auto-group should fail when no valid grouping exists (all same value)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["primary_group", "constant_col", "Basename"]
            rows = [
                {"primary_group": "", "constant_col": "same", "Basename": "S1"},
                {"primary_group": "", "constant_col": "same", "Basename": "S2"},
                {"primary_group": "", "constant_col": "same", "Basename": "S3"},
            ]
            write_config(config_path, headers, rows)

            with self.assertRaises(ValueError):
                auto_group_config(
                    config_path=config_path,
                    group_con="Control",
                    group_test="Case",
                )

    def test_auto_group_rejects_labels_outside_control_test(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "configure.tsv")
            headers = ["primary_group", "condition", "Basename"]
            rows = [
                {"primary_group": "", "condition": "control", "Basename": "S1"},
                {"primary_group": "", "condition": "unknown", "Basename": "S2"},
                {"primary_group": "", "condition": "case", "Basename": "S3"},
            ]
            write_config(config_path, headers, rows)

            with self.assertRaises(ValueError) as cm:
                auto_group_config(
                    config_path=config_path,
                    group_con="Control",
                    group_test="Case",
                    group_column="condition",
                )
            self.assertIn("outside control/test", str(cm.exception))
            self.assertIn("--group-map", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
