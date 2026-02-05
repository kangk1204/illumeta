import csv
import importlib.util
import json
import os
import tempfile
import unittest


def load_benchmark_module():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    module_path = os.path.join(repo_root, "scripts", "build_benchmark_table.py")
    spec = importlib.util.spec_from_file_location("build_benchmark_table", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path, fieldnames, rows, delimiter=","):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


class BenchmarkTableTests(unittest.TestCase):
    def test_build_benchmark_table_fields(self):
        bbt = load_benchmark_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = os.path.join(tmpdir, "GSE_TEST")
            os.makedirs(result_dir, exist_ok=True)

            summary = {
                "n_con": 5,
                "n_test": 5,
                "primary_result_mode": "tier3_stratified_meta",
                "primary_tier3_meta_method": "random",
                "primary_tier3_meta_i2_median": 0.6,
                "primary_dmr_status": "ran_tier3_meta",
                "primary_branch_override": "Minfi",
                "primary_branch_reason": "sesame_fallback:noob_only",
                "sesame_dyebias_mode": "noob_only",
                "sesame_dyebias_note": "missing_controls",
            }
            with open(os.path.join(result_dir, "summary.json"), "w", encoding="utf-8") as handle:
                json.dump(summary, handle)

            with open(os.path.join(result_dir, "analysis_parameters.json"), "w", encoding="utf-8") as handle:
                json.dump({"array_type": "EPIC", "preset": "conservative"}, handle)

            write_csv(
                os.path.join(result_dir, "QC_Summary.csv"),
                ["metric", "value"],
                [
                    {"metric": "Total_samples_input", "value": "10"},
                    {"metric": "Samples_passed_QC", "value": "8"},
                    {"metric": "Total_probes_raw", "value": "100"},
                    {"metric": "Probes_final", "value": "50"},
                ],
            )

            write_csv(
                os.path.join(result_dir, "Minfi_Metrics.csv"),
                ["metric", "value"],
                [
                    {"metric": "n_samples", "value": "10"},
                    {"metric": "n_cpgs", "value": "100"},
                    {"metric": "batch_sig_p_lt_0.05_before", "value": "10"},
                    {"metric": "batch_sig_p_lt_0.05_after", "value": "4"},
                    {"metric": "batch_min_p_before", "value": "0.001"},
                    {"metric": "batch_min_p_after", "value": "0.2"},
                    {"metric": "group_min_p_before", "value": "0.05"},
                    {"metric": "group_min_p_after", "value": "0.001"},
                ],
            )

            write_csv(
                os.path.join(result_dir, "Input_Group_Distribution.csv"),
                ["group", "n"],
                [{"group": "Control", "n": "5"}, {"group": "Case", "n": "5"}],
            )

            write_csv(
                os.path.join(result_dir, "decision_ledger.tsv"),
                ["stage", "decision", "value"],
                [
                    {"stage": "consensus", "decision": "primary_branch", "value": "Minfi"},
                    {"stage": "confounding", "decision": "tier3_batch", "value": "Sentrix_ID"},
                ],
                delimiter="\t",
            )

            rows = [{"gse_id": "GSE_TEST", "analysis_dir": result_dir}]
            out_rows, missing = bbt.build_rows(rows, projects_root="")

            self.assertEqual(len(missing), 0)
            self.assertEqual(len(out_rows), 1)
            row = out_rows[0]
            self.assertEqual(row["qc_pass_rate"], "0.8")
            self.assertEqual(row["qc_probe_retention"], "0.5")
            self.assertEqual(row["batch_sig_reduction"], "6.0")
            self.assertEqual(row["primary_tier3_meta_method"], "random")
            self.assertEqual(str(row["primary_tier3_meta_i2_median"]), "0.6")


if __name__ == "__main__":
    unittest.main()
