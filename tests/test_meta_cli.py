import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


BASE_DIR = Path(__file__).resolve().parents[1]
ILLUMETA = BASE_DIR / "illumeta.py"


def write_dmp_table(path: Path, effect_scale: float = 1.0) -> None:
    rows = [
        {
            "CpG": "cg_good",
            "Gene": "GENE1",
            "chr": "chr1",
            "pos": "100",
            "Region": "Body",
            "Island_Context": "Island",
            "logFC": str(0.20 * effect_scale),
            "t": "10",
            "P.Value": "1e-8",
            "adj.P.Val": "1e-5",
            "Delta_Beta": str(0.025 * effect_scale),
        },
        {
            "CpG": "cg_noise",
            "Gene": "GENE2",
            "chr": "chr2",
            "pos": "200",
            "Region": "TSS1500",
            "Island_Context": "OpenSea",
            "logFC": "0.01",
            "t": "0.1",
            "P.Value": "0.9",
            "adj.P.Val": "0.95",
            "Delta_Beta": "0.001",
        },
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def rewrite_p_values(path: Path, value: str) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["P.Value"] = value
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_result_dir(root: Path, name: str, effect_scale: float = 1.0) -> Path:
    result_dir = root / name / "Case_vs_Control_results"
    result_dir.mkdir(parents=True)
    (result_dir / "summary.json").write_text(
        json.dumps(
            {
                "n_con": 5,
                "n_test": 5,
                "primary_result_mode": "standard",
                "primary_lambda_guard_status": "ok",
            }
        ),
        encoding="utf-8",
    )
    write_dmp_table(result_dir / "Minfi_DMPs_full.csv", effect_scale)
    write_dmp_table(result_dir / "Sesame_DMPs_full.csv", effect_scale * 0.9)
    return result_dir


def make_result_dir_without_branch(root: Path, name: str) -> Path:
    result_dir = root / name / "Case_vs_Control_results"
    result_dir.mkdir(parents=True)
    (result_dir / "summary.json").write_text(
        json.dumps({"n_con": 5, "n_test": 5, "primary_result_mode": "standard"}),
        encoding="utf-8",
    )
    return result_dir


def gzip_branch_table(path: Path) -> Path:
    gz_path = path.with_suffix(path.suffix + ".gz")
    with path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
        dst.write(src.read())
    path.unlink()
    return gz_path


class MetaCliTests(unittest.TestCase):
    def run_illumeta(self, *args, cwd=None):
        env = os.environ.copy()
        env["ILLUMETA_ALLOW_NON_CONDA"] = "1"
        return subprocess.run(
            [sys.executable, str(ILLUMETA), *args],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=cwd or BASE_DIR,
            env=env,
        )

    def test_meta_cli_runs_branch_separated_cross_cohort_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
                make_result_dir(root, "GSE_D", 1.05),
            ]
            out_dir = root / "meta_out"

            result = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi,sesame_strict",
                "--min-cohorts",
                "3",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            summary = (out_dir / "meta_branch_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("minfi", summary)
            self.assertIn("sesame_strict", summary)
            concordant = (out_dir / "branch_concordant_core_candidates.tsv").read_text(encoding="utf-8")
            self.assertIn("cg_good", concordant)
            self.assertIn("minfi;sesame_strict", concordant)
            report = (out_dir / "meta_analysis_report.md").read_text(encoding="utf-8")
            self.assertIn("Branches are analyzed separately", report)

    def test_meta_cli_requires_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.run_illumeta("meta", "--output", str(Path(tmp) / "out"))

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Provide result directories or --manifest", result.stderr)

    def test_meta_cli_accepts_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
                make_result_dir(root, "GSE_D", 1.05),
            ]
            manifest = root / "manifest.tsv"
            with manifest.open("w", encoding="utf-8") as handle:
                handle.write("cohort\tresult_dir\tplatform\ttissue\n")
                for idx, result_dir in enumerate(result_dirs, start=1):
                    handle.write(f"C{idx}\t{result_dir}\t450K\tbrain\n")
            out_dir = root / "meta_manifest_out"

            result = self.run_illumeta(
                "meta",
                "--manifest",
                str(manifest),
                "--branches",
                "minfi",
                "--min-cohorts",
                "3",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual([c["cohort_id"] for c in payload["cohorts"]], ["C1", "C2", "C3", "C4"])

    def test_meta_cli_accepts_gzipped_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            manifest = root / "manifest.tsv.gz"
            with gzip.open(manifest, "wt", encoding="utf-8", newline="") as handle:
                handle.write("cohort\tresult_dir\tplatform\ttissue\n")
                for idx, result_dir in enumerate(result_dirs, start=1):
                    handle.write(f"C{idx}\t{result_dir}\t450K\tbrain\n")
            out_dir = root / "meta_manifest_gz_out"

            result = self.run_illumeta(
                "meta",
                "--manifest",
                str(manifest),
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual([c["cohort_id"] for c in payload["cohorts"]], ["C1", "C2", "C3"])

    def test_meta_cli_keeps_duplicate_cohort_columns_distinct(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            manifest = root / "manifest.tsv"
            with manifest.open("w", encoding="utf-8") as handle:
                handle.write("cohort\tresult_dir\n")
                for result_dir in result_dirs:
                    handle.write(f"DUP\t{result_dir}\n")
            out_dir = root / "meta_duplicate_ids"

            result = self.run_illumeta(
                "meta",
                "--manifest",
                str(manifest),
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(out_dir / "minfi_meta_full.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(len(reader.fieldnames or []), len(set(reader.fieldnames or [])))
                self.assertIn("effect_DUP", reader.fieldnames)
                self.assertIn("effect_DUP_2", reader.fieldnames)
                self.assertIn("effect_DUP_3", reader.fieldnames)
                row = next(row for row in reader if row["CpG"] == "cg_good")
            self.assertEqual(row["effect_DUP"], "0.2")
            self.assertEqual(row["effect_DUP_2"], "0.22")
            self.assertEqual(row["effect_DUP_3"], "0.19")

    def test_meta_cli_reads_gzipped_branch_file_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            for result_dir in result_dirs:
                gzip_branch_table(result_dir / "Minfi_DMPs_full.csv")
            out_dir = root / "meta_gzip_input"

            result = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--branch-file",
                "minfi=Minfi_DMPs_full.csv.gz",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            concordant = (out_dir / "branch_concordant_core_candidates.tsv").read_text(encoding="utf-8")
            self.assertIn("cg_good", concordant)

    def test_meta_cli_reused_output_removes_stale_owned_branch_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            out_dir = root / "meta_reused"

            first = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi,sesame_strict",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--top-n",
                "1",
                "--output",
                str(out_dir),
            )
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            self.assertTrue((out_dir / "sesame_strict_meta_full.tsv.gz").exists())
            self.assertTrue((out_dir / "minfi_meta_top1.tsv").exists())

            second = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--top-n",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            self.assertFalse((out_dir / "sesame_strict_meta_full.tsv.gz").exists())
            self.assertFalse((out_dir / "sesame_strict_meta_top1.tsv").exists())
            self.assertFalse((out_dir / "sesame_strict_core_candidates.tsv").exists())
            self.assertFalse((out_dir / "sesame_strict_pc_candidates.tsv").exists())
            self.assertFalse((out_dir / "minfi_meta_top1.tsv").exists())
            self.assertTrue((out_dir / "minfi_meta_top2.tsv").exists())
            self.assertFalse(any(path.name.startswith(".illumeta_meta_tmp_") for path in out_dir.iterdir()))
            payload = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["branches"], {"minfi": "Minfi_DMPs_full.csv"})

    def test_meta_cli_failed_reuse_keeps_previous_success_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            out_dir = root / "meta_failed_reuse"

            first = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--top-n",
                "1",
                "--output",
                str(out_dir),
            )
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            before = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))

            second = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi,sesame_native",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--top-n",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertNotEqual(second.returncode, 0)
            after = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(after["branches"], before["branches"])
            self.assertTrue((out_dir / "minfi_meta_top1.tsv").exists())
            self.assertFalse((out_dir / "minfi_meta_top2.tsv").exists())
            self.assertFalse(any(path.name.startswith(".illumeta_meta_tmp_") for path in out_dir.iterdir()))
            failure = json.loads((out_dir / "failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(failure["code"], "META_ANALYSIS_FAILED")

    def test_meta_cli_success_removes_stale_failure_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            out_dir = root / "meta_stale_failure"
            out_dir.mkdir()
            (out_dir / "failure_summary.json").write_text('{"code":"OLD"}\n', encoding="utf-8")
            (out_dir / "failure_reason.txt").write_text("OLD\n", encoding="utf-8")

            result = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse((out_dir / "failure_summary.json").exists())
            self.assertFalse((out_dir / "failure_reason.txt").exists())
            self.assertTrue((out_dir / "meta_input_manifest.json").exists())

    def test_meta_publish_rolls_back_when_file_replace_fails(self):
        import illumeta_meta

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_dir = root / "out"
            stage_dir = root / "stage"
            out_dir.mkdir()
            stage_dir.mkdir()
            (out_dir / "meta_input_manifest.json").write_text('{"branches":{"minfi":"old"}}\n', encoding="utf-8")
            (out_dir / "minfi_meta_top1.tsv").write_text("old\n", encoding="utf-8")
            (stage_dir / "meta_input_manifest.json").write_text('{"branches":{"minfi":"new"}}\n', encoding="utf-8")
            (stage_dir / "minfi_meta_top2.tsv").write_text("new\n", encoding="utf-8")
            real_replace = illumeta_meta.os.replace

            def flaky_replace(src, dst):
                if Path(src).name == "minfi_meta_top2.tsv":
                    raise OSError("simulated publish failure")
                return real_replace(src, dst)

            with mock.patch.object(illumeta_meta.os, "replace", side_effect=flaky_replace):
                with self.assertRaises(OSError):
                    illumeta_meta._publish_meta_outputs(stage_dir, out_dir, ["minfi"])

            manifest = (out_dir / "meta_input_manifest.json").read_text(encoding="utf-8")
            self.assertIn('"old"', manifest)
            self.assertEqual((out_dir / "minfi_meta_top1.tsv").read_text(encoding="utf-8"), "old\n")
            self.assertFalse((out_dir / "minfi_meta_top2.tsv").exists())
            self.assertFalse(any(path.name.startswith(".illumeta_meta_backup_") for path in out_dir.iterdir()))

    def test_meta_cli_rejects_invalid_thresholds(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.run_illumeta(
                "meta",
                "--top-n",
                "0",
                "--output",
                str(Path(tmp) / "out"),
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--top-n must be at least 1", result.stderr)
            payload = json.loads((Path(tmp) / "out" / "failure_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["code"], "META_ANALYSIS_FAILED")
            self.assertEqual(payload["stage"], "meta")

    def test_meta_cli_ignores_out_of_range_input_p_values_for_pc_screen(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            for result_dir in result_dirs:
                rewrite_p_values(result_dir / "Minfi_DMPs_full.csv", "-0.01")
            out_dir = root / "meta_bad_p"

            result = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(out_dir / "minfi_meta_full.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                row = next(row for row in reader if row["CpG"] == "cg_good")
            self.assertEqual(row["pc_directional_p"], "")
            self.assertEqual(row["pc_directional_fdr"], "")
            self.assertEqual(row["replicable_pc_candidate"], "false")

    def test_meta_cli_warns_when_summary_json_is_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir(root, "GSE_A", 1.0),
                make_result_dir(root, "GSE_B", 1.1),
                make_result_dir(root, "GSE_C", 0.95),
            ]
            (result_dirs[0] / "summary.json").unlink()
            out_dir = root / "meta_missing_summary"

            result = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(any("summary.json missing" in warning for warning in payload["warnings"]))
            report = (out_dir / "meta_analysis_report.md").read_text(encoding="utf-8")
            self.assertIn("summary.json missing", report)

    def test_meta_cli_warns_when_allowed_branch_is_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dirs = [
                make_result_dir_without_branch(root, "GSE_A"),
                make_result_dir_without_branch(root, "GSE_B"),
                make_result_dir_without_branch(root, "GSE_C"),
            ]
            out_dir = root / "meta_empty_branch"

            result = self.run_illumeta(
                "meta",
                *[str(path) for path in result_dirs],
                "--branches",
                "minfi",
                "--allow-missing-branches",
                "--min-cohorts",
                "2",
                "--partial-conjunction-r",
                "2",
                "--output",
                str(out_dir),
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads((out_dir / "meta_input_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(any("branch output is empty" in warning for warning in payload["warnings"]))


if __name__ == "__main__":
    unittest.main()
