import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


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


if __name__ == "__main__":
    unittest.main()
