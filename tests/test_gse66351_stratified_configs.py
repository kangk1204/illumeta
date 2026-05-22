"""Tests for GSE66351 stratified config generation guards."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[1]
SCRIPT = BASE_DIR / "scripts" / "build_gse66351_stratified_configs.py"
FIELD_MAP = {
    "cell_type": "characteristics_ch1",
    "diagnosis": "characteristics_ch1.1",
    "braak_stage": "characteristics_ch1.2",
    "brain_region": "characteristics_ch1.3",
    "age": "characteristics_ch1.4",
    "sex": "characteristics_ch1.5",
    "donor_id": "characteristics_ch1.6",
    "Sentrix_ID": "characteristics_ch1.7",
    "Sentrix_Position": "characteristics_ch1.8",
}
HEADERS = ["geo_accession", "source_name_ch1", *FIELD_MAP.values(), "description"]


def make_row(idx: int, cell_type: str, diagnosis: str, brain_region: str) -> dict[str, str]:
    return {
        "geo_accession": f"GSM{idx:07d}",
        "source_name_ch1": brain_region,
        FIELD_MAP["cell_type"]: cell_type,
        FIELD_MAP["diagnosis"]: diagnosis,
        FIELD_MAP["braak_stage"]: "1",
        FIELD_MAP["brain_region"]: brain_region,
        FIELD_MAP["age"]: "70",
        FIELD_MAP["sex"]: "F",
        FIELD_MAP["donor_id"]: f"Dnr{idx:02d}",
        FIELD_MAP["Sentrix_ID"]: "1234567890",
        FIELD_MAP["Sentrix_Position"]: "R01C01",
        "description": "test",
    }


def write_source(path: Path, rows: list[dict[str, str]], headers: list[str] = HEADERS) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows([{key: row.get(key, "") for key in headers} for row in rows])


class GSE66351StratifiedConfigTests(unittest.TestCase):
    def run_script(self, source: Path, out_dir: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--source", str(source), "--out-dir", str(out_dir)],
            cwd=BASE_DIR,
            capture_output=True,
            text=True,
            check=False,
        )

    def valid_rows(self) -> list[dict[str, str]]:
        rows = []
        idx = 1
        for cell_type, brain_region in [
            ("Neuron", "Occipital cortex"),
            ("Glia", "Occipital cortex"),
            ("bulk", "Frontal cortex"),
            ("bulk", "Temporal cortex"),
        ]:
            for diagnosis in ("AD", "CTRL"):
                rows.append(make_row(idx, cell_type, diagnosis, brain_region))
                idx += 1
        return rows

    def test_missing_required_column_fails_without_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source = tmp_path / "configure.tsv"
            out_dir = tmp_path / "out"
            headers = [h for h in HEADERS if h != FIELD_MAP["diagnosis"]]
            write_source(source, self.valid_rows(), headers=headers)

            result = self.run_script(source, out_dir)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("missing required column", result.stderr)
            self.assertFalse((out_dir / "configure_AD_vs_CTRL_neuron_occipital.tsv").exists())

    def test_empty_expected_stratum_fails_without_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source = tmp_path / "configure.tsv"
            out_dir = tmp_path / "out"
            rows = [row for row in self.valid_rows() if row[FIELD_MAP["cell_type"]] != "Neuron"]
            write_source(source, rows)

            result = self.run_script(source, out_dir)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Occipital cortex / Neuron: no rows matched", result.stderr)
            self.assertFalse((out_dir / "configure_AD_vs_CTRL_neuron_occipital.tsv").exists())


if __name__ == "__main__":
    unittest.main()
