"""Integrity tests for the deterministic synthetic dashboard artifact."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_dashboard_demo.py"
DEMO_DIR = ROOT / "docs" / "dashboard_demo"
RESULTS_DIR = DEMO_DIR / "Case_vs_Control_results"
DEMO = DEMO_DIR / "Case_vs_Control_results_index.html"
MANIFEST = DEMO_DIR / "manifest.json"
PLOT = RESULTS_DIR / "Minfi_Volcano.html"
TABLE = RESULTS_DIR / "Minfi_Top_DMPs.html"
CONSENSUS_TABLE = RESULTS_DIR / "Intersection_Native_Consensus_DMPs.html"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_dashboard_demo_is_deterministic_offline_and_disclosed():
    subprocess.run([sys.executable, str(BUILDER)], cwd=ROOT, check=True)
    first = _sha256(DEMO)
    subprocess.run([sys.executable, str(BUILDER)], cwd=ROOT, check=True)
    second = _sha256(DEMO)

    assert first == second
    html = DEMO.read_text(encoding="utf-8")
    assert "Synthetic interface fixture." in html
    assert "not study results" in html
    assert 'const RESULTS_PATH = "Case_vs_Control_results";' in html
    assert str(ROOT) not in html
    assert "/tmp/" not in html
    assert re.search(r"https?://", html) is None
    assert 'role="tablist"' in html
    assert ":focus-visible" in html
    assert '<link rel="icon" href="data:,"' in html
    assert "Cross-reactive probe filtering was skipped" not in html
    assert "N/A" not in html
    assert "0.05 / 0.5 / 0.0" in html
    assert "12 / 0" in html
    assert "8830" in html
    assert "Minfi_Volcano.html" in html
    assert "Minfi_Top_DMPs.html" in html
    assert "Intersection_Native_Consensus_DMPs.html" in html

    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert payload["fixture_only"] is True
    assert payload["study_result"] is False
    assert payload["artifact_sha256"] == second
    declared = set(payload["inputs"])
    for relative, expected_sha256 in payload["inputs"].items():
        assert _sha256(ROOT / relative) == expected_sha256
    for path in (
        PLOT,
        TABLE,
        CONSENSUS_TABLE,
        RESULTS_DIR / "Intersection_Native_Consensus_DMPs.csv",
        RESULTS_DIR / "Minfi_Metrics.csv",
        RESULTS_DIR / "CRF_MMC_Summary.csv",
        RESULTS_DIR / "CRF_RSS_Summary.csv",
        RESULTS_DIR / "CRF_NCS_Summary.csv",
        RESULTS_DIR / "CRF_Sample_Tier.csv",
    ):
        assert path.relative_to(ROOT).as_posix() in declared

    parameters = json.loads(
        (RESULTS_DIR / "analysis_parameters.json").read_text(encoding="utf-8")
    )
    assert parameters["pval_threshold"] == 0.05
    assert parameters["lfc_threshold"] == 0.5
    assert parameters["cross_reactive_active"] is True
    assert parameters["cross_reactive_count"] > 0

    qc_rows = (RESULTS_DIR / "QC_Summary.csv").read_text(encoding="utf-8")
    assert "Samples_passed_QC,12" in qc_rows
    assert "Probes_final,8830" in qc_rows

    for child in (PLOT, TABLE, CONSENSUS_TABLE):
        child_html = child.read_text(encoding="utf-8")
        assert re.search(r"""(?:src|href)=["']https?://""", child_html) is None
        assert "Synthetic interface fixture" in child_html
    assert 'data-illumeta-fixture="interactive-plot"' in PLOT.read_text(
        encoding="utf-8"
    )
    table_html = TABLE.read_text(encoding="utf-8")
    assert 'data-illumeta-fixture="interactive-table"' in table_html
    assert "sortRows" in table_html
    assert "downloadCsv" in table_html
