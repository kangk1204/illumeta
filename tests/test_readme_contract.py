"""README checks for the beginner path and documented public interface."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"


def test_fresh_ubuntu_block_is_copy_paste_complete():
    text = README.read_text(encoding="utf-8")
    for marker in (
        "sudo apt-get install -y ca-certificates curl git",
        'case "$(uname -m)" in',
        "Miniforge3-Linux-${MINIFORGE_ARCH}.sh",
        "git clone https://github.com/kangk1204/illumeta.git",
        "./scripts/install_full.sh --preflight",
        "./scripts/install_full.sh",
        "./scripts/illumeta doctor",
        "./scripts/illumeta demo",
    ):
        assert marker in text


def test_normal_cli_examples_use_the_environment_wrapper():
    text = README.read_text(encoding="utf-8")
    direct_commands = re.findall(
        r"(?m)^\s*(?:[A-Z][A-Z0-9_]*=[^\s]+\s+)*python illumeta\.py\b.*$",
        text,
    )
    assert direct_commands == []
    assert "use `./scripts/illumeta ...`" in text


def test_documented_core_commands_exist_in_cli_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "illumeta.py"), "--help"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for command in ("download", "search", "analysis", "meta", "demo", "doctor"):
        assert re.search(rf"\b{command}\b", result.stdout)


def test_dashboard_capabilities_and_output_contract_are_documented():
    text = README.read_text(encoding="utf-8")
    for marker in (
        "search, sort, and export CSV",
        "analysis_parameters.json",
        "QC_Summary.csv",
        "methods.md",
        "sessionInfo.txt",
        "Case_vs_Control_results_index.html",
        "docs/dashboard_demo/Case_vs_Control_results_index.html",
    ):
        assert marker in text
