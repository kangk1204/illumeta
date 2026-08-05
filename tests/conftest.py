"""Shared pytest configuration.

Turns the R-dependent skips into hard failures when the environment declares that
R must be present.

Why this exists: every test that exercises the R statistical core (probe QC, design
invariants, bacon recalibration, genomic inflation, stratified meta-analysis) calls
`self.skipTest("Rscript not installed")` when Rscript is absent. On a runner without
R that is 16+ silent skips and the suite still reports green, so "all tests passed"
did not imply that any statistical code had been executed. Set ILLUMETA_REQUIRE_R=1
(CI does) and those skips become failures instead.
"""

from __future__ import annotations

import os
import shutil

import pytest

#: Substrings identifying a skip that only happens because R is missing.
_R_SKIP_MARKERS = ("rscript is not installed", "rscript not installed")


def _require_r() -> bool:
    return os.environ.get("ILLUMETA_REQUIRE_R", "").strip().lower() in {"1", "true", "yes", "on"}


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "requires_r: test that shells out to Rscript"
    )
    if _require_r() and shutil.which("Rscript") is None:
        raise pytest.UsageError(
            "ILLUMETA_REQUIRE_R is set but Rscript is not on PATH. The R statistical "
            "tests would silently skip and the suite would report a false green. "
            "Install R (see scripts/install_full.sh) or unset ILLUMETA_REQUIRE_R."
        )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if not _require_r() or not report.skipped:
        return
    reason = ""
    if isinstance(report.longrepr, tuple) and len(report.longrepr) == 3:
        reason = str(report.longrepr[2])
    else:
        reason = str(report.longrepr or "")
    if any(marker in reason.lower() for marker in _R_SKIP_MARKERS):
        report.outcome = "failed"
        report.longrepr = (
            f"{reason}\n\n"
            "ILLUMETA_REQUIRE_R=1 is set, so an R-dependent test may not skip. "
            "This test covers the R statistical core; skipping it would let the suite "
            "report success without executing any statistical code."
        )
