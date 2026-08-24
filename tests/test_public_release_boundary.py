"""Fail closed when private research or generated deliverables enter the public tree."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).relative_to(ROOT).as_posix()
MAX_TRACKED_SIZE = 10 * 1024 * 1024
FORBIDDEN_ROOTS = {
    "archive",
    "benchmarks",
    "manuscript",
    "outputs",
    "paper_assets",
    "references",
    "reports",
    "submission",
}
FORBIDDEN_COMPONENTS = {".claude", ".codex", ".omc", ".omx"}
FORBIDDEN_PREFIXES = {
    ("docs", "submission_package"),
}
FORBIDDEN_NAME_TOKENS = (
    "ai_use",
    "agent_replay",
    "chatgpt",
    "claude",
    "codex",
    "openai",
    "subagent",
)
FORBIDDEN_BINARY_SUFFIXES = {".docx", ".idat", ".pdf", ".rds", ".tar", ".tgz", ".xlsx", ".zip"}
FORBIDDEN_CONTENT = (
    "dementia_special_issue",
    "illumeta_manuscript_package",
    "requirements-paper.txt",
    "prepare_geo_submission.py",
    "readme_dashboard_gse125605.png",
    "cg11979621",
    "RBM33",
    "AI use log",
    "agent-replay",
    "suggested reviewers",
)
SECRET_PATTERNS = (
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"gh[pousr]_[A-Za-z0-9_]{36,}"),
    re.compile(r"xox[baprs]-[A-Za-z0-9-]{20,}"),
    re.compile(r"sk-[A-Za-z0-9_-]{32,}"),
    re.compile(r"BEGIN [A-Z ]*PRIVATE KEY"),
    re.compile(r"(?i)(?:api[_-]?key|password|secret|token)\s*[:=]\s*['\"][^'\"]{16,}"),
)


def _tracked_paths() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    paths = [Path(raw.decode()) for raw in result.stdout.split(b"\0") if raw]
    return [path for path in paths if (ROOT / path).is_file()]


def test_tracked_paths_respect_public_release_boundary():
    problems = []
    for path in _tracked_paths():
        parts = path.parts
        lowered_parts = tuple(part.lower() for part in parts)
        name = path.name.lower()
        if parts and parts[0].lower() in FORBIDDEN_ROOTS:
            problems.append(f"forbidden root: {path}")
        if any(part in FORBIDDEN_COMPONENTS for part in lowered_parts):
            problems.append(f"forbidden workspace state: {path}")
        if any(lowered_parts[: len(prefix)] == prefix for prefix in FORBIDDEN_PREFIXES):
            problems.append(f"forbidden generated deliverable: {path}")
        if any(token in name for token in FORBIDDEN_NAME_TOKENS):
            problems.append(f"forbidden filename token: {path}")
        if path.suffix.lower() in FORBIDDEN_BINARY_SUFFIXES:
            problems.append(f"forbidden binary deliverable: {path}")
        absolute = ROOT / path
        if absolute.is_file() and absolute.stat().st_size > MAX_TRACKED_SIZE:
            problems.append(f"tracked file exceeds 10 MiB: {path}")
    assert problems == []


def test_tracked_text_has_no_private_identifiers_or_secrets():
    problems = []
    for path in _tracked_paths():
        relative = path.as_posix()
        absolute = ROOT / path
        if relative == SELF or not absolute.is_file() or absolute.stat().st_size > 2_000_000:
            continue
        try:
            text = absolute.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for marker in FORBIDDEN_CONTENT:
            if marker in text:
                problems.append(f"{path}: forbidden marker {marker!r}")
        for pattern in SECRET_PATTERNS:
            if pattern.search(text):
                problems.append(f"{path}: possible secret matching {pattern.pattern!r}")
    assert problems == []


def test_public_readme_has_no_private_workflow_language():
    readme = (ROOT / "README.md").read_text(encoding="utf-8").lower()
    for marker in (
        "application note",
        "cover letter",
        "manuscript",
        "requirements-paper",
        "submission helper",
        "supplementary materials",
        "--require-publication-artifacts",
    ):
        assert marker not in readme
