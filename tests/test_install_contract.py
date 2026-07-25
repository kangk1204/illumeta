"""Contract tests for the beginner installer entrypoint."""

from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts" / "install_full.sh"
DEFAULT_ENV = ROOT / "environment.yml"


def _fake_home(tmp_path: Path) -> Path:
    home = tmp_path / "home"
    conda = home / "miniforge3" / "bin" / "conda"
    conda.parent.mkdir(parents=True)
    conda.write_text("#!/usr/bin/env sh\nexit 0\n", encoding="utf-8")
    conda.chmod(conda.stat().st_mode | stat.S_IXUSR)
    return home


def _retrying_fake_home(tmp_path: Path) -> tuple[Path, Path]:
    home = tmp_path / "home"
    conda = home / "miniforge3" / "bin" / "conda"
    conda.parent.mkdir(parents=True)
    call_log = tmp_path / "conda-calls.log"
    conda.write_text(
        """#!/usr/bin/env sh
set -eu
printf '%s\n' "$*" >> "$FAKE_CONDA_CALL_LOG"
if [ "${1:-} ${2:-}" = "env list" ]; then
    if [ -f "$FAKE_CONDA_CREATED" ]; then
        printf '%s %s\n' "$FAKE_CONDA_ENV_NAME" "$HOME/miniforge3/envs/$FAKE_CONDA_ENV_NAME"
    fi
    exit 0
fi
if [ "${1:-} ${2:-}" = "env create" ]; then
    if [ ! -f "$FAKE_CONDA_FAILED_ONCE" ]; then
        : > "$FAKE_CONDA_FAILED_ONCE"
        exit 42
    fi
    : > "$FAKE_CONDA_CREATED"
fi
exit 0
""",
        encoding="utf-8",
    )
    conda.chmod(conda.stat().st_mode | stat.S_IXUSR)
    return home, call_log


def _preflight(tmp_path: Path, *args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    home = _fake_home(tmp_path)
    env["HOME"] = str(home)
    env["PATH"] = "/usr/bin:/bin"
    return subprocess.run(
        [str(INSTALLER), "--preflight", *args],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def test_default_preflight_uses_supported_environment(tmp_path):
    result = _preflight(tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert str(tmp_path / "home" / "miniforge3" / "bin" / "conda") in result.stdout
    assert f"Env file: {DEFAULT_ENV}" in result.stdout
    assert "Env name: illumeta" in result.stdout


def test_missing_env_file_argument_is_reported(tmp_path):
    result = _preflight(tmp_path, "--env-file")
    assert result.returncode == 2
    assert "--env-file requires a path" in result.stdout


def test_installer_retries_transient_conda_transaction_failure(tmp_path):
    env = os.environ.copy()
    home, call_log = _retrying_fake_home(tmp_path)
    env.update(
        {
            "HOME": str(home),
            "PATH": "/usr/bin:/bin",
            "FAKE_CONDA_CALL_LOG": str(call_log),
            "FAKE_CONDA_FAILED_ONCE": str(tmp_path / "failed-once"),
            "FAKE_CONDA_CREATED": str(tmp_path / "created"),
            "FAKE_CONDA_ENV_NAME": "illumeta-retry-test",
            "ILLUMETA_CONDA_RETRY_DELAY": "0",
            "ILLUMETA_LOG_DIR": str(tmp_path / "install-logs"),
        }
    )
    result = subprocess.run(
        [
            str(INSTALLER),
            "--env",
            "illumeta-retry-test",
            "--minimal",
            "--skip-doctor",
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    calls = call_log.read_text(encoding="utf-8").splitlines()
    assert sum(call.startswith("env create ") for call in calls) == 2
    assert "failed (attempt 1/3)" in result.stdout
    assert len(list((tmp_path / "install-logs").glob("illumeta_install_full_*.log"))) == 1
    assert "[*] Full install completed." in result.stdout
