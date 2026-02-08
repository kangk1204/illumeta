#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/environment.yml"
ENV_NAME=""
RUN_DOCTOR=1
INSTALL_PAPER=0
INSTALL_MINIMAL=0
INSTALL_CLOCKS=0
INSTALL_DEVTOOLS=0
INSTALL_EPICV2=0

usage() {
  cat <<'USAGE'
Usage: scripts/install_full.sh [options]

Options:
  --r45              Use environment-r45.yml (R 4.5)
  --paper            Also install Python deps for paper/figure generation (requirements-paper.txt)
  --minimal          Core-only R install (skip optional cell refs/RefFreeEWAS/planet; some features disabled)
  --clocks           Install optional epigenetic clock packages (methylclock/wateRmelon)
  --devtools         Install optional devtools/tidyverse (for development)
  --epicv2           Install optional EPIC v2 manifest/annotation packages
  --full             Install all optional features: --epicv2 --clocks --devtools
  --env-file PATH    Use a custom conda env file
  --env NAME         Override env name (default: name: from env file)
  --skip-doctor      Skip `illumeta.py doctor` after install
  -h, --help         Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --r45)
      ENV_FILE="${ROOT_DIR}/environment-r45.yml"
      shift
      ;;
    --env-file)
      ENV_FILE="$2"
      shift 2
      ;;
    --env)
      ENV_NAME="$2"
      shift 2
      ;;
    --paper)
      INSTALL_PAPER=1
      shift
      ;;
    --minimal)
      INSTALL_MINIMAL=1
      shift
      ;;
    --clocks)
      INSTALL_CLOCKS=1
      shift
      ;;
    --devtools)
      INSTALL_DEVTOOLS=1
      shift
      ;;
    --epicv2)
      INSTALL_EPICV2=1
      shift
      ;;
    --full)
      INSTALL_CLOCKS=1
      INSTALL_DEVTOOLS=1
      INSTALL_EPICV2=1
      shift
      ;;
    --skip-doctor)
      RUN_DOCTOR=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[!] Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "[!] Missing env file: ${ENV_FILE}"
  exit 1
fi

if [[ -z "${ENV_NAME}" ]]; then
  ENV_NAME="$(awk -F: '/^name:/ {gsub(/ /,"",$2); print $2; exit}' "${ENV_FILE}")"
fi

CONDA_BIN=""
if command -v conda >/dev/null 2>&1; then
  CONDA_BIN="conda"
elif command -v mamba >/dev/null 2>&1; then
  CONDA_BIN="mamba"
else
  echo "[!] conda/mamba not found in PATH."
  echo "    Install Miniforge first: https://github.com/conda-forge/miniforge"
  exit 1
fi

# macOS needs Xcode Command Line Tools even when using conda, because compilers
# rely on the system SDK headers. Fail early with a clear instruction.
if [[ "$(uname -s)" == "Darwin" ]]; then
  if ! xcode-select -p >/dev/null 2>&1; then
    echo "[!] Xcode Command Line Tools not found."
    echo "    Install them first: xcode-select --install"
    exit 1
  fi
fi

LOG_DIR="${ROOT_DIR}/projects"
mkdir -p "${LOG_DIR}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/illumeta_install_full_${RUN_ID}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[*] IlluMeta full install (conda)"
echo "[*] Env file: ${ENV_FILE}"
echo "[*] Env name: ${ENV_NAME}"
echo "[*] R options: minimal=${INSTALL_MINIMAL} epicv2=${INSTALL_EPICV2} clocks=${INSTALL_CLOCKS} devtools=${INSTALL_DEVTOOLS}"
echo "[*] Log file: ${LOG_FILE}"

ENV_EXISTS=0
if "${CONDA_BIN}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  ENV_EXISTS=1
fi

if [[ "${ENV_EXISTS}" -eq 1 ]]; then
  echo "[*] Updating conda env..."
  "${CONDA_BIN}" env update -f "${ENV_FILE}" --prune
else
  echo "[*] Creating conda env..."
  if ! "${CONDA_BIN}" env create -f "${ENV_FILE}"; then
    if [[ "${CONDA_BIN}" == "conda" ]]; then
      echo "[!] env create failed; retrying with --solver=classic"
      "${CONDA_BIN}" env create -f "${ENV_FILE}" --solver=classic
    else
      exit 1
    fi
  fi
fi

echo "[*] Ensuring Python deps..."
"${CONDA_BIN}" run -n "${ENV_NAME}" python -m pip install --upgrade pip
"${CONDA_BIN}" run -n "${ENV_NAME}" python -m pip install -r "${ROOT_DIR}/requirements.txt"
if [[ "${INSTALL_PAPER}" -eq 1 ]]; then
  echo "[*] Installing paper/figure Python deps (requirements-paper.txt)..."
  "${CONDA_BIN}" run -n "${ENV_NAME}" python -m pip install -r "${ROOT_DIR}/requirements-paper.txt"
fi

export ILLUMETA_USE_CONDA_LIBS="${ILLUMETA_USE_CONDA_LIBS:-1}"
export ILLUMETA_CLEAN_MISMATCHED_RLIB="${ILLUMETA_CLEAN_MISMATCHED_RLIB:-1}"
export ILLUMETA_DOWNLOAD_RETRIES="${ILLUMETA_DOWNLOAD_RETRIES:-3}"
export ILLUMETA_FORCE_SETUP="${ILLUMETA_FORCE_SETUP:-1}"
export ILLUMETA_INSTALL_MINIMAL="${ILLUMETA_INSTALL_MINIMAL:-${INSTALL_MINIMAL}}"
export ILLUMETA_INSTALL_DEVTOOLS="${ILLUMETA_INSTALL_DEVTOOLS:-${INSTALL_DEVTOOLS}}"
export ILLUMETA_INSTALL_CLOCKS="${ILLUMETA_INSTALL_CLOCKS:-${INSTALL_CLOCKS}}"
export ILLUMETA_REQUIRE_EPICV2="${ILLUMETA_REQUIRE_EPICV2:-${INSTALL_EPICV2}}"

echo "[*] Running R setup_env.R..."
"${CONDA_BIN}" run -n "${ENV_NAME}" Rscript "${ROOT_DIR}/r_scripts/setup_env.R"

if [[ "${RUN_DOCTOR}" -eq 1 ]]; then
  echo "[*] Running illumeta doctor..."
  "${CONDA_BIN}" run -n "${ENV_NAME}" python "${ROOT_DIR}/illumeta.py" doctor
fi

echo "[*] Full install completed."
