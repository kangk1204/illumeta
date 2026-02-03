#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/environment.yml"
ENV_NAME=""
RUN_DOCTOR=1

usage() {
  cat <<'USAGE'
Usage: scripts/install_full.sh [options]

Options:
  --r45              Use environment-r45.yml (R 4.5)
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

LOG_DIR="${ROOT_DIR}/projects"
mkdir -p "${LOG_DIR}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/illumeta_install_full_${RUN_ID}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[*] IlluMeta full install (conda)"
echo "[*] Env file: ${ENV_FILE}"
echo "[*] Env name: ${ENV_NAME}"
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

export ILLUMETA_USE_CONDA_LIBS="${ILLUMETA_USE_CONDA_LIBS:-1}"
export ILLUMETA_CLEAN_MISMATCHED_RLIB="${ILLUMETA_CLEAN_MISMATCHED_RLIB:-1}"
export ILLUMETA_DOWNLOAD_RETRIES="${ILLUMETA_DOWNLOAD_RETRIES:-3}"
export ILLUMETA_FORCE_SETUP="${ILLUMETA_FORCE_SETUP:-1}"
export ILLUMETA_INSTALL_DEVTOOLS="${ILLUMETA_INSTALL_DEVTOOLS:-1}"
export ILLUMETA_INSTALL_CLOCKS="${ILLUMETA_INSTALL_CLOCKS:-1}"
export ILLUMETA_REQUIRE_EPICV2="${ILLUMETA_REQUIRE_EPICV2:-1}"

echo "[*] Running R setup_env.R (full install)..."
"${CONDA_BIN}" run -n "${ENV_NAME}" Rscript "${ROOT_DIR}/r_scripts/setup_env.R"

if [[ "${RUN_DOCTOR}" -eq 1 ]]; then
  echo "[*] Running illumeta doctor..."
  "${CONDA_BIN}" run -n "${ENV_NAME}" python "${ROOT_DIR}/illumeta.py" doctor
fi

echo "[*] Full install completed."
