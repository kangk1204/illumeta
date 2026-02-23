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
PRECHECK_ONLY=0

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
  --preflight        Validate prerequisites only (no env/package install)
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
    --preflight)
      PRECHECK_ONLY=1
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

miniforge_installer_url() {
  local os arch
  os="$(uname -s)"
  arch="$(uname -m)"
  case "${os}:${arch}" in
    Linux:x86_64) echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" ;;
    Linux:aarch64|Linux:arm64) echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh" ;;
    Darwin:arm64) echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh" ;;
    Darwin:x86_64) echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh" ;;
    *) echo "" ;;
  esac
}

print_conda_install_hint() {
  local url shell_name
  url="$(miniforge_installer_url)"
  shell_name="$(basename "${SHELL:-bash}")"
  if [[ "${shell_name}" != "bash" && "${shell_name}" != "zsh" ]]; then
    shell_name="bash"
  fi

  echo "    Install Miniforge, then rerun this script."
  if [[ -n "${url}" ]]; then
    if [[ "$(uname -s)" == "Linux" ]]; then
      echo "    (Ubuntu/WSL) sudo apt-get update && sudo apt-get install -y curl"
    fi
    echo "    curl -L -o Miniforge3.sh ${url}"
    echo "    bash Miniforge3.sh -b -p \"\$HOME/miniforge3\""
    echo "    \"\$HOME/miniforge3/bin/conda\" init ${shell_name}"
    echo "    source \"\$HOME/.${shell_name}rc\""
  else
    echo "    https://github.com/conda-forge/miniforge"
  fi
}

CONDA_BIN=""
if command -v conda >/dev/null 2>&1; then
  CONDA_BIN="conda"
elif command -v mamba >/dev/null 2>&1; then
  CONDA_BIN="mamba"
else
  for candidate in \
    "${HOME}/miniforge3/bin/conda" \
    "${HOME}/mambaforge/bin/conda" \
    "${HOME}/anaconda3/bin/conda"; do
    if [[ -x "${candidate}" ]]; then
      CONDA_BIN="${candidate}"
      break
    fi
  done
fi

if [[ -z "${CONDA_BIN}" ]]; then
  echo "[!] conda/mamba not found in PATH."
  print_conda_install_hint
  exit 1
fi

if [[ "${CONDA_BIN}" == "${HOME}"/*"/bin/conda" ]]; then
  echo "[*] Using detected conda binary: ${CONDA_BIN}"
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

if [[ "${PRECHECK_ONLY}" -eq 1 ]]; then
  echo "[*] Preflight check passed."
  echo "[*] Env file: ${ENV_FILE}"
  echo "[*] Env name: ${ENV_NAME}"
  echo "[*] Conda command: ${CONDA_BIN}"
  echo "[*] Next: ./scripts/install_full.sh"
  exit 0
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

# macOS: pre-install R packages that frequently fail to compile from source
# (stringi needs ICU, textshaping/ragg/systemfonts need system graphics libs).
# Conda prebuilt binaries bypass these compilation issues entirely.
if [[ "$(uname -s)" == "Darwin" ]]; then
  echo "[*] macOS detected — installing prebuilt R packages via conda to avoid compilation errors..."
  "${CONDA_BIN}" install -n "${ENV_NAME}" -c conda-forge -y \
    r-stringi r-stringr r-tidyr r-rvest r-selectr \
    r-textshaping r-ragg r-systemfonts r-plotly 2>&1 || {
    echo "[!] Warning: some prebuilt R packages failed to install via conda."
    echo "    R setup_env.R will attempt source compilation as fallback."
  }
fi

echo "[*] Running R setup_env.R..."
"${CONDA_BIN}" run -n "${ENV_NAME}" Rscript "${ROOT_DIR}/r_scripts/setup_env.R"

if [[ "${RUN_DOCTOR}" -eq 1 ]]; then
  echo "[*] Running illumeta doctor..."
  "${CONDA_BIN}" run -n "${ENV_NAME}" python "${ROOT_DIR}/illumeta.py" doctor
fi

echo "[*] Full install completed."
