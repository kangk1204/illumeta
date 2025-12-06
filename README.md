# IlluMeta
*IlluMeta: Illuminating Methylation Analytics*

A compact, ready-to-run toolkit to go from GEO accession to interpretable methylation results in one command set. IlluMeta (Illuminating Methylation Analytics) automates IDAT retrieval, metadata cleanup, minfi + sesame processing, and generates an HTML dashboard ready to share. Works the same on Ubuntu, macOS (Intel/Apple Silicon), and Windows 11 via WSL2. Currently scoped to human (Homo sapiens) methylation arrays.

## What you get
- GEO download helper that pulls raw IDATs (or uses existing `idat/` if present) and writes `configure.tsv` + `configure_original.tsv`.
- Minfi + sesame processing, QC, DMP/DMR tables, and interactive plots.
- A one-page dashboard: `<Test>_vs_<Control>_results_index.html` linking to all outputs inside `<Test>_vs_<Control>_results/`.

## Prerequisites (pick your OS)
- **Ubuntu 22.04 / WSL2 (recommended)**  
  ```bash
  sudo apt-get update
  sudo apt-get install -y python3 python3-pip r-base pandoc \
      libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev
  ```
- **Windows 11 + WSL2**  
  Install WSL with Ubuntu, then run the same Ubuntu commands above inside the Ubuntu shell.
- **macOS (Apple Silicon or Intel)**  
  ```bash
  xcode-select --install                          # one-time build tools
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"  # if brew is missing
  brew install python r pandoc libxml2 openssl@3
  ```
  Make sure `R` and `Rscript` are on your PATH (restart the terminal after installing).
- **Python packages (use a venv to avoid PEP 668 errors)**  
  Homebrew’s Python on macOS is marked as “externally managed”, so install packages in a virtual environment:
  ```bash
  python3 -m venv .venv
  source .venv/bin/activate
  python -m pip install --upgrade pip
  python -m pip install requests
  ```
  If you insist on the system interpreter, use `python3 -m pip install --user requests` (or `--break-system-packages` at your own risk).
- **R packages**  
  First run of `illumeta.py` calls `r_scripts/setup_env.R` to install Bioconductor/CRAN deps (minfi, sesame, dmrff, etc.) and cache sesameData under `cache/`.  
  If the R library path is not writable, set a user library:
  ```bash
  export R_LIBS_USER="$HOME/R/library"
  mkdir -p "$R_LIBS_USER"
  ```
  Re-run setup manually with `Rscript r_scripts/setup_env.R`, or force a reinstall by prefixing commands with `ILLUMETA_FORCE_SETUP=1`.

## Quick start (same flow on Ubuntu/macOS/WSL2)
1. **Get the code**  
   ```bash
   git clone https://github.com/kangk1204/illumeta.git
   cd illumeta
   # Recommended for macOS/Homebrew Python to dodge PEP 668
   python3 -m venv .venv && source .venv/bin/activate
   ```
2. **Find a GSE with IDATs (optional)**  
   ```bash
   python3 illusearch.py --keywords "breast cancer" --email you@example.com -o geo_idat_methylation.tsv
   ```
   검색은 인간(Homo sapiens) 데이터에 한정됩니다(illumeta 분석 파이프라인과 일치).
3. **Download a dataset** (creates `idat/` + configs under the chosen folder):  
   ```bash
   python3 illumeta.py download GSE12345 -o projects/GSE12345
   ```
   If IDATs are already in `projects/GSE12345/idat/`, the download step is skipped.
4. **Label your groups**  
   Open `projects/GSE12345/configure.tsv` and fill the `primary_group` column (e.g., `Control`, `Case`). Use consistent spelling/casing—you will reference these labels in the next step. Optional covariates already present in the file can be kept; degenerate columns are dropped automatically.
5. **Run the analysis**  
   ```bash
   python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case
   ```
   Add `--tmp-dir /path/on/large/disk` if you want temp files off the system drive.
6. **View results**  
   Open `projects/GSE12345/Case_vs_Control_results_index.html` in a browser. All plots/tables live in `projects/GSE12345/Case_vs_Control_results/`.

## Using your own IDATs (no GEO download)
1. Place paired files as `*_Grn.idat` + `*_Red.idat` (compressed `.gz` also works) under `my_project/idat/`.
2. Copy `configure_original.tsv` (from any run) or create a small TSV with at least: `primary_group`, a sample ID column matching the IDAT basename (GSM/geo_accession works), and any covariates you want considered.
3. Run the analysis command pointing at `my_project` with your control/test labels.

## CLI options you might care about
- `illumeta.py analysis ... --pval 0.01 --lfc 1` to tighten thresholds.
- `--disable-auto-covariates` to skip automatic covariate selection from PCs.
- `--disable-sva` to skip surrogate variable analysis.
- `--include-covariates age,sex,batch` to force specific covariates (if present).
- `--permutations 1000` to estimate null DMP counts; `--vp-top 10000` to change variancePartition probe count.
- `--tmp-dir /mnt/drive/tmp` to write temp files to a larger disk.

## Notes and tips
- The pipeline caches data in `cache/` (kept alongside the repo). Delete it only if you want to reclaim space or force a recache.
- Large GEO folders (`GSE*/`) and `cache/` are ignored by git; keep raw data out of commits.
- Re-running with the same input will re-use installed R packages and cached annotations unless `ILLUMETA_FORCE_SETUP=1` is set.

## Troubleshooting (common fixes)
- **`R: command not found`**: install R (see prerequisites) and open a new shell.
- **`pandoc: command not found`**: install pandoc (`sudo apt-get install pandoc` or `brew install pandoc`).
- **`error: externally-managed-environment` from pip**: activate the repo venv (`source .venv/bin/activate`), then `python -m pip install --upgrade pip requests`. Homebrew’s system Python blocks global installs without a venv.
- **Missing R packages (e.g., `library(optparse)` failing)**: rerun setup with a user library:
  ```bash
  export R_LIBS_USER="$HOME/R/library"
  mkdir -p "$R_LIBS_USER"
  ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
  ```
  Then retry the download/analysis command.
- **R package install errors about libxml2/ssl/curl/icu**: on Ubuntu/WSL, install `libxml2-dev libcurl4-openssl-dev libssl-dev libicu-dev`; on macOS, `brew install libxml2 openssl@3` and retry with `R_LIBS_USER` set.
- **`library path not writable`**: set `R_LIBS_USER` as shown above, then rerun `Rscript r_scripts/setup_env.R`.
- **No samples matched your groups**: check spelling/case of `primary_group` values in `configure.tsv` and that IDAT basenames include the GSM/accession strings used in the file.
