# IlluMeta
*IlluMeta: Illuminating Methylation Analytics*

A compact, ready-to-run toolkit to go from GEO accession to interpretable methylation results in one command set. IlluMeta (Illuminating Methylation Analytics) automates IDAT retrieval, metadata cleanup, minfi + sesame processing, and generates an HTML dashboard ready to share. Works the same on Ubuntu, macOS (Intel/Apple Silicon), and Windows 11 via WSL2. Currently scoped to human (Homo sapiens) methylation arrays.

## What you get
- GEO download helper that pulls raw IDATs (or uses existing `idat/` if present) and writes `configure.tsv` + `configure_original.tsv`.
- Minfi + sesame processing, QC, DMP/DMR tables, and interactive plots.
- A one-page dashboard: `<Test>_vs_<Control>_results_index.html` linking to all outputs inside `<Test>_vs_<Control>_results/`.

## Setup: Ubuntu / WSL2 (recommended)
1) Install system deps:
   ```bash
   sudo apt-get update
   sudo apt-get install -y python3 python3-pip r-base pandoc cmake \
       libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev
   ```
2) Clone the repo and enter it:
   ```bash
   git clone https://github.com/kangk1204/illumeta.git
   cd illumeta
   ```
3) Create a Python venv in the repo (do this once, then just `source .venv/bin/activate` in future shells):
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install requests
   ```
   If you really want to use the system interpreter, use `python3 -m pip install --user requests` instead of global installs.
4) Optional: point R to a user-writable library to avoid permission issues:
   ```bash
   export R_LIBS_USER="$HOME/R/library"
   mkdir -p "$R_LIBS_USER"
   ```
5) First-time R deps (optional; `illumeta.py` will auto-run this on first use):
   ```bash
   ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
   ```

## Setup: macOS (Apple Silicon or Intel)
1) Install tools and Homebrew packages:
   ```bash
   xcode-select --install || true                        # one-time build tools
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"  # if brew is missing
   brew install python r pandoc gettext libxml2 openssl@3 curl pkg-config cmake
   ```
2) Help R/clang find Homebrew headers/libs (libintl.h, libxml2, OpenSSL, curl):
   ```bash
   export PATH="/opt/homebrew/opt/gettext/bin:$PATH"
   export PKG_CONFIG_PATH="/opt/homebrew/opt/libxml2/lib/pkgconfig:/opt/homebrew/opt/openssl@3/lib/pkgconfig:/opt/homebrew/opt/curl/lib/pkgconfig:$PKG_CONFIG_PATH"
   export CPPFLAGS="-I/opt/homebrew/opt/gettext/include -I/opt/homebrew/opt/libxml2/include ${CPPFLAGS}"
   export LDFLAGS="-L/opt/homebrew/opt/gettext/lib -L/opt/homebrew/opt/libxml2/lib -L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/curl/lib ${LDFLAGS}"
   ```
   Restart the shell so `R`/`Rscript` pick up your PATH.
3) Clone the repo and enter it:
   ```bash
   git clone https://github.com/kangk1204/illumeta.git
   cd illumeta
   ```
4) Create a Python venv in the repo (do this once; later just `source .venv/bin/activate`):
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install requests
   ```
   Homebrew’s Python is “externally managed”, so a venv avoids `externally-managed-environment` errors.
5) Optional: point R to a user-writable library:
   ```bash
   export R_LIBS_USER="$HOME/R/library"
   mkdir -p "$R_LIBS_USER"
   ```
6) First-time R deps (optional; `illumeta.py` will auto-run this on first use):
   ```bash
   ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
   ```

## Quick start (after setup; run from repo root with `.venv` activated)
1. **Enter the repo and activate your venv (don’t recreate it)**  
   ```bash
   cd illumeta
   source .venv/bin/activate
   ```
2. **Find a GSE with IDATs (optional)**  
   ```bash
   python3 illusearch.py --keywords "breast cancer" --email you@example.com -o geo_idat_methylation.tsv
   ```
   Search is limited to human (Homo sapiens) datasets to match the IlluMeta pipeline.
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
   `--group_con` and `--group_test` are required; use the exact labels from `configure.tsv`.
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
  If a download/analysis command fails mid-run, rerun the same command with `ILLUMETA_FORCE_SETUP=1` set.
  Then retry the download/analysis command.
- **Bioconductor install failed or partial install**: remove the repo-local library and rerun setup:
  ```bash
  rm -rf .r-lib
  ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
  ```
  IlluMeta defaults to `.r-lib`, so it will be recreated automatically on the next run.
- **`variancePartition`/`nloptr`/`lme4` complains about CMake**: install CMake so `nloptr` can build (Ubuntu/WSL: `sudo apt-get install cmake`, macOS: `brew install cmake`, or without sudo: `python3 -m pip install cmake` and ensure `~/.local/bin` or your venv’s `bin` is on `PATH`). Then rerun `Rscript r_scripts/setup_env.R`.
- **`failed to lock directory '.r-lib/00LOCK-...'`**: delete stale locks and retry:
  ```bash
  rm -rf .r-lib/00LOCK*
  ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
  ```
- **Segfault while installing `h5mread`/`Rhdf5lib` on macOS**: start with a clean user library and rerun setup to avoid stale compiled binaries:
  ```bash
  export R_LIBS_USER="$PWD/.r-lib"
  rm -rf "$R_LIBS_USER"
  mkdir -p "$R_LIBS_USER"
  ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
  ```
  Then re-run your `download`/`analysis` command from the repo root with the same `R_LIBS_USER` exported.
- **R package install errors about libxml2/ssl/curl/icu**: on Ubuntu/WSL, install `libxml2-dev libcurl4-openssl-dev libssl-dev libicu-dev`; on macOS, `brew install libxml2 openssl@3` and retry with `R_LIBS_USER` set.
- **macOS `fatal error: 'libintl.h' file not found` or `C++ compiler cannot create executables (BiocParallel)`**: install build tools and point R to Homebrew headers/libs before running setup:
  ```bash
  xcode-select --install || true
  brew install gettext libxml2 openssl@3 curl pkg-config
  export PATH="/opt/homebrew/opt/gettext/bin:$PATH"
  export PKG_CONFIG_PATH="/opt/homebrew/opt/libxml2/lib/pkgconfig:/opt/homebrew/opt/openssl@3/lib/pkgconfig:/opt/homebrew/opt/curl/lib/pkgconfig:$PKG_CONFIG_PATH"
  export CPPFLAGS="-I/opt/homebrew/opt/gettext/include -I/opt/homebrew/opt/libxml2/include ${CPPFLAGS}"
  export LDFLAGS="-L/opt/homebrew/opt/gettext/lib -L/opt/homebrew/opt/libxml2/lib -L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/curl/lib ${LDFLAGS}"
  ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
  ```
- **`library path not writable`**: set `R_LIBS_USER` as shown above, then rerun `Rscript r_scripts/setup_env.R`.
- **No samples matched your groups**: check spelling/case of `primary_group` values in `configure.tsv` and that IDAT basenames include the GSM/accession strings used in the file.
