# IlluMeta

IlluMeta (Illuminating Methylation Analytics) is an end-to-end pipeline that transforms a GEO accession (or raw IDAT files) into publication-ready DNA methylation results and an interactive HTML dashboard.
IlluMeta currently supports Illumina array IDATs (450k/EPIC/EPIC v2). Sequencing-based methylation (WGBS/methyl-capture) is on the roadmap but not yet implemented.

## Highlights
- **One-command workflow**: `download` → (auto-group or edit `configure.tsv`) → `analysis`.
- **Two independent pipelines**: **minfi (Noob)** and **sesame** run side-by-side; Sesame is reported in both **strict (Minfi-aligned)** and **native (pOOBAH-preserving)** views.
- **High-confidence consensus set (intersection)**: CpGs significant in **both** pipelines with the **same direction**. This is a conservative subset and may miss signals present in only one pipeline; review Minfi/Sesame results for sensitivity.
- **Sesame-native covariates**: when available, Sesame uses EpiDISH-based cell composition; otherwise Minfi-derived covariates are reused for stability.
- **Batch handling**: evaluates correction strategies (SVA/ComBat/limma) when a batch factor exists.
- **Primary-branch selection**: uses a documented multi-criteria scoring heuristic (not a global optimum); weights are saved in `analysis_parameters.json`.
- **CRF v2.1 robustness**: sample-size–adaptive robustness report (MMC/NCS/SSS) with tiered warnings plus lambda CI, MMC Spearman concordance, and SSS sign-consistency metrics.
- **Defensive stats**: guards against low-variance or single-group covariates to prevent hard failures in small studies.
- **Auto-group helper**: can populate `primary_group` from metadata when a reliable group column exists.
- **Auto-group caution**: auto-grouping is heuristic; always verify group labels/counts (see `preflight_report.json`, `decision_ledger.tsv`).
- **Paper-ready artifacts**: interactive HTML + static PNG figures, plus `methods.md`, `summary.json`, `analysis_parameters.json`, `sessionInfo.txt`, and `code_version.txt`.

## Quick start (conda, recommended)
1. Clone:
```bash
git clone https://github.com/kangk1204/illumeta.git
cd illumeta
```
2. Full install (EPIC v2 + devtools + clocks):
```bash
./scripts/install_full.sh
```
If you need R 4.5:
```bash
./scripts/install_full.sh --r45
```
3. Activate the environment:
```bash
conda activate illumeta
```
If you used `--r45`, activate `illumeta-r45` instead.
4. Run:
```bash
python3 illumeta.py download GSEXXXXX -o projects/GSEXXXXX
python3 illumeta.py analysis -c projects/GSEXXXXX/configure.tsv --auto-group
```

## Scope and supported data
- Supported: Illumina array IDATs (450k/EPIC/EPIC v2).
- Not yet supported: WGBS or methyl-capture sequencing (planned extension).
- Requirement: raw IDAT files must be available from GEO or provided manually.

## Auto-group limitations (must-read)
Auto-group is a convenience feature and **not** a substitute for manual group curation.
- It relies on **heuristics/metadata text** and can mislabel samples if metadata is incomplete or ambiguous.
- It is **not designed for multi-factor designs** (e.g., time-series + treatment + batch).
- Always verify `primary_group` counts in `Input_Group_Distribution.csv` and spot-check sample metadata before interpreting results.

## Installation

IlluMeta is easiest to install with **conda**. If you are on **Windows**, we strongly recommend using **WSL2 (Ubuntu)** for the most reliable R package installation.

### Quick full install (conda, all OS)
This script creates or updates the conda env, installs Python deps, runs `setup_env.R`, then finishes with `illumeta.py doctor`.
It requires `conda` or `mamba` in your PATH (Miniforge recommended).
```bash
./scripts/install_full.sh
```
Options:
- `--r45` uses `environment-r45.yml` (R 4.5).
- `--env-file PATH` uses a custom conda env file.
- `--env NAME` overrides the env name.
- `--skip-doctor` skips the final `illumeta.py doctor` check.

Logs are saved to `projects/illumeta_install_full_*.log`.
If `illumeta.py doctor` reports missing **optional** packages, core features still work; those optional features will be skipped.

Before you start (important):
- **Pick one environment**: conda **or** system R. Do **not** mix them.
- If using conda, always run `conda activate illumeta` before `Rscript`.
- If you want to use system R instead, run `conda deactivate` first.
- Quick check (should point inside your chosen environment):
```bash
which R
R -q -e 'cat(R.version.string, "\n"); cat(R.home(), "\n")'
```
- If you previously installed with a different R version, run the **clean install** command below (it removes mismatched packages automatically).

Choose your OS section below and follow steps 0-4. The commands are intentionally repeated for clarity.

### Ubuntu (native or WSL2)

#### 0) Install prerequisites (Git + Conda)
```bash
sudo apt-get update
sudo apt-get install -y git curl
```
If you plan to use system R or install `devtools`/`roxygen2` (xml2/XML), install system libraries once:
```bash
sudo apt-get install -y \
  build-essential cmake git pandoc pkg-config gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev \
  libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
  libgit2-dev
```
Skipping this often causes `xml2`/`roxygen2` install errors.

Install **Miniforge (conda)**:
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source ~/.bashrc
conda --version
```
If your Ubuntu is ARM64, use `Miniforge3-Linux-aarch64.sh` instead.

#### 1) Clone the repository
```bash
git clone https://github.com/kangk1204/illumeta.git
cd illumeta
```

#### 2) Create an environment (choose one)

##### Option A: Conda (recommended for beginners)
This uses conda to provide R/Python plus the system libraries needed by many R packages.

Recommended (default conda-forge):
```bash
conda env create -f environment.yml
conda activate illumeta
```
If you need the latest GitHub `sesame` (requires R >= 4.5), use:
```bash
conda env create -f environment-r45.yml
conda activate illumeta-r45
```
If conda solve fails, try:
```bash
conda env create -f environment.yml --solver=classic
```
If you already created the environment before updating IlluMeta, refresh it:
```bash
conda env update -f environment.yml --prune
```
Optional sanity check (conda):
```bash
which R
R -q -e 'cat(R.home(), "\n")'
```
`which R` should point inside the `illumeta` conda env.

##### Option B: venv + system R
###### Prerequisites
- **Python** 3.8+
- **R** 4.4+ recommended (required for EPIC v2; R 4.3 works for 450k/EPIC only)
- **pandoc** (required for self-contained HTML reports via `htmlwidgets`)
- **System libraries** (for compiling common R packages)

Install system libraries:
```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake git pandoc pkg-config gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev \
  libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
  libgit2-dev
```

###### Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

#### 3) Install R dependencies (first run / CI / reproducible setup)
Before running `setup_env.R` (important for `devtools`/`roxygen2`/`xml2`):
Conda (recommended):
```bash
conda env update -f environment.yml --prune
# Or, minimal install without full env update:
conda install -c conda-forge libxml2-devel zlib xz pkg-config
# devtools/tidyverse extras:
conda install -c conda-forge libgit2 harfbuzz fribidi fontconfig expat
```
System R (Ubuntu/WSL):
```bash
sudo apt-get install -y libxml2-dev zlib1g-dev liblzma-dev pkg-config
# devtools/tidyverse extras:
sudo apt-get install -y libgit2-dev libharfbuzz-dev libfribidi-dev libfontconfig1-dev libexpat1-dev
```

This installs required R/Bioconductor packages into the repo-local library (`.r-lib/R-<major.minor>`) and may take ~10-30 minutes.

Recommended (Beginner, core features):
- Works on R 4.3+ (EPIC v2 is optional).
```bash
# Ensure your environment is activated first:
# - Conda: conda activate illumeta
# - venv:  source .venv/bin/activate

ILLUMETA_CLEAN_MISMATCHED_RLIB=1 \
ILLUMETA_DOWNLOAD_RETRIES=3 \
ILLUMETA_FORCE_SETUP=1 \
Rscript r_scripts/setup_env.R
```

Full install (all optional features, including EPIC v2, devtools, and clocks):
- Requires R 4.4+ for EPIC v2 and extra system libraries for devtools/tidyverse.
```bash
ILLUMETA_CLEAN_MISMATCHED_RLIB=1 \
ILLUMETA_DOWNLOAD_RETRIES=3 \
ILLUMETA_FORCE_SETUP=1 ILLUMETA_INSTALL_DEVTOOLS=1 ILLUMETA_INSTALL_CLOCKS=1 ILLUMETA_REQUIRE_EPICV2=1 \
Rscript r_scripts/setup_env.R
```

Optional toggles:
- `ILLUMETA_REQUIRE_EPICV2=1` (requires R 4.4+ / Bioconductor 3.19+)
- `ILLUMETA_INSTALL_DEVTOOLS=1` (devtools/tidyverse; requires `libgit2` + `pkg-config` and harfbuzz/fribidi/fontconfig/expat for textshaping/ragg; conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge libgit2 harfbuzz fribidi fontconfig expat`)
- `ILLUMETA_INSTALL_CLOCKS=1` (methylclock/planet/wateRmelon)
- `ILLUMETA_CLEAN_MISMATCHED_RLIB=1` (use after switching R versions)
- `ILLUMETA_DOWNLOAD_RETRIES=2` (retry big downloads if they fail; increase to 3+ on flaky networks)

If you see errors like `libxml-2.0`, `xml2`, `lzma.h`, or `zlib.h` not found:
- Conda (Linux/WSL): `conda env update -f environment.yml --prune` (or `conda install -c conda-forge libxml2-devel zlib xz pkg-config`)
- System R (Ubuntu/WSL): `sudo apt-get install -y libxml2-dev zlib1g-dev liblzma-dev pkg-config`
Then rerun the setup command.

#### 4) Check your environment (recommended)
```bash
python3 illumeta.py doctor
```

### macOS (Apple Silicon M1-M4)

#### 0) Install prerequisites (Git + Conda)
Install Xcode Command Line Tools (includes git and compilers):
```bash
xcode-select --install
```

Install Homebrew (if you do not already have it):
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
If you plan to use system R or install `devtools`/`roxygen2` (xml2/XML), install system libraries once:
```bash
brew install cmake git pandoc pkg-config openssl@3 libxml2 freetype libpng libtiff jpeg webp harfbuzz fribidi fontconfig libgit2 libomp gcc
```
Skipping this often causes `xml2`/`roxygen2` install errors.

Install **Miniforge (conda)**:
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```
If you use bash, replace `zsh` with `bash` and source `~/.bashrc` instead.

#### 1) Clone the repository
```bash
git clone https://github.com/kangk1204/illumeta.git
cd illumeta
```

#### 2) Create an environment (choose one)

##### Option A: Conda (recommended for beginners)
This uses conda to provide R/Python plus the system libraries needed by many R packages.

Recommended (default conda-forge):
```bash
conda env create -f environment.yml
conda activate illumeta
```
If you need the latest GitHub `sesame` (requires R >= 4.5), use:
```bash
conda env create -f environment-r45.yml
conda activate illumeta-r45
```
If conda solve fails, try:
```bash
conda env create -f environment.yml --solver=classic
```
If you already created the environment before updating IlluMeta, refresh it:
```bash
conda env update -f environment.yml --prune
```
Optional sanity check (conda):
```bash
which R
R -q -e 'cat(R.home(), "\n")'
```
`which R` should point inside the `illumeta` conda env.

##### Option B: venv + system R
###### Prerequisites
- **Python** 3.8+
- **R** 4.4+ recommended (required for EPIC v2; R 4.3 works for 450k/EPIC only)
- **pandoc** (required for self-contained HTML reports via `htmlwidgets`)
- **System libraries** (for compiling common R packages)

Install system libraries:
```bash
brew install cmake git pandoc pkg-config openssl@3 libxml2 freetype libpng libtiff jpeg webp harfbuzz fribidi fontconfig libgit2 libomp gcc
```

###### Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

#### 3) Install R dependencies (first run / CI / reproducible setup)
Before running `setup_env.R` (important for `devtools`/`roxygen2`/`xml2`):
Conda (recommended):
```bash
conda env update -f environment.yml --prune
# Or, minimal install without full env update:
conda install -c conda-forge libxml2-devel zlib xz pkg-config
# devtools/tidyverse extras:
conda install -c conda-forge libgit2 harfbuzz fribidi fontconfig expat
```
System R (macOS):
```bash
brew install libxml2 zlib xz pkg-config
# devtools/tidyverse extras:
brew install libgit2 harfbuzz fribidi fontconfig expat
```

This installs required R/Bioconductor packages into the repo-local library (`.r-lib/R-<major.minor>`) and may take ~10-30 minutes.

Recommended (Beginner, core features):
- Works on R 4.3+ (EPIC v2 is optional).
```bash
# Ensure your environment is activated first:
# - Conda: conda activate illumeta
# - venv:  source .venv/bin/activate

ILLUMETA_CLEAN_MISMATCHED_RLIB=1 \
ILLUMETA_DOWNLOAD_RETRIES=3 \
ILLUMETA_FORCE_SETUP=1 \
Rscript r_scripts/setup_env.R
```

Full install (all optional features, including EPIC v2, devtools, and clocks):
- Requires R 4.4+ for EPIC v2 and extra system libraries for devtools/tidyverse.
```bash
ILLUMETA_CLEAN_MISMATCHED_RLIB=1 \
ILLUMETA_DOWNLOAD_RETRIES=3 \
ILLUMETA_FORCE_SETUP=1 ILLUMETA_INSTALL_DEVTOOLS=1 ILLUMETA_INSTALL_CLOCKS=1 ILLUMETA_REQUIRE_EPICV2=1 \
Rscript r_scripts/setup_env.R
```

Optional toggles:
- `ILLUMETA_REQUIRE_EPICV2=1` (requires R 4.4+ / Bioconductor 3.19+)
- `ILLUMETA_INSTALL_DEVTOOLS=1` (devtools/tidyverse; requires `libgit2` + `pkg-config` and harfbuzz/fribidi/fontconfig/expat for textshaping/ragg; conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge libgit2 harfbuzz fribidi fontconfig expat`)
- `ILLUMETA_INSTALL_CLOCKS=1` (methylclock/planet/wateRmelon)
- `ILLUMETA_CLEAN_MISMATCHED_RLIB=1` (use after switching R versions)
- `ILLUMETA_DOWNLOAD_RETRIES=2` (retry big downloads if they fail; increase to 3+ on flaky networks)

If you see errors like `libxml-2.0`, `xml2`, `lzma.h`, or `zlib.h` not found:
- Conda (macOS): `conda env update -f environment.yml --prune` (or `conda install -c conda-forge libxml2-devel zlib xz pkg-config`)
- System R (macOS): `brew install libxml2 zlib xz pkg-config`
Then rerun the setup command.

#### 4) Check your environment (recommended)
```bash
python3 illumeta.py doctor
```

### Windows 11 (WSL2 + Ubuntu)

#### 0) Install prerequisites (WSL2 + Git + Conda)
In **PowerShell (Admin)**:
```powershell
wsl --install -d Ubuntu
```
Reboot if prompted, then open **Ubuntu (WSL2)**.

In the Ubuntu terminal, install Git + curl:
```bash
sudo apt-get update
sudo apt-get install -y git curl
```
If you plan to use system R or install `devtools`/`roxygen2` (xml2/XML), install system libraries once:
```bash
sudo apt-get install -y \
  build-essential cmake git pandoc pkg-config gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev \
  libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
  libgit2-dev
```
Skipping this often causes `xml2`/`roxygen2` install errors.

Install **Miniforge (conda)** inside WSL:
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source ~/.bashrc
conda --version
```

#### 1) Clone the repository
```bash
git clone https://github.com/kangk1204/illumeta.git
cd illumeta
```

#### 2) Create an environment (choose one)

##### Option A: Conda (recommended for beginners)
This uses conda to provide R/Python plus the system libraries needed by many R packages.

Recommended (default conda-forge):
```bash
conda env create -f environment.yml
conda activate illumeta
```
If you need the latest GitHub `sesame` (requires R >= 4.5), use:
```bash
conda env create -f environment-r45.yml
conda activate illumeta-r45
```
If conda solve fails, try:
```bash
conda env create -f environment.yml --solver=classic
```
If you already created the environment before updating IlluMeta, refresh it:
```bash
conda env update -f environment.yml --prune
```
Optional sanity check (conda):
```bash
which R
R -q -e 'cat(R.home(), "\n")'
```
`which R` should point inside the `illumeta` conda env.

##### Option B: venv + system R
###### Prerequisites
- **Python** 3.8+
- **R** 4.4+ recommended (required for EPIC v2; R 4.3 works for 450k/EPIC only)
- **pandoc** (required for self-contained HTML reports via `htmlwidgets`)
- **System libraries** (for compiling common R packages)

Install system libraries:
```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake git pandoc pkg-config gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev \
  libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
  libgit2-dev
```

###### Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

#### 3) Install R dependencies (first run / CI / reproducible setup)
Before running `setup_env.R` (important for `devtools`/`roxygen2`/`xml2`):
Conda (recommended):
```bash
conda env update -f environment.yml --prune
# Or, minimal install without full env update:
conda install -c conda-forge libxml2-devel zlib xz pkg-config
# devtools/tidyverse extras:
conda install -c conda-forge libgit2 harfbuzz fribidi fontconfig expat
```
System R (WSL/Ubuntu):
```bash
sudo apt-get install -y libxml2-dev zlib1g-dev liblzma-dev pkg-config
# devtools/tidyverse extras:
sudo apt-get install -y libgit2-dev libharfbuzz-dev libfribidi-dev libfontconfig1-dev libexpat1-dev
```

This installs required R/Bioconductor packages into the repo-local library (`.r-lib/R-<major.minor>`) and may take ~10-30 minutes.

Recommended (Beginner, core features):
- Works on R 4.3+ (EPIC v2 is optional).
```bash
# Ensure your environment is activated first:
# - Conda: conda activate illumeta
# - venv:  source .venv/bin/activate

ILLUMETA_CLEAN_MISMATCHED_RLIB=1 \
ILLUMETA_DOWNLOAD_RETRIES=3 \
ILLUMETA_FORCE_SETUP=1 \
Rscript r_scripts/setup_env.R
```

Full install (all optional features, including EPIC v2, devtools, and clocks):
- Requires R 4.4+ for EPIC v2 and extra system libraries for devtools/tidyverse.
```bash
ILLUMETA_CLEAN_MISMATCHED_RLIB=1 \
ILLUMETA_DOWNLOAD_RETRIES=3 \
ILLUMETA_FORCE_SETUP=1 ILLUMETA_INSTALL_DEVTOOLS=1 ILLUMETA_INSTALL_CLOCKS=1 ILLUMETA_REQUIRE_EPICV2=1 \
Rscript r_scripts/setup_env.R
```

Optional toggles:
- `ILLUMETA_REQUIRE_EPICV2=1` (requires R 4.4+ / Bioconductor 3.19+)
- `ILLUMETA_INSTALL_DEVTOOLS=1` (devtools/tidyverse; requires `libgit2` + `pkg-config` and harfbuzz/fribidi/fontconfig/expat for textshaping/ragg; conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge libgit2 harfbuzz fribidi fontconfig expat`)
- `ILLUMETA_INSTALL_CLOCKS=1` (methylclock/planet/wateRmelon)
- `ILLUMETA_CLEAN_MISMATCHED_RLIB=1` (use after switching R versions)
- `ILLUMETA_DOWNLOAD_RETRIES=2` (retry big downloads if they fail; increase to 3+ on flaky networks)

If you see errors like `libxml-2.0`, `xml2`, `lzma.h`, or `zlib.h` not found:
- Conda (Linux/WSL): `conda env update -f environment.yml --prune` (or `conda install -c conda-forge libxml2-devel zlib xz pkg-config`)
- System R (Ubuntu/WSL): `sudo apt-get install -y libxml2-dev zlib1g-dev liblzma-dev pkg-config`
Then rerun the setup command.

#### 4) Check your environment (recommended)
```bash
python3 illumeta.py doctor
```

Notes on R libraries:
- IlluMeta stores packages in per-R-version subfolders (e.g., `.r-lib/R-4.4` or `~/.illumeta/r-lib/R-4.4`) to avoid mixing binaries across R versions. Switching R versions will create a new folder; rerun setup (you can delete old folders to reclaim space).
- IlluMeta **defaults to a repo-local** R library under `.r-lib/R-<major.minor>` for reproducibility.
- To keep an externally-set `R_LIBS_USER`, run with `ILLUMETA_RESPECT_R_LIBS_USER=1`.
- On macOS, if the project is on an external volume (e.g. `/Volumes/...`), IlluMeta installs to `~/.illumeta/r-lib/R-<major.minor>` to avoid I/O errors. Set `ILLUMETA_ALLOW_EXTERNAL_LIB=1` to force `.r-lib` on the external drive.
- If you are using a conda R and need conda system libs, set `ILLUMETA_USE_CONDA_LIBS=1`.
- To clean mismatched packages after switching R versions, run:
  `ILLUMETA_CLEAN_MISMATCHED_RLIB=1 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.

## Docker (optional)
Build the container:
```bash
docker build -t illumeta .
```
Run a quick environment check:
```bash
docker run --rm -it -v "$PWD":/app illumeta doctor
```
Run analysis (mount your project directory):
```bash
docker run --rm -it -v "$PWD":/app illumeta analysis \
  -i projects/GSE121633 \
  --group_con Control \
  --group_test Case
```

## Quick start

### 0) Choose your input type
- If you have **GEO**: follow steps 1–4 below.
- If you have **your own IDATs**: see “Analyze your own IDATs (non-GEO)” right after this section.

### 1) Download a dataset from GEO
```bash
# Activate your environment (choose one)
# - Conda: conda activate illumeta
# - venv:  source .venv/bin/activate

python3 illumeta.py download GSE121633 -o projects/GSE121633
# If a GEO series has multiple platforms, force one by GPL ID:
python3 illumeta.py download GSE121633 -o projects/GSE121633 --platform GPL21145
```

### 2) Assign groups (manual or auto)
Option A (manual):
Edit `projects/GSE121633/configure.tsv` and fill in `primary_group` (e.g., `Control` / `Case`).
`configure.tsv` must be **tab-delimited (TSV)**; CSV is not supported.

Option B (auto-group on analysis):
```bash
python3 illumeta.py analysis \
  -i projects/GSE121633 \
  --group_con Control \
  --group_test Case \
  --auto-group \
  --group-column disease_state
```
If your grouping is encoded in GEO characteristics, you can use `--group-key` (e.g., `--group-key disease`).

### 3) Run analysis
```bash
python3 illumeta.py analysis \
  -i projects/GSE121633 \
  --group_con Control \
  --group_test Case
```
If you chose auto-grouping above, include the same `--auto-group` flags here instead of editing the TSV.
Note: the default output folder name is derived from the group labels. If it contains non-ASCII characters, IlluMeta normalizes it to a safe ASCII name for filesystem compatibility (the dashboard filename follows the folder name).
Optional (signal preservation checks):
```bash
# Provide a CpG marker list (TSV/CSV with CpG column or one CpG per line)
python3 illumeta.py analysis \
  -i projects/GSE121633 \
  --group_con Control \
  --group_test Case \
  --marker-list markers.tsv
```
This generates `*_Signal_Preservation.csv` and (if provided) `*_Known_Marker_Summary.csv`.

### 4) Open the dashboard
Open the generated HTML:
`projects/GSE121633/Case_vs_Control_results_index.html`

### 5) Interpret results (beginner checklist)
1. **QC first**: open `QC_Summary.csv`, `Sample_QC_DetectionP_FailFraction.html`, and `Sample_QC_Intensity_Medians.html`. If many samples fail QC, stop and fix QC before interpreting DMP/DMR.
2. **Batch check**: compare `*_Batch_Evaluation_Before.html` vs `*_Batch_Evaluation_After.html`. After-correction should reduce batch association without erasing group signal. The chosen method (SVA/ComBat/limma/none) and covariates are logged in `decision_ledger.tsv`.
3. **Consensus vs pipeline results**:  
   - **Consensus (Strict)** = Minfi ∩ Sesame (strict alignment); most conservative.  
   - **Consensus (Native)** = Minfi ∩ Sesame native; more sensitive to Sesame’s masking philosophy.  
4. **DMP tables**: use `*_Top_DMPs.html` for top hits (sorted by P.Value). Full tables are `*_DMPs_full.csv`.
5. **DMR tables**: use `*_DMRs_Table.html` / `*_DMRs.csv` (sorted by p.value). Check volcano/manhattan for signal distribution.

## Usage

### Search GEO for IDAT-enabled datasets
Beginner-friendly flow:
1) Activate your environment (conda or venv).
2) Run a search with simple keywords.
3) Open the TSV and pick a GSE ID to analyze.

```bash
# Required: --keywords (use quotes for multi-word queries)
python3 illumeta.py search --keywords "breast cancer" --email your_email@example.com -o search_results.tsv

# Faster (skip FTP supplement checks)
python3 illumeta.py search --keywords "breast cancer" --no-check-suppl -o search_results.tsv
```
Output columns:
- `gse_id`: GEO Series ID to use with `illumeta.py download`
- `platform_type`: 450k / EPIC (850k) / EPIC v2 (~936k; sometimes called “950k”)
- `suppl_has_idat`: yes/no/error (only if supplement checks are enabled)
IlluMeta requires raw IDATs; `download` will stop if no IDATs are available.

If you see `Error: Python package 'requests' is missing`, install it via:
`pip install -r requirements.txt` (inside your env).

### Check installation (illumeta doctor)
Use this before a first run or when something looks wrong. It does not install anything.
```bash
python3 illumeta.py doctor
```
How to read the output:
- **Core R packages: OK** = ready to run analysis.
- **Optional R packages missing** = some optional features (e.g., clocks) are disabled; this is safe for most users.
Note: For EPIC v2 (~936k) datasets, IlluMeta requires EPICv2 manifest/annotation packages and will stop if they are missing. If you plan to run EPIC v2, set `ILLUMETA_REQUIRE_EPICV2=1` before running setup.


### Analyze your own IDATs (non-GEO)
1. Create `my_project/idat/` and place `_Grn.idat` / `_Red.idat` pairs there.
2. Create `my_project/configure.tsv` with at least (tab-delimited TSV; CSV is not supported):
   - `Basename` (e.g., `idat/Sample1_R01C01`)
   - `primary_group` (e.g., `Treated`, `Untreated`)
   - Optional: `tissue` (e.g., `Blood`, `CordBlood`, `Placenta`) to auto-select reference (Placenta uses planet reference when available)
   - Example:
     ```tsv
     Basename	primary_group	SampleID	sex	age
     idat/S1_R01C01	control	S1	M	52
     idat/S1_R02C01	control	S2	F	47
     idat/S2_R01C01	case	S3	M	61
     ```
   - If you prefer auto-grouping, you may leave `primary_group` empty and pass `--auto-group` with `--group-column` when running analysis.
3. Run:
```bash
python3 illumeta.py analysis -i my_project --group_con Untreated --group_test Treated
```

### Auto-grouping (optional)
If your metadata already contains a reliable group column, IlluMeta can populate `primary_group` automatically:
```bash
python3 illumeta.py analysis -i projects/GSE12345 \
  --group_con Control --group_test Case \
  --auto-group --group-column disease_state \
  --group-map "normal=Control,tumor=Case"
```
For GEO characteristics, use a key (parsed from `key: value` patterns):
```bash
python3 illumeta.py analysis -i projects/GSE12345 \
  --group_con Control --group_test Case \
  --auto-group --group-key disease
```
Auto-detection is conservative and will stop if no clear group signal is found; in that case specify `--group-column` or edit `configure.tsv`.
Auto-group decisions and warnings are recorded in `preflight_report.json` and `decision_ledger.tsv`.
Auto-grouping is heuristic; always verify group counts and labels before interpreting results, especially when multiple categorical columns exist.
Auto-group prioritizes columns with high coverage and low category counts; numeric-coded categories (e.g., 0/1) are supported, while highly fragmented columns or heavy missingness are rejected.

### Common analysis options
```bash
# Auto-group from a metadata column
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case \
  --auto-group --group-column disease_state --group-map "normal=Control,tumor=Case"

# Tighten thresholds
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --pval 0.01 --lfc 1.0

# Add an absolute Delta Beta filter
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --delta-beta 0.05

# Beginner-safe mode (stricter group checks + conservative thresholds)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --beginner-safe

# Disable SVA (very small n)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --disable-sva

# Provide covariates to always try (if present in configure.tsv)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --include-covariates age,sex

# Use epigenetic clocks as candidate covariates (auto-selected if relevant)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --include-clock-covariates

# Placenta reference cell composition + clocks (planet)
python3 illumeta.py analysis -i projects/GSE307314 --group_con control --group_test test --tissue Placenta

# Auto tissue inference from configure.tsv (falls back to RefFreeEWAS)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --tissue Auto

# Custom cell reference (package, package::object, or .rds/.rda)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case \
  --tissue Blood --cell-reference FlowSorted.Blood.EPIC

# Custom cell reference with explicit platform string
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case \
  --tissue Blood --cell-reference refs/blood_ref.rds --cell-reference-platform IlluminaHumanMethylationEPIC

# Enable sesame dyeBiasCorrTypeINorm (disabled by default for stability)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --sesame-typeinorm

# Mixed-array safeguard override (only if you know what you're doing)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --force-idat

# Permutation-based null test ("perm")
# Note: there is no `illumeta perm` command; use --permutations with analysis.
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --permutations 50
```

### Permutation test ("perm") for beginners
This checks how many significant DMPs you would get by chance if labels were shuffled.

How to run it:
```bash
python3 illumeta.py analysis \
  -i projects/GSE12345 \
  --group_con Control \
  --group_test Case \
  --permutations 50
```

### Benchmarking & Robustness (Recommended for papers)
Use these scripts to report objective QC, batch correction, and robustness metrics across datasets.

```bash
# 1) Build a benchmark summary table across multiple GSE runs
python3 scripts/build_benchmark_table.py \
  --input-tsv benchmarks/benchmark_inputs.tsv \
  --projects-root projects \
  --out benchmarks/benchmark_summary.tsv

# 2) Generate a primary-branch summary table + figure
python3 scripts/build_benchmark_figures.py \
  --input-tsv benchmarks/benchmark_summary.tsv \
  --out-dir benchmarks

# 3) Quantify batch correction (PVCA before/after)
python3 scripts/build_correction_summary.py \
  --root projects \
  --out-dir benchmarks/correction_summary

# 4) Run ablation variants (SVA on/off, auto-covariates on/off, etc.)
python3 scripts/ablation_runner.py \
  -i projects/GSE12345 \
  --group-con Control \
  --group-test Case \
  --out-root ablation_runs/GSE12345

# Optional: include sesame TypeINorm variant (experimental)
python3 scripts/ablation_runner.py \
  -i projects/GSE12345 \
  --group-con Control \
  --group-test Case \
  --variants baseline,sesame_typeinorm
```

Notes:
- **TypeINorm is disabled by default** for stability (thread errors observed on some systems). Enable with `--sesame-typeinorm` when needed.
- **Sesame runs single-thread by default** for stability. If you want multi-threading, set `ILLUMETA_SESAME_SINGLE_THREAD=0` in your shell.
- **Auto-grouping is heuristic**: always verify group labels/counts before interpreting results.
- **Small n**: consider `--disable-sva` and report this choice explicitly.
Tips:
- Start small (e.g., 10–50). Larger values take longer.
- Results are saved as `*_Permutation_Results.csv` and `*_Permutation_Summary.csv`.
- The dashboard shows **Perm mean/max sig (null)** so you can compare against real results.

### Robustness / ablation (for manuscripts)
IlluMeta writes per-pipeline metrics to help assess robustness:
- `Minfi_Metrics.csv`, `Sesame_Metrics.csv`, `Sesame_Native_Metrics.csv` (lambda, batch stats, SV/covariate usage, permutation stats)
- `*_Ablation_Summary.csv` (raw vs corrected batch metrics)
- `*_Permutation_Results.csv` and `*_Permutation_Summary.csv` (null calibration)
- `*_LambdaGuard_*` (only when lambda exceeds threshold; group-only check on pre-correction betas)

To run a standard ablation suite:
```bash
python3 scripts/ablation_runner.py -i projects/GSE12345 \
  --group-con Control --group-test Case
```
Outputs:
- `ablation_manifest.tsv`: run status per variant
- `ablation_counts.tsv`: DMP counts across pipelines
- `ablation_metrics_long.tsv`: long-form metrics table
- `ablation_parameters.json`: exact parameters per variant

Lambda guard + variancePartition autoscale settings live in `config.yaml` (copy from `config.yaml.template`):
```yaml
lambda_guard:
  enabled: true
  threshold: 1.5
  min_samples: 8
  action: warn_simplify

variance_partition:
  autoscale_numeric: true
  autoscale_on_fail: true
```

### Multi-GEO benchmark table
Build a long-form benchmark table (one row per pipeline per dataset):
```bash
python3 scripts/build_benchmark_table.py \
  --input-tsv geo_idat_methylation.tsv \
  --projects-root projects/bench_top5 \
  --output-tsv benchmarks/benchmark_summary.tsv \
  --output-md benchmarks/benchmark_summary.md
```
Outputs:
- `benchmarks/benchmark_summary.tsv`: metrics table for downstream plotting
- `benchmarks/benchmark_summary.md`: quick stats + run metadata
The summary table now includes objective QC/robustness indicators such as:
- QC pass rate (samples) and probe retention
- Batch signal reduction before/after correction
- Lambda / permutation calibration (KS, lambda)
- Variance partition signal for primary_group
- Tier3/DMR status and primary-branch overrides (if any)

### Paper-ready summary figure
Generate a compact table + figure from the benchmark TSV (primary branch only by default):
```bash
python3 scripts/build_benchmark_figures.py \
  --input-tsv benchmarks/benchmark_summary.tsv \
  --out-dir benchmarks
```
Outputs:
- `benchmarks/benchmark_primary_summary.tsv`
- `benchmarks/benchmark_primary_summary.md`
- `benchmarks/benchmark_primary_summary.png`
Tip: add `--all-rows` to include every pipeline row instead of collapsing to the primary branch.

### Smoke tests (multi-run)
Prepare a TSV with columns `config`, `group_con`, `group_test` (optional: `name`, `output`, `extra_args`) and run:
```bash
python3 scripts/run_smoke_pipeline.py \
  --jobs benchmarks/smoke_jobs.tsv \
  --out-dir benchmarks/smoke_runs \
  --illumeta illumeta.py
```
Outputs:
- `benchmarks/smoke_runs/smoke_report.tsv`: per-run status summary
- `benchmarks/smoke_runs/logs/*.log`: full stdout/stderr for each run

## Outputs (what to use in a paper)

IlluMeta writes results to the analysis output directory (default: `[input]/[test]_vs_[con]_results/`).

### Reproducibility
- `methods.md`: auto-generated methods summary for the run.
- `summary.json`: primary result mode (tier3 vs standard), Tier3 meta method (fixed/random/auto), lambda guard status, DMR status, and sesame normalization notes.
- `analysis_parameters.json`: thresholds, presets, and scoring weights used for optimization.
- `*_BetaMatrix.tsv.gz`: processed beta matrix used for modeling (one per pipeline: Minfi, Sesame, Sesame_Native).
- `*_MvalueMatrix.tsv.gz`: processed M-value matrix (logit of the modeling betas) used for statistical tests.
- `*_BetaMatrix_PreFilter.tsv.gz`: pre-filter beta matrix (after normalization, before sample/probe QC; all samples; Minfi + Sesame).
- `*_MvalueMatrix_PreFilter.tsv.gz`: pre-filter M-value matrix (logit of pre-filter betas; before sample/probe QC; all samples; Minfi + Sesame).
- `Minfi_MethylatedSignal_PreFilter.tsv.gz`: pre-filter methylated signal matrix (after normalization, before sample/probe QC; all samples).
- `Minfi_UnmethylatedSignal_PreFilter.tsv.gz`: pre-filter unmethylated signal matrix (after normalization, before sample/probe QC; all samples).
- `Minfi_DetectionP_PreFilter.tsv.gz`: pre-filter detection P-value matrix (all samples; GEO processed/normalized requirement).

### GEO submission helper
Generate a manifest (and optional bundle) of processed files for GEO:
```bash
python3 scripts/prepare_geo_submission.py \
  --result-dir projects/GSE12345/Case_vs_Control_results \
  --copy
```
Outputs:
- `geo_submission/geo_submission_manifest.tsv` (file list + descriptions)
- `geo_submission/*` (copied processed files when `--copy` is set)
- `sessionInfo.txt`: full R session and package versions.
- `code_version.txt`: git commit hash (when available).
- `config_used.yaml`: resolved config and preset details.
- `decision_ledger.tsv`: automated decision log (covariates, batch strategy, consensus branch).
- `preflight_report.json`: preflight summary (includes auto-group info when used).
- `Correction_Adequacy_Report.txt`: Correction Adequacy Framework (CAF) report for the primary branch.
- `Correction_Adequacy_Summary.csv`: CAF component scores for the primary branch.

### QC
- `Preflight_Summary.csv`: preflight checks and group counts before processing.
- `Preflight_IDAT_Pairs.csv`: per-sample IDAT pair existence check.
- `QC_Summary.csv`: sample/probe QC counts.
- `Sample_QC_Metrics.csv`: per-sample QC metrics.
- `Sample_QC_DetectionP_FailFraction.png` and `Sample_QC_Intensity_Medians.png`: static QC figures (HTML versions are also saved).
- Defaults: detection P threshold 0.05, sample fail fraction 0.20, intensity median threshold log2 9.0 (override with `--qc-detection-p-threshold`, `--qc-sample-fail-frac`, `--qc-intensity-threshold`).

### Differential methylation (per pipeline)
For each pipeline (`Minfi`, `Sesame` = strict/Minfi-aligned, `Sesame_Native` = native Sesame):
- `*_DMPs_full.csv`: full differential methylation results table.
- `*_Volcano.html/.png`, `*_Manhattan.html/.png`, `*_QQPlot.html/.png`
- `*_Top100_Heatmap.html/.png`
- `*_DMR_Volcano.html/.png`, `*_DMR_Manhattan.html/.png`, `*_Top_DMRs_Heatmap.html/.png`
- If Tier3 confounding is detected, DMRs are emitted from the stratified/meta primary results as `*_Tier3_Primary_DMRs.*`.

### Consensus (intersection) call set
- `Intersection_Consensus_DMPs.csv` and `Intersection_Consensus_DMPs.html` (strict)
- `Intersection_Native_Consensus_DMPs.csv` and `Intersection_Native_Consensus_DMPs.html` (native)
- `Intersection*_LogFC_Concordance.html/.png` (minfi vs sesame logFC concordance)
- `Intersection*_Significant_Overlap.html/.png` (significant counts and overlap)
> Tip: use `Intersection_Native_*` as a high-confidence subset and review Minfi/Sesame outputs for additional signals. The intersection is intentionally conservative.

## Troubleshooting

Start here:
```bash
python3 illumeta.py doctor
```

Common issues:
- **Missing system libraries** (R packages fail to compile): install the OS prerequisites above, then rerun `ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **`xml2` / `libxml-2.0` errors**: install `libxml2-devel` + `zlib` + `pkg-config` (conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge libxml2-devel zlib pkg-config`; Ubuntu/WSL: `sudo apt-get install -y libxml2-dev zlib1g-dev pkg-config`; macOS: `brew install libxml2 zlib pkg-config`), then rerun setup.
- **`lzma.h` / `liblzma` / `Rhtslib` errors**: install xz (conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge xz`; Ubuntu/WSL: `sudo apt-get install -y liblzma-dev`; macOS: `brew install xz`), then rerun setup.
- **`zlib.h` errors**: install zlib (conda: `conda install -c conda-forge zlib`; Ubuntu/WSL: `sudo apt-get install -y zlib1g-dev`; macOS: `brew install zlib`), then rerun setup.
- **`gert` / `git2.h` / `libgit2` errors** (devtools install): install `libgit2` + `pkg-config` (conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge libgit2 pkg-config`; Ubuntu/WSL: `sudo apt-get install -y libgit2-dev pkg-config`; macOS: `brew install libgit2 pkg-config`), then rerun setup.
- **`textshaping` / `ragg` errors** (devtools/tidyverse): install harfbuzz + fribidi + fontconfig + expat (conda: `conda env update -f environment.yml --prune` or `conda install -c conda-forge harfbuzz fribidi fontconfig expat`; Ubuntu/WSL: `sudo apt-get install -y libharfbuzz-dev libfribidi-dev libfontconfig1-dev libexpat1-dev`; macOS: `brew install harfbuzz fribidi fontconfig expat`), then rerun setup.
- **`gfortran` / Fortran errors**: install a Fortran compiler (Ubuntu/WSL: `sudo apt-get install -y gfortran`; macOS: `brew install gcc`; conda: `conda install -c conda-forge gfortran`), then rerun setup.
- **OpenMP / `libomp` errors** (macOS): install `libomp` (`brew install libomp`) and rerun setup.
- **`clang: error: unsupported option '-fopenmp'`** (macOS, often `data.table`): install LLVM and rerun setup with Homebrew clang:
  `brew install llvm && CC=/opt/homebrew/opt/llvm/bin/clang CXX=/opt/homebrew/opt/llvm/bin/clang++ Rscript r_scripts/setup_env.R`.
- **`stringi` configure: `cannot run C++ compiled programs`** (macOS conda): install prebuilt R packages in the conda env, then rerun setup:
  `conda install -c conda-forge r-stringi r-stringr r-tidyr r-plotly r-selectr r-rvest`.
- **`C17 standard requested but CC17 is not defined`**: update to the latest IlluMeta and rerun setup. If it persists, reinstall compilers (macOS: `xcode-select --install`, conda: `conda install -c conda-forge c-compiler cxx-compiler`).
- **Bioconductor version mismatch** (e.g., R 4.4 but Bioc 3.22): rerun with an explicit Bioc version:
  `ILLUMETA_BIOC_VERSION=3.20 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **`lzma decoding result 10` / “read error from connection”** (large annotation downloads corrupted): rerun with retries:
  `ILLUMETA_DOWNLOAD_RETRIES=3 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **GEO download warnings** (`getGEOSuppFiles` returned NULL / 404 for `*_RAW.tar.gz`): some GEO series no longer host `.tar.gz` on FTP. IlluMeta retries and falls back to a direct HTTPS `.tar` download. If it still fails, manually download `GSE*_RAW.tar` from GEO and place it under `projects/<GSE>/<GSE>/`, then rerun `python3 illumeta.py download <GSE> -o projects/<GSE>`.
- **Segfaults or “package built under R x.y” warnings**: you likely switched R versions. Rerun setup and, if needed, clean mismatched packages:
  `ILLUMETA_CLEAN_MISMATCHED_RLIB=1 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **`pandoc: command not found`**: install `pandoc` (Ubuntu: `sudo apt-get install pandoc`, macOS: `brew install pandoc`).
- **Too few samples after QC**: IlluMeta stops if total n is too small for reliable stats; inspect `QC_Summary.csv` and consider adjusting `--qc-intensity-threshold` (or disable by setting `--qc-intensity-threshold 0`).
- **ComBat covariate confounding** (`At least one covariate is confounded with batch`): IlluMeta now auto-drops batch-confounded covariates for ComBat and falls back to limma/none if needed. Check `*_BatchMethodComparison.csv` and `*_Metrics.csv` to see the applied method.
- **Model matrix errors** (`contrasts can be applied only to factors with 2 or more levels`): IlluMeta drops single-level covariates after NA filtering and during batch evaluation. Check `decision_ledger.tsv` for dropped covariates; remove or merge constant columns in `configure.tsv` if the issue persists.
- **`eBayes` failures** (`No finite residual standard deviations`): this can happen with very small n or near-zero variance after correction. IlluMeta skips stability scoring in those cases; consider reducing covariates or disabling batch correction for tiny cohorts.
- **Mixed array sizes**: by default, IlluMeta drops samples that deviate from the modal array size; use `--force-idat` only when appropriate.
- **Missing IDAT pairs**: IlluMeta now auto-filters samples with incomplete `_Grn/_Red` pairs during preflight and logs the decision in `preflight_report.json`. If you prefer to fail instead, use `--keep-missing-idat`.
- **Reference package unavailable** (e.g., FlowSorted.* not in your Bioconductor): IlluMeta falls back to RefFreeEWAS; consider `--cell-reference` or upgrading R/Bioconductor.
- **Sesame: `No normalization control probes found!`**: some EPIC IDATs do not include normalization controls or the cache is stale. IlluMeta now attempts a one-time refresh of `EPIC.1.SigDF` and retries automatically; you can also run `R -q -e 'library(sesame); sesameDataCache("EPIC.1.SigDF")'` manually. If it persists, IlluMeta continues without dye bias correction (or use `--skip-sesame`).
- **CSV configure file**: IlluMeta requires `configure.tsv` to be **tab-delimited**. Convert CSV to TSV and retry.
- **Auto-group failed**: IlluMeta could not find a reliable grouping column; pass `--group-column` or `--group-key`, or edit `configure.tsv` manually.

## Citation
See `CITATION.cff`.

## License
Apache-2.0 (see `LICENSE`).
