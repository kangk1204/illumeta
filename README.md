# IlluMeta

IlluMeta (Illuminating Methylation Analytics) is an end-to-end pipeline that transforms a GEO accession (or raw IDAT files) into publication-ready DNA methylation results and an interactive HTML dashboard.
In addition to array-based DNA methylation analysis, IlluMeta is designed to support sequencing-based methylation data, including methyl-capture sequencing and whole-genome bisulfite sequencing (WGBS), enabling a unified and extensible framework for comprehensive DNA methylation analysis across multiple experimental platforms.

## Highlights
- **One-command workflow**: `download` → edit `configure.tsv` → `analysis`.
- **Two independent pipelines**: **minfi (Noob)** and **sesame** run side-by-side.
- **Consensus (intersection) call set**: CpGs significant in **both** pipelines with the **same direction**.
- **Batch handling**: evaluates correction strategies (SVA/ComBat/limma) when a batch factor exists.
- **Paper-ready artifacts**: interactive HTML + static PNG figures, plus `methods.md`, `analysis_parameters.json`, `sessionInfo.txt`, and `code_version.txt`.

## Installation

IlluMeta is easiest to install with **conda**. If you are on **Windows**, we strongly recommend using **WSL2 (Ubuntu)** for the most reliable R package installation.

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

Recommended (default conda-forge; `setup_env.sh` is not included in this repo):
```bash
conda env create -f environment.yml
conda activate illumeta
```
If conda solve fails, try:
```bash
conda env create -f environment.yml --solver=classic
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
- `ILLUMETA_INSTALL_DEVTOOLS=1` (devtools/tidyverse; requires `libgit2` + `pkg-config` and harfbuzz/fribidi for textshaping/ragg; conda: `conda install -c conda-forge libgit2 harfbuzz fribidi`)
- `ILLUMETA_INSTALL_CLOCKS=1` (methylclock/planet/wateRmelon)
- `ILLUMETA_CLEAN_MISMATCHED_RLIB=1` (use after switching R versions)
- `ILLUMETA_DOWNLOAD_RETRIES=2` (retry big downloads if they fail; increase to 3+ on flaky networks)

If you see errors like `libxml-2.0`, `xml2`, `lzma.h`, or `zlib.h` not found:
- Conda (Linux/WSL): `conda install -c conda-forge libxml2-devel zlib xz pkg-config`
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

Recommended (default conda-forge; `setup_env.sh` is not included in this repo):
```bash
conda env create -f environment.yml
conda activate illumeta
```
If conda solve fails, try:
```bash
conda env create -f environment.yml --solver=classic
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
brew install cmake git pandoc pkg-config openssl@3 libxml2 freetype libpng libtiff jpeg webp harfbuzz fribidi libgit2 libomp gcc
```

###### Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

#### 3) Install R dependencies (first run / CI / reproducible setup)
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
- `ILLUMETA_INSTALL_DEVTOOLS=1` (devtools/tidyverse; requires `libgit2` + `pkg-config` and harfbuzz/fribidi for textshaping/ragg; conda: `conda install -c conda-forge libgit2 harfbuzz fribidi`)
- `ILLUMETA_INSTALL_CLOCKS=1` (methylclock/planet/wateRmelon)
- `ILLUMETA_CLEAN_MISMATCHED_RLIB=1` (use after switching R versions)
- `ILLUMETA_DOWNLOAD_RETRIES=2` (retry big downloads if they fail; increase to 3+ on flaky networks)

If you see errors like `libxml-2.0`, `xml2`, `lzma.h`, or `zlib.h` not found:
- Conda (macOS): `conda install -c conda-forge libxml2 zlib xz pkg-config`
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

Recommended (default conda-forge; `setup_env.sh` is not included in this repo):
```bash
conda env create -f environment.yml
conda activate illumeta
```
If conda solve fails, try:
```bash
conda env create -f environment.yml --solver=classic
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
- `ILLUMETA_INSTALL_DEVTOOLS=1` (devtools/tidyverse; requires `libgit2` + `pkg-config` and harfbuzz/fribidi for textshaping/ragg; conda: `conda install -c conda-forge libgit2 harfbuzz fribidi`)
- `ILLUMETA_INSTALL_CLOCKS=1` (methylclock/planet/wateRmelon)
- `ILLUMETA_CLEAN_MISMATCHED_RLIB=1` (use after switching R versions)
- `ILLUMETA_DOWNLOAD_RETRIES=2` (retry big downloads if they fail; increase to 3+ on flaky networks)

If you see errors like `libxml-2.0`, `xml2`, `lzma.h`, or `zlib.h` not found:
- Conda (Linux/WSL): `conda install -c conda-forge libxml2-devel zlib xz pkg-config`
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

## Quick start

### 1) Download a dataset from GEO
```bash
# Activate your environment (choose one)
# - Conda: conda activate illumeta
# - venv:  source .venv/bin/activate

python3 illumeta.py download GSE121633 -o projects/GSE121633
```

### 2) Assign groups
Edit `projects/GSE121633/configure.tsv` and fill in `primary_group` (e.g., `Control` / `Case`).

### 3) Run analysis
```bash
python3 illumeta.py analysis \
  -i projects/GSE121633 \
  --group_con Control \
  --group_test Case
```
Note: the default output folder name is derived from the group labels. If it contains non-ASCII characters, IlluMeta normalizes it to a safe ASCII name for filesystem compatibility (the dashboard filename follows the folder name).

### 4) Open the dashboard
Open the generated HTML:
`projects/GSE121633/Case_vs_Control_results_index.html`

## Usage

### Search GEO for IDAT-enabled datasets
```bash
python3 illusearch.py --keywords "breast cancer" --email your_email@example.com -o search_results.tsv
```

### Analyze your own IDATs (non-GEO)
1. Create `my_project/idat/` and place `_Grn.idat` / `_Red.idat` pairs there.
2. Create `my_project/configure.tsv` with at least:
   - `Basename` (e.g., `idat/Sample1_R01C01`)
   - `primary_group` (e.g., `Treated`, `Untreated`)
3. Run:
```bash
python3 illumeta.py analysis -i my_project --group_con Untreated --group_test Treated
```

### Common analysis options
```bash
# Tighten thresholds
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --pval 0.01 --lfc 1.0

# Disable SVA (very small n)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --disable-sva

# Provide covariates to always try (if present in configure.tsv)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --include-covariates age,sex

# Use epigenetic clocks as candidate covariates (auto-selected if relevant)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --include-clock-covariates

# Placenta clocks (planet)
python3 illumeta.py analysis -i projects/GSE307314 --group_con control --group_test test --tissue Placenta

# Mixed-array safeguard override (only if you know what you're doing)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --force-idat
```

## Outputs (what to use in a paper)

IlluMeta writes results to the analysis output directory (default: `[input]/[test]_vs_[con]_results/`).

### Reproducibility
- `methods.md`: auto-generated methods summary for the run.
- `analysis_parameters.json`: thresholds and key settings.
- `sessionInfo.txt`: full R session and package versions.
- `code_version.txt`: git commit hash (when available).

### QC
- `QC_Summary.csv`: sample/probe QC counts.
- `Sample_QC_Metrics.csv`: per-sample QC metrics.
- `Sample_QC_DetectionP_FailFraction.png` and `Sample_QC_Intensity_Medians.png`: static QC figures (HTML versions are also saved).

### Differential methylation (per pipeline)
For each pipeline (`Minfi`, `Sesame`):
- `*_DMPs_full.csv`: full differential methylation results table.
- `*_Volcano.html/.png`, `*_Manhattan.html/.png`, `*_QQPlot.html/.png`
- `*_Top100_Heatmap.html/.png`
- `*_DMR_Volcano.html/.png`, `*_DMR_Manhattan.html/.png`, `*_Top_DMRs_Heatmap.html/.png`

### Consensus (intersection) call set
- `Intersection_Consensus_DMPs.csv` and `Intersection_Consensus_DMPs.html`
- `Intersection_LogFC_Concordance.html/.png` (minfi vs sesame logFC concordance)
- `Intersection_Significant_Overlap.html/.png` (significant counts and overlap)

## Troubleshooting

Start here:
```bash
python3 illumeta.py doctor
```

Common issues:
- **Missing system libraries** (R packages fail to compile): install the OS prerequisites above, then rerun `ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **`xml2` / `libxml-2.0` errors**: install `libxml2` + `pkg-config` (conda Linux/WSL: `conda install -c conda-forge libxml2-devel pkg-config`; conda macOS: `conda install -c conda-forge libxml2 pkg-config`; Ubuntu/WSL: `sudo apt-get install -y libxml2-dev pkg-config`; macOS: `brew install libxml2 pkg-config`), then rerun setup.
- **`lzma.h` / `liblzma` / `Rhtslib` errors**: install xz (conda: `conda install -c conda-forge xz`; Ubuntu/WSL: `sudo apt-get install -y liblzma-dev`; macOS: `brew install xz`), then rerun setup.
- **`zlib.h` errors**: install zlib (conda: `conda install -c conda-forge zlib`; Ubuntu/WSL: `sudo apt-get install -y zlib1g-dev`; macOS: `brew install zlib`), then rerun setup.
- **`gert` / `git2.h` / `libgit2` errors** (devtools install): install `libgit2` + `pkg-config` (conda: `conda install -c conda-forge libgit2 pkg-config`; Ubuntu/WSL: `sudo apt-get install -y libgit2-dev pkg-config`; macOS: `brew install libgit2 pkg-config`), then rerun setup.
- **`textshaping` / `ragg` errors** (devtools/tidyverse): install harfbuzz + fribidi (conda: `conda install -c conda-forge harfbuzz fribidi`; Ubuntu/WSL: `sudo apt-get install -y libharfbuzz-dev libfribidi-dev`; macOS: `brew install harfbuzz fribidi`), then rerun setup.
- **`gfortran` / Fortran errors**: install a Fortran compiler (Ubuntu/WSL: `sudo apt-get install -y gfortran`; macOS: `brew install gcc`; conda: `conda install -c conda-forge gfortran`), then rerun setup.
- **OpenMP / `libomp` errors** (macOS): install `libomp` (`brew install libomp`) and rerun setup.
- **`C17 standard requested but CC17 is not defined`**: update to the latest IlluMeta and rerun setup. If it persists, reinstall compilers (macOS: `xcode-select --install`, conda: `conda install -c conda-forge c-compiler cxx-compiler`).
- **Bioconductor version mismatch** (e.g., R 4.4 but Bioc 3.22): rerun with an explicit Bioc version:
  `ILLUMETA_BIOC_VERSION=3.20 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **`lzma decoding result 10` / “read error from connection”** (large annotation downloads corrupted): rerun with retries:
  `ILLUMETA_DOWNLOAD_RETRIES=3 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **Segfaults or “package built under R x.y” warnings**: you likely switched R versions. Rerun setup and, if needed, clean mismatched packages:
  `ILLUMETA_CLEAN_MISMATCHED_RLIB=1 ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`.
- **`pandoc: command not found`**: install `pandoc` (Ubuntu: `sudo apt-get install pandoc`, macOS: `brew install pandoc`).
- **Too few samples after QC**: IlluMeta stops if total n is too small for reliable stats; inspect `QC_Summary.csv` and consider adjusting `--qc-intensity-threshold` (or disable by setting `--qc-intensity-threshold 0`).
- **Mixed array sizes**: by default, IlluMeta drops samples that deviate from the modal array size; use `--force-idat` only when appropriate.

## Citation
See `CITATION.cff`.

## License
Apache-2.0 (see `LICENSE`).
