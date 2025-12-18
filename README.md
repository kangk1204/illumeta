# IlluMeta

IlluMeta (Illuminating Methylation Analytics) is an end-to-end pipeline to go from a GEO accession (or raw IDATs) to publication-ready DNA methylation results and an interactive HTML dashboard.

## Highlights
- **One-command workflow**: `download` → edit `configure.tsv` → `analysis`.
- **Two independent pipelines**: **minfi (Noob)** and **sesame** run side-by-side.
- **Consensus (intersection) call set**: CpGs significant in **both** pipelines with the **same direction**.
- **Batch handling**: evaluates correction strategies (SVA/ComBat/limma) when a batch factor exists.
- **Paper-ready artifacts**: interactive HTML + static PNG figures, plus `methods.md`, `analysis_parameters.json`, `sessionInfo.txt`, and `code_version.txt`.

## Installation

### 0) Install prerequisites (Git + Conda)
IlluMeta is easiest to install with **conda**. If you are on **Windows**, we strongly recommend using **WSL2 (Ubuntu)** for the most reliable R package installation.

#### Windows 11 (recommended: WSL2 + Ubuntu)
In **PowerShell (Admin)**:
```powershell
wsl --install -d Ubuntu
```

Then open **Ubuntu** and install Git + curl:
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

#### Git (if not installed yet)
- Ubuntu/WSL: `sudo apt-get install -y git`
- macOS: `xcode-select --install` (or `brew install git`)
- Windows (native): https://git-scm.com/download/win

#### Conda (if not installed yet)
- Recommended: **Miniforge** (conda-forge): https://github.com/conda-forge/miniforge
  - Linux/WSL (x86_64): `Miniforge3-Linux-x86_64.sh`
  - macOS (Apple Silicon): `Miniforge3-MacOSX-arm64.sh`
  - macOS (Intel): `Miniforge3-MacOSX-x86_64.sh`
  - Windows: `Miniforge3-Windows-x86_64.exe`
- After installation, open a new terminal and confirm: `conda --version`

### 1) Clone the repository
```bash
git clone https://github.com/kangk1204/illumeta.git
cd illumeta
```

### 2) Create an environment (choose one)

#### Option A: Conda (recommended for beginners)
This uses conda to provide R/Python plus the system libraries needed by many R packages.
```bash
conda env create -f environment.yml
conda activate illumeta
```

#### Option B: venv + system R
##### Prerequisites
- **Python** 3.8+
- **R** 4.2+
- **pandoc** (required for self-contained HTML reports via `htmlwidgets`)
- **Git**
- **System libraries** (for compiling common R packages)

Ubuntu/WSL:
```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake git pandoc \
  libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev
```

macOS (Homebrew):
```bash
brew install cmake git pandoc openssl@3 libxml2 freetype libpng libtiff jpeg webp
```

##### Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### 3) Install R dependencies (first run / CI / reproducible setup)
This installs required R/Bioconductor packages into the repo-local library (`.r-lib/`) and may take ~10–30 minutes.
```bash
# Ensure your environment is activated first:
# - Conda: conda activate illumeta
# - venv:  source .venv/bin/activate

ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
```

### 4) Check your environment (recommended)
```bash
python3 illumeta.py doctor
```

Notes on R libraries:
- IlluMeta **defaults to a repo-local** `R_LIBS_USER=.r-lib/` for reproducibility.
- To keep an externally-set `R_LIBS_USER`, run with `ILLUMETA_RESPECT_R_LIBS_USER=1`.

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
- **`pandoc: command not found`**: install `pandoc` (Ubuntu: `sudo apt-get install pandoc`, macOS: `brew install pandoc`).
- **Too few samples after QC**: IlluMeta stops if total n is too small for reliable stats; inspect `QC_Summary.csv` and consider adjusting `--qc-intensity-threshold` (or disable by setting `--qc-intensity-threshold 0`).
- **Mixed array sizes**: by default, IlluMeta drops samples that deviate from the modal array size; use `--force-idat` only when appropriate.

## Citation
See `CITATION.cff`.

## License
Apache-2.0 (see `LICENSE`).
