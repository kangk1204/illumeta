# ✨ IlluMeta: Automated DNA Methylation Analysis Pipeline

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Focus](https://img.shields.io/badge/focus-Epigenetics-green.svg)]()

> **"From IDATs to a publication-ready dashboard in one command."**
> IlluMeta transforms a GEO accession (or raw IDATs) into reproducible methylation results plus an interactive HTML report.

---

## 🎯 What is IlluMeta?

**IlluMeta** is a user-friendly tool for analyzing DNA methylation data from Illumina arrays (450K, EPIC, EPIC v2).

> **Note:** IlluMeta currently supports **human (Homo sapiens)** Infinium methylation arrays only. Non-human arrays are not supported yet.

**In simple terms:** DNA methylation is a chemical modification that can turn genes "on" or "off". Scientists study methylation differences between groups (e.g., healthy vs. disease) to understand biological processes and discover disease biomarkers.

### Who is this for?

| You are... | IlluMeta helps you... |
|------------|----------------------|
| 🔬 **Biologist** with methylation data | Get publication-ready results without coding |
| 🎓 **Graduate student** learning EWAS | Understand the analysis workflow step-by-step |
| 👨‍💻 **Bioinformatician** | Save time with automated dual-pipeline analysis |
| 📊 **Researcher** writing a paper | Generate reproducible figures and methods |

### How does it work?

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  1. INPUT   │ ──► │  2. PROCESS │ ──► │ 3. ANALYZE  │ ──► │  4. OUTPUT  │
│             │     │             │     │             │     │             │
│ • GEO ID or │     │ • Quality   │     │ • Find DMPs │     │ • Dashboard │
│ • IDAT files│     │   control   │     │   (changed  │     │ • Tables    │
│ • Sample    │     │ • Normalize │     │   sites)    │     │ • Figures   │
│   groups    │     │   data      │     │ • Find DMRs │     │ • Reports   │
│             │     │ • Batch     │     │   (regions) │     │             │
│             │     │   correct   │     │ • Consensus │     │             │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
```

---

## 📚 Key Terms (Glossary)

New to DNA methylation analysis? Here are the key terms you'll encounter:

| Term | What it means | Simple analogy |
|------|---------------|----------------|
| **IDAT** | Raw data files from Illumina arrays | Like a "RAW" photo file |
| **CpG** | A location in DNA where methylation occurs | An address on a street |
| **Beta value** | Methylation level (0-1, where 1 = fully methylated) | Dimmer switch (0=off, 1=full) |
| **DMP** | Differentially Methylated Position (single CpG) | One changed address |
| **DMR** | Differentially Methylated Region (multiple CpGs) | A changed neighborhood |
| **FDR** | False Discovery Rate (adjusted p-value) | Confidence after checking many sites |
| **Batch effect** | Technical variation between experiments | Different lighting in photos |
| **Normalization** | Removing technical variation | Color-correcting photos |

<details>
<summary><strong>📖 More terms (click to expand)</strong></summary>

| Term | What it means |
|------|---------------|
| **GEO** | Gene Expression Omnibus - public database for genomic data |
| **Minfi** | R package for methylation analysis (method 1) |
| **Sesame** | R package for methylation analysis (method 2) |
| **Consensus** | CpGs found significant by BOTH methods (high confidence) |
| **Covariate** | Variables like age, sex that might affect results |
| **SVA** | Surrogate Variable Analysis - finds hidden batch effects |
| **ComBat** | Method to remove known batch effects |
| **Lambda (λ)** | Genomic inflation factor - checks for bias (ideal ≈ 1.0) |

</details>

---

## 🚀 Quick Start (Copy & Paste)

> **Prerequisites:** You need `conda` (recommended) or `mamba`. If you don't have it yet, install **Miniforge (conda)** first (macOS/Ubuntu) using the snippet below.

<details>
<summary><strong>Install conda (Miniforge) - macOS / Ubuntu (fresh machine)</strong></summary>

Ubuntu / WSL:
```bash
sudo apt-get update && sudo apt-get install -y curl git
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source ~/.bashrc
conda --version
```
If your Ubuntu is ARM64, use `Miniforge3-Linux-aarch64.sh` instead.

macOS:
```bash
xcode-select --install
```
Apple Silicon (arm64):
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```
Intel (x86_64):
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```
Tip: run `uname -m` to confirm your Mac architecture (`arm64` vs `x86_64`).
If you use bash on macOS, replace `zsh` with `bash` and source `~/.bashrc` instead.
</details>

### Option A: One-line installer (Recommended for beginners)
```bash
git clone https://github.com/kangk1204/illumeta.git && cd illumeta && ./scripts/install_full.sh
```
This takes **30-60 minutes** (downloads R packages). Grab a coffee! ☕

### Option B: Step-by-step
```bash
# 1) Clone and create conda environment
git clone https://github.com/kangk1204/illumeta.git
cd illumeta
conda env create -f environment.yml
conda activate illumeta

# 2) Install R packages (takes 10-30 min)
Rscript r_scripts/setup_env.R

# 3) Verify installation
python3 illumeta.py doctor
```

### Run an example
```bash
# Download a GEO dataset (DCIS vs Adjacent Normal breast tissue, 450K)
python3 illumeta.py download GSE66313 -o projects/GSE66313

# Run analysis with auto-grouping
python3 illumeta.py analysis -i projects/GSE66313 \
  --group_con Control --group_test Case \
  --auto-group --group-column source_name_ch1 \
  --group-map "Adjacent-Normal=Control;DCIS=Case" \
  --tier3-on-fail skip
```

Open the dashboard:
`projects/GSE66313/Case_vs_Control_results_index.html`

---

## 📊 What You Get

<!-- Dashboard preview screenshot will be added before publication -->
<!-- ![Dashboard preview](docs/assets/dashboard_preview.gif) -->

*Interactive dashboard generated automatically. No additional coding required.*

---

## 🚑 If Something Breaks

Run the self-check doctor:

```bash
python3 illumeta.py doctor
```

---

## 📦 OS-Specific Install (Collapsed)

<details>
<summary><strong>Windows (WSL2 required)</strong></summary>

1. Install WSL2 (Ubuntu).
2. Follow the **Ubuntu** section in Installation below.
</details>

<details>
<summary><strong>Mac (Apple Silicon)</strong></summary>

1. Install Miniforge (conda).
2. Follow the **Mac** section in Installation below.
</details>

<details>
<summary><strong>Linux</strong></summary>

1. Install Git + Conda.
2. Follow the **Ubuntu** section in Installation below.
</details>

---

## ✨ Why IlluMeta?

### For Beginners
- 🎯 **Simple workflow**: Just 3 commands from data to results
- 🔍 **Automatic quality control**: Catches problems before they affect your analysis
- 📊 **Interactive dashboard**: Explore results without coding
- 📝 **Auto-generated methods**: Ready-to-use text for your paper

### For Experts
- **Dual-pipeline design**: Runs both **Minfi** and **Sesame** independently
- **Consensus calling**: High-confidence results where both methods agree
- **Adaptive batch correction**: Automatically selects optimal method (SVA/ComBat/limma)
- **CRF robustness framework**: Sample-size-adaptive quality assessment
- **Full reproducibility**: All parameters and decisions logged

<details>
<summary><strong>🔬 Technical highlights (click to expand)</strong></summary>

- **Two independent pipelines**: Minfi (Noob) and Sesame run side-by-side; Sesame reports both strict (Minfi-aligned) and native (pOOBAH-preserving) views.
- **High-confidence consensus**: CpGs significant in BOTH pipelines with the SAME direction.
- **Batch handling**: Evaluates correction strategies (SVA/ComBat/limma) when a batch factor exists.
- **CRF**: Sample-size-adaptive robustness report (MMC/NCS/SSS) with tiered warnings.
- **Defensive stats**: Guards against low-variance or single-group covariates.
- **Executive dashboard**: Verdict + CRF quick stats + warnings at the top.
- **Paper-ready outputs**: HTML + PNG figures, methods.md, summary.json, sessionInfo.txt.

</details>

## Analysis Pipeline Overview

IlluMeta runs a fully automated dual-pipeline analysis. The diagram below shows the complete flow from raw IDAT files to the interactive dashboard:

```
Input: IDAT files + sample sheet
        │
        ├── Preflight ── validate samples, detect array (450k / EPIC / EPIC v2)
        │
        ├── Sample QC ── median intensity check, sex mismatch detection
        │
        ├── Cell Composition ── reference-based (blood/saliva) or reference-free (RefFreeEWAS)
        │
        ├── Probe QC ── detection P, cross-reactive, SNP (MAF), sex chromosome removal
        │
        │       ┌──────────────────────┬────────────────────────┐
        │       │                      │                        │
        ▼       ▼                      ▼                        ▼
    ┌────────────────┐   ┌──────────────────┐   ┌───────────────────────┐
    │     Minfi      │   │  SeSAMe (Strict) │   │   SeSAMe (Native)     │
    │  Noob (minfi)  │   │  Same probe set  │   │  pOOBAH-preserved     │
    │                │   │  as Minfi        │   │  probe set (wider)    │
    └───────┬────────┘   └────────┬─────────┘   └──────────┬────────────┘
            │                     │                        │
            ├── Auto-covariate selection (PCA + PVCA)      │
            ├── Batch evaluation (SVA / ComBat / limma)    │
            ├── DMP analysis (limma) + lambda QC           │
            ├── DMR analysis (dmrff)                       │
            ├── Permutation testing (FPR)                  │
            ├── variancePartition                          │
            └── Epigenetic clocks (methylclock)            │
                     │                                     │
                     ▼                                     │
        ┌─────────────────────────┐                        │
        │  Consensus Intersection │◄───────────────────────┘
        │  (Strict & Native)      │
        │  Same direction + both  │
        │  pipelines significant  │
        └────────┬────────────────┘
                 │
                 ▼
        ┌───────────────────────────────────────┐
        │  CRF Robustness Assessment            │
        │  + CAF Score + Signal Preservation    │
        │  + summary.json + decision_ledger.tsv │
        └────────┬──────────────────────────────┘
                 │
                 ▼
        Interactive HTML Dashboard + methods.md
```

**Key design principle**: Minfi and SeSAMe use fundamentally different preprocessing algorithms. By requiring significance in **both** pipelines (consensus intersection), IlluMeta minimizes false positives and increases reproducibility.

## Dashboard Navigation

IlluMeta generates an interactive HTML dashboard as the primary interface for exploring results. The dashboard includes:

| Section | What it shows |
|---------|---------------|
| **Executive Summary** | Verdict badge, DMP counts, key metrics (lambda, stability, CRF tier) |
| **Beginner Path** | 3-step guide: (1) Check QC, (2) Review consensus intersection, (3) Explore pipelines |
| **Run Controls & QC** | Analysis parameters, sample/probe QC summary |
| **Pipeline Tabs** | 5 tabs: Intersection Native, Intersection Strict, Minfi, Sesame Strict, Sesame Native |
| **Run Documentation** | Methods text, decision ledger, logs, and parameter files for reproducibility |

Each pipeline tab contains expandable sections for Primary Results (volcano, Manhattan, heatmap, DMP table), Region-Level Results (DMRs), Quality & Diagnostics (PCA, Q-Q, PVCA), Batch & Covariates, and Clocks & Age.

For a complete dashboard guide with interpretation instructions, see the **Supplementary Data** document (generate via `python3 scripts/build_supplementary_data_docx.py`; requires `./scripts/install_full.sh --paper`).

## Supplementary materials (recommended for manuscripts)
For transparency and reproducibility, we recommend including the following files in your Supplementary Materials:
1. `preflight_report.json` - input data validation record
2. `analysis_parameters.json` - full analysis parameters
3. `decision_ledger.tsv` - automated decision log
4. `CRF_Sample_Tier.csv` - sample size tier assessment
5. `methods.md` - auto-generated Methods text

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

> **Most users:** Just run the [Quick Start](#-quick-start-copy--paste) above. This section is for troubleshooting or advanced setup.

IlluMeta is easiest to install with **conda**. Windows users should use **WSL2 (Ubuntu)**.

### Quick full install (all OS)

Prerequisites (fresh machine):
- `conda` or `mamba` in `PATH` (Miniforge recommended).
- macOS: `xcode-select --install` (includes git + macOS SDK headers)
- Ubuntu/WSL: `sudo apt-get update && sudo apt-get install -y git curl`

<details>
<summary><strong>❓ No conda yet? Install Miniforge (recommended)</strong></summary>

Linux (x86_64):
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source ~/.bashrc
conda --version
```
If your Linux is ARM64, use `Miniforge3-Linux-aarch64.sh` instead.

macOS:
Apple Silicon (arm64):
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```

Intel (x86_64):
```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```
Tip: run `uname -m` to confirm your Mac architecture (`arm64` vs `x86_64`).
If you use bash on macOS, replace `zsh` with `bash` and source `~/.bashrc` instead.
</details>

```bash
git clone https://github.com/kangk1204/illumeta.git && cd illumeta && ./scripts/install_full.sh
```

If you also want to build the Application Note DOCX / paper figures on the same machine:
```bash
./scripts/install_full.sh --paper
```

Optional (manuscript assets; requires local benchmark outputs):
```bash
# Figures/tables for the Application Note
python3 scripts/generate_application_note_figures.py \
  --results-dir benchmarks/application_note \
  --out-dir benchmarks/paper_figures

# Build DOCX drafts (NOT required to run IlluMeta itself)
python3 scripts/build_application_note_docx.py \
  --results-dir benchmarks/application_note \
  --figures-dir benchmarks/paper_figures \
  --out benchmarks/application_note/IlluMeta_Application_Note.docx \
  --contact-email you@example.com

python3 scripts/build_supplementary_data_docx.py \
  --results-dir benchmarks/application_note \
  --figures-dir benchmarks/paper_figures \
  --out benchmarks/application_note/IlluMeta_Supplementary_Data.docx
```

<details>
<summary><strong>📋 install_full.sh options</strong></summary>

- `--r45` uses `environment-r45.yml` (R 4.5).
- `--env-file PATH` uses a custom conda env file.
- `--env NAME` overrides the env name.
- `--skip-doctor` skips the final `illumeta.py doctor` check.

Logs are saved to `projects/illumeta_install_full_*.log`.
</details>

<details>
<summary><strong>⚠️ Important notes before installing</strong></summary>

- **Pick one environment**: conda **or** system R. Do **not** mix them.
- If using conda, always run `conda activate illumeta` before `Rscript`.
- Quick check (should point inside your chosen environment):
```bash
which R
R -q -e 'cat(R.version.string, "\n"); cat(R.home(), "\n")'
```
- If you previously installed with a different R version, the installer auto-cleans mismatched packages.
</details>

---

<details>
<summary><strong>🐧 Ubuntu / WSL2 - Detailed Installation</strong></summary>

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

This installs required R/Bioconductor packages into the repo-local library (`.r-lib/R-<major.minor>`).

> **Expected time:** 10-30 minutes (core), 30-60 minutes (full install with clocks/devtools). The terminal may appear "stuck" while compiling - this is normal.

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

</details>

<details>
<summary><strong>🍎 macOS (Apple Silicon & Intel) - Detailed Installation</strong></summary>

### macOS (Apple Silicon M1-M4 and Intel)

#### 0) Install prerequisites (Git + Conda)

**Step 1: Install Xcode Command Line Tools** (includes git and a macOS SDK):
```bash
xcode-select --install
```
> A popup will appear. Click **"Install"** and wait for completion (may take 5-10 minutes).

**Step 2: Install Miniforge (conda)**:

<details>
<summary><strong>Apple Silicon (M1/M2/M3/M4)</strong></summary>

```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```
</details>

<details>
<summary><strong>Intel Mac (x86_64)</strong></summary>

```bash
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
source ~/.zshrc
conda --version
```
</details>

> **Tip:** Not sure which Mac you have? Run `uname -m`. If it says `arm64`, you have Apple Silicon. If it says `x86_64`, you have Intel.

If you use bash instead of zsh, replace `zsh` with `bash` and source `~/.bashrc` instead.

**Optional (only if you use system R or install `devtools`/`roxygen2` outside conda): Homebrew + system libraries**
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install cmake git pandoc pkg-config openssl@3 libxml2 freetype libpng libtiff jpeg webp harfbuzz fribidi fontconfig libgit2 libomp gcc
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

This installs required R/Bioconductor packages into the repo-local library (`.r-lib/R-<major.minor>`).

> **Expected time:** 10-30 minutes (core), 30-60 minutes (full install with clocks/devtools). The terminal may appear "stuck" while compiling - this is normal.

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

</details>

<details>
<summary><strong>🪟 Windows 11 (WSL2) - Detailed Installation</strong></summary>

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

This installs required R/Bioconductor packages into the repo-local library (`.r-lib/R-<major.minor>`).

> **Expected time:** 10-30 minutes (core), 30-60 minutes (full install with clocks/devtools). The terminal may appear "stuck" while compiling - this is normal.

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

</details>

---

### Notes on R libraries
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
  -i projects/GSE66313 \
  --group_con Control \
  --group_test Case \
  --auto-group --group-column source_name_ch1 \
  --group-map "Adjacent-Normal=Control;DCIS=Case" \
  --tier3-on-fail skip
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

python3 illumeta.py download GSE66313 -o projects/GSE66313
# If a GEO series has multiple platforms, force one by GPL ID:
python3 illumeta.py download GSE66313 -o projects/GSE66313 --platform GPL13534
```

### 2) Assign groups (manual or auto)
Option A (manual):
Edit `projects/GSE66313/configure.tsv` and fill in `primary_group` (e.g., `Control` / `Case`).
`configure.tsv` must be **tab-delimited (TSV)**; CSV is not supported.

Option B (auto-group on analysis):
If your dataset's `primary_group` is empty, IlluMeta can populate it from metadata:
```bash
python3 illumeta.py analysis \
  -i projects/GSE12345 \
  --group_con Control \
  --group_test Case \
  --auto-group \
  --group-column <your_column_name>
```
If your grouping is encoded in GEO characteristics, you can use `--group-key` (e.g., `--group-key disease`).

### 3) Run analysis
```bash
python3 illumeta.py analysis \
  -i projects/GSE66313 \
  --group_con Control \
  --group_test Case \
  --auto-group --group-column source_name_ch1 \
  --group-map "Adjacent-Normal=Control;DCIS=Case" \
  --tier3-on-fail skip
```
If you filled `primary_group` manually (Option A), you can omit the `--auto-group` flags.
Note: the default output folder name is derived from the group labels. If it contains non-ASCII characters, IlluMeta normalizes it to a safe ASCII name for filesystem compatibility (the dashboard filename follows the folder name).
Optional (signal preservation checks):
```bash
# Provide a CpG marker list (TSV/CSV with CpG column or one CpG per line)
python3 illumeta.py analysis \
  -i projects/GSE66313 \
  --group_con Control \
  --group_test Case \
  --auto-group --group-column source_name_ch1 \
  --group-map "Adjacent-Normal=Control;DCIS=Case" \
  --tier3-on-fail skip \
  --marker-list markers.tsv
```
This generates `*_Signal_Preservation.csv` and (if provided) `*_Known_Marker_Summary.csv`.

### 4) Open the dashboard
Open the generated HTML:
`projects/GSE66313/Case_vs_Control_results_index.html`

### 5) Interpret results (beginner checklist)

**🎯 Quick interpretation guide:**

| Step | What to check | Where to look | What's good? |
|------|---------------|---------------|--------------|
| 1️⃣ | **Quality Control** | `QC_Summary.csv` | >90% samples pass |
| 2️⃣ | **Batch Effects** | Dashboard "Batch" section | Reduced after correction |
| 3️⃣ | **Lambda (inflation)** | Dashboard header | λ ≈ 0.9 - 1.1 |
| 4️⃣ | **DMP counts** | Dashboard summary | Depends on your study |
| 5️⃣ | **Top hits** | `*_Top_DMPs.html` | Biologically relevant genes |

**📂 Key output files explained:**

| File | What it contains | When to use it |
|------|------------------|----------------|
| `*_index.html` | Interactive dashboard | First look at results |
| `Intersection_Consensus_DMPs.csv` | High-confidence DMPs | Primary results for paper |
| `*_Volcano.html` | Effect size vs significance plot | Visualize overall signal |
| `*_Manhattan.html` | Genomic distribution of hits | Check for regional enrichment |
| `methods.md` | Auto-generated methods text | Copy to your paper |
| `decision_ledger.tsv` | All automated decisions | Reproducibility/debugging |

**⚠️ Red flags to watch for:**
- Lambda (λ) > 1.5 or < 0.8 → possible bias or batch issues
- Many samples failing QC → data quality problems
- Zero consensus DMPs but many pipeline-specific → methods disagree (check data)

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
  --group-map "normal=Control;tumor=Case"
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
  --auto-group --group-column disease_state --group-map "normal=Control;tumor=Case"

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

### Cell Composition Estimation

IlluMeta estimates cell type composition to adjust for cellular heterogeneity, which is critical for accurate differential methylation analysis.

#### Supported Methods

| `--tissue` | Method | Cell Types | Best For |
|------------|--------|------------|----------|
| `Auto` | **RefFreeEWAS** | Latent1-5 (data-driven) | Unknown tissue, tumors, mixed samples |
| `Blood` | FlowSorted.Blood.EPIC | CD8T, CD4T, NK, Bcell, Mono, Neu | Blood/PBMC samples |
| `CordBlood` | FlowSorted.CordBlood | Nucleated RBC, CD8T, CD4T, etc. | Cord blood samples |
| `Placenta` | planet | Trophoblast, Stromal, etc. | Placenta samples |
| `DLPFC` | FlowSorted.DLPFC | NeuN+, NeuN- | Brain (prefrontal cortex) |

#### Which method should I use?

<details>
<summary><strong>Decision guide (click to expand)</strong></summary>

**Use `--tissue Auto` (Reference-Free) when:**
- You don't know the exact tissue type
- Analyzing tumor samples (disease may alter cell methylation signatures)
- Mixed or heterogeneous samples
- You want a method that works across all tissues

**Use `--tissue Blood` (Reference-Based) when:**
- Analyzing blood or PBMC samples
- You need interpretable cell type names (CD8T, CD4T, NK, etc.)
- Clinical/epidemiological studies where cell types matter

**Pros and Cons:**

| Aspect | Reference-Free (Auto) | Reference-Based (Blood, etc.) |
|--------|----------------------|------------------------------|
| **Interpretability** | Low (Latent1, Latent2...) | High (CD8T, NK, Mono...) |
| **Tissue flexibility** | Any tissue | Tissue-specific only |
| **Disease bias** | None | May be biased if disease alters cell methylation |
| **Precision** | Moderate | High for matched tissues |

</details>

#### Examples
```bash
# Reference-free (works for any tissue)
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --tissue Auto

# Blood samples with reference-based deconvolution
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --tissue Blood

# Placenta samples
python3 illumeta.py analysis -i projects/GSE12345 --group_con Control --group_test Case --tissue Placenta
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

## ❓ Frequently Asked Questions (FAQ)

<details>
<summary><strong>🆕 I'm completely new to methylation analysis. Where do I start?</strong></summary>

1. **Install IlluMeta** using the Quick Start section above
2. **Run the example** with GSE66313 (DCIS vs Normal, 450K, 55 samples)
3. **Open the dashboard** and explore the interactive plots
4. **Read the methods.md** file - it explains what was done in plain English
5. **Check QC first** - look at `QC_Summary.csv` to ensure data quality

</details>

<details>
<summary><strong>📊 Which results should I use - Minfi, Sesame, or Consensus?</strong></summary>

| Use this... | When you want... | Trade-off |
|-------------|------------------|-----------|
| **Consensus** | Highest confidence, conservative | May miss some real signals |
| **Minfi** | Well-established, widely used | Single method perspective |
| **Sesame** | Latest algorithms, masking low-quality probes | Different from traditional minfi |

**Recommendation for papers**: Report Consensus as primary results, mention Minfi/Sesame counts for comparison.

</details>

<details>
<summary><strong>⚠️ I got very few (or zero) significant DMPs. Is something wrong?</strong></summary>

This is common! Check these in order:

1. **Sample size**: Small n (< 20) often yields few significant results
2. **QC**: Look at `QC_Summary.csv` - are many samples/probes failing?
3. **Effect size**: Your groups may have subtle differences
4. **Lambda**: Check the λ value in dashboard - should be near 1.0
5. **Try relaxing thresholds**: Use `--pval 0.1` or `--delta-beta 0` temporarily to see if signals exist

Zero DMPs doesn't mean your experiment failed - it might mean no large methylation differences exist.

</details>

<details>
<summary><strong>🔢 How many samples do I need?</strong></summary>

| Sample size | What to expect |
|-------------|----------------|
| **6-12 total** | Exploratory only; very limited power |
| **12-24 total** | Can detect large effects; validate findings independently |
| **24-50 total** | Reasonable power for moderate effects |
| **50+ total** | Good power; robust statistical analysis |

**Rule of thumb**: Aim for at least 10-15 samples per group for publishable results.

</details>

<details>
<summary><strong>📁 What files should I include in my paper's supplementary?</strong></summary>

**Essential (for reproducibility)**:
- `methods.md` - Auto-generated methods text
- `analysis_parameters.json` - All parameters used
- `decision_ledger.tsv` - All automated decisions

**Recommended (for transparency)**:
- `QC_Summary.csv` - Quality control summary
- `*_Metrics.csv` - Lambda, batch correction details
- `Intersection_Consensus_DMPs.csv` - Full consensus results

</details>

<details>
<summary><strong>🐛 The analysis is taking very long / seems stuck</strong></summary>

Normal analysis times:
| Sample count | Expected time |
|--------------|---------------|
| 20 samples | 10-20 minutes |
| 50 samples | 30-60 minutes |
| 100+ samples | 1-3 hours |

If stuck for >2x expected time:
1. Check memory usage (`htop` or Activity Monitor)
2. Look at the log file in the output directory
3. Try running with `--skip-sesame` to isolate issues

</details>

<details>
<summary><strong>🔄 Can I rerun with different parameters?</strong></summary>

Yes! Just specify a different output folder:

```bash
# Stricter thresholds
python3 illumeta.py analysis -i projects/GSE12345 \
  --group_con Control --group_test Case \
  --pval 0.01 --lfc 1.0 \
  --output projects/GSE12345/strict_results

# Different batch correction
python3 illumeta.py analysis -i projects/GSE12345 \
  --group_con Control --group_test Case \
  --disable-sva \
  --output projects/GSE12345/no_sva_results
```

</details>

---

## Known Limitations

### Sesame pthread errors on non-standard R builds

On certain R installations (particularly with OpenBLAS or non-standard threading configurations), Sesame's internal parallelization may trigger `pthread` errors or segmentation faults. IlluMeta automatically:
1. Detects these errors during Sesame normalization
2. Falls back to single-threaded mode (`SESAME_NTHREAD=1`)
3. Retries the normalization

If you see warnings like `"pthread_create: Resource temporarily unavailable"` or `"caught segfault"`, IlluMeta handles these gracefully. To force single-threaded mode from the start:

```bash
SESAME_NTHREAD=1 python3 illumeta.py analysis ...
```

### Small sample size limitations

IlluMeta implements a sample-size-adaptive Correction Robustness Framework (CRF):

| Tier | Total n | Per-group min | Key limitations |
|------|---------|---------------|-----------------|
| **Minimal** | < 12 | 3 | Exploratory only; severe power limitation; SVA disabled |
| **Small** | 12-23 | 6 | Independent validation required; SVA disabled |
| **Moderate** | 24-49 | 12 | Full pipeline available; moderate power |
| **Large** | >= 50 | 25 | All robustness checks fully powered |

For minimal/small tiers, interpret results cautiously and plan for replication in independent cohorts.

### EPIC v2 (EPICv2) array support

IlluMeta supports EPIC v2 arrays with automatic manifest detection. However:
- Some probe annotations may differ between EPIC v1 and v2
- Cross-reactive probe lists are v1-based; v2-specific lists are being incorporated
- DMR analysis uses EPIC v1 annotations; v2-specific boundaries may vary slightly

### Multi-batch studies

When analyzing studies with multiple batches:
- ComBat assumes balanced designs; highly imbalanced batches may introduce bias
- For >=3 batches with small per-batch n, consider stratified analysis (`--tier3-batch`)
- Cross-platform (450K + EPIC) studies should subset to overlapping probes first

## Citing IlluMeta

If you use IlluMeta in your research, please cite:

> Kang, K. (2026). IlluMeta: an automated dual-pipeline framework for robust DNA methylation analysis. *Bioinformatics* (submitted).

See `CITATION.cff` for machine-readable citation metadata.

## License
Apache-2.0 (see `LICENSE`).
