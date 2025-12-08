# IlluMeta
*IlluMeta: Illuminating Methylation Analytics*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.0000000.svg)](https://doi.org/10.5281/zenodo.0000000)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/kangk1204/illumeta)

A publication-ready, fully automated toolkit to go from GEO accession (or raw IDAT files) to interpretable methylation results in a single command. IlluMeta automates data retrieval, quality control (QC), normalization (Minfi/Sesame), batch effect correction (Combat/SVA), and generates an interactive HTML dashboard.

## Key Features
- **One-Command Workflow:** From `GSEID` to final results.
- **Double Validation:** Intersects results from **Minfi (Noob)** and **Sesame** pipelines to minimize false positives.
- **Advanced Batch Correction:** Automatically detects and corrects batch effects using SVA (Surrogate Variable Analysis) and ComBat.
- **Interactive Reporting:** Generates a full HTML dashboard with Volcano plots, Heatmaps, Manhattan plots, and QC metrics.
- **Cross-Platform:** Runs seamlessly on **Ubuntu/Linux**, **macOS (Intel & Apple Silicon)**, and **Windows 11 (via WSL2)**.

---

## 🚀 Installation & Setup

### 1. Prerequisites
Ensure you have the following installed:
- **Python** (3.8+)
- **R** (4.2+)
- **Git**
- **System Libraries** (for R packages):
  - **Ubuntu/WSL:** `sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev cmake pandoc`
  - **macOS:** `brew install openssl@3 libxml2 cmake pandoc`

### 2. Clone & Install
```bash
# 1. Clone the repository
git clone https://github.com/kangk1204/illumeta.git
cd illumeta

# 2. Set up Python environment
python3 -m venv .venv
source .venv/bin/activate

# 3. Install Python dependencies
pip install --upgrade pip
pip install -r requirements.txt

# 4. (Optional but recommended) First-time R setup
# This pre-installs all R packages and references (takes ~10-15 mins).
# If skipped, IlluMeta will attempt to install them during the first run.
ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R
```

---

## 🏃 Quick Start (Demo)

Try the pipeline with a small public dataset (**GSE121633**, ~4 samples) to verify your installation.

### Step 1: Download Data
```bash
# Activate venv if not already active
source .venv/bin/activate

# Download IDAT files for the demo dataset
python3 illumeta.py download GSE121633 -o projects/GSE121633
```

### Step 2: Configure Groups
The download creates a `configure.tsv` file in `projects/GSE121633/`. You need to tell IlluMeta which samples are "Control" and which are "Test".

1. Open `projects/GSE121633/configure.tsv`.
2. Locate the `primary_group` column (it will be empty).
3. Fill it in. For this demo, let's pretend:
   - First 2 samples -> `Control`
   - Last 2 samples -> `Case`
4. Save the file.

### Step 3: Run Analysis
```bash
python3 illumeta.py analysis \
    -i projects/GSE121633 \
    --group_con Control \
    --group_test Case
```

### Step 4: View Results
Open the generated dashboard in your browser:
`projects/GSE121633/Case_vs_Control_results_index.html`

---

## 🛠️ Usage Guide

### 1. Search for Datasets
Find human methylation datasets on GEO that have raw IDAT files.
```bash
python3 illusearch.py --keywords "breast cancer" --email your_email@example.com -o search_results.tsv
```

### 2. Analyze Your Own Data (Non-GEO)
If you have your own `.idat` files:
1. Create a folder: `my_project/idat/`
2. Put your `_Grn.idat` and `_Red.idat` files there.
3. Create a `configure.tsv` in `my_project/`. It **must** have:
   - `Basename`: Path prefix to the IDATs (e.g., `idat/Sample1_R01C01`).
   - `primary_group`: Your study groups (e.g., `Treated`, `Untreated`).
4. Run the analysis:
   ```bash
   python3 illumeta.py analysis -i my_project --group_con Untreated --group_test Treated
   ```

### 3. Advanced Options
```bash
# Tighten statistical thresholds
python3 illumeta.py analysis -i projects/GSE12345 --pval 0.01 --lfc 1.0 ...

# Disable SVA (if sample size is very small)
python3 illumeta.py analysis -i projects/GSE12345 --disable-sva ...

# Use a custom temporary directory (for large datasets)
python3 illumeta.py analysis -i projects/GSE12345 --tmp-dir /mnt/big_drive/tmp ...
```

---

## ⚠️ Troubleshooting

**Common Issues:**

- **`ERROR: Some required R packages failed to install`**:
  Usually due to missing system libraries.
  - **Ubuntu:** Run `sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev libicu-dev`
  - **macOS:** Run `brew install libxml2 openssl@3 curl`
  - Then force a retry: `ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R`

- **`pandoc: command not found`**:
  Required for report generation.
  - **Ubuntu:** `sudo apt-get install pandoc`
  - **macOS:** `brew install pandoc`

- **`C compiler cannot create executables` (macOS/Conda)**:
  If you have Conda installed, it might interfere with R's compilers.
  Try running this before setup:
  ```bash
  export R_ENVIRON_USER=/dev/null
  unset CC CXX CXX11 CXX14 CXX17 CXX20
  ```

---

## 📄 Citation
If you use IlluMeta in your research, please cite:

> Kang, K. (2025). IlluMeta: Automated Pipeline for Illumina DNA Methylation Analysis. GitHub. https://github.com/kangk1204/illumeta

Or use the `CITATION.cff` file provided in the repository.

## 📜 License
This project is licensed under the **Apache 2.0 License**.