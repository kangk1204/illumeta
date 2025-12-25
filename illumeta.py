#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import csv
import json
import re
import unicodedata
from datetime import datetime

__version__ = "1.0.0"

# Configuration for R script paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def safe_path_component(name: str, fallback: str = "results") -> str:
    """Return an ASCII-only path component safe across filesystems/locales."""
    if not name:
        return fallback
    normalized = unicodedata.normalize("NFKD", name)
    ascii_name = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_name = re.sub(r"\s+", "_", ascii_name)
    ascii_name = re.sub(r"[^A-Za-z0-9._-]+", "_", ascii_name)
    ascii_name = re.sub(r"_+", "_", ascii_name).strip("._-")
    return ascii_name or fallback

def detect_r_major_minor(env=None):
    override = (env or {}).get("ILLUMETA_R_LIB_VERSION") or os.environ.get("ILLUMETA_R_LIB_VERSION")
    if override:
        return override
    cmd = [
        "Rscript",
        "-e",
        "cat(paste0(R.version$major, '.', sub('\\\\..*', '', R.version$minor)))",
    ]
    try:
        res = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except FileNotFoundError:
        res = None
    if res and res.returncode == 0:
        version = res.stdout.strip()
        if version:
            return version
    try:
        res = subprocess.run(
            ["R", "--version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except FileNotFoundError:
        return None
    if res.returncode != 0:
        return None
    match = re.search(r"R version (\\d+)\\.(\\d+)", res.stdout)
    if not match:
        return None
    return f"{match.group(1)}.{match.group(2)}"

def default_r_lib_base(base_dir: str, env=None) -> str:
    allow_external = (env or {}).get("ILLUMETA_ALLOW_EXTERNAL_LIB") == "1" or (
        os.environ.get("ILLUMETA_ALLOW_EXTERNAL_LIB") == "1"
    )
    if sys.platform == "darwin" and base_dir.startswith("/Volumes/") and not allow_external:
        return os.path.join(os.path.expanduser("~"), ".illumeta", "r-lib")
    return os.path.join(base_dir, ".r-lib")

def default_r_lib(base_dir: str, env=None) -> str:
    lib_root = default_r_lib_base(base_dir, env)
    r_version = detect_r_major_minor(env)
    if r_version:
        return os.path.join(lib_root, f"R-{r_version}")
    return lib_root

def read_recorded_r_lib(base_dir: str):
    record_path = os.path.join(base_dir, ".illumeta_r_lib_path")
    if not os.path.isfile(record_path):
        return None
    try:
        with open(record_path, "r", encoding="utf-8") as handle:
            value = handle.read().strip()
    except OSError:
        return None
    return value or None
R_SCRIPTS_DIR = os.path.join(BASE_DIR, "r_scripts")
DOWNLOAD_SCRIPT = os.path.join(R_SCRIPTS_DIR, "download.R")
ANALYZE_SCRIPT = os.path.join(R_SCRIPTS_DIR, "analyze.R")
SETUP_MARKER = os.path.join(BASE_DIR, ".r_setup_done")
SETUP_SCRIPT = os.path.join(R_SCRIPTS_DIR, "setup_env.R")
DEFAULT_CONDA_PREFIX = os.environ.get("CONDA_PREFIX")
CORE_R_PACKAGES = [
    "xml2",
    "XML",
    "optparse",
    "GEOquery",
    "minfi",
    "sesame",
    "limma",
    "dmrff",
    "sesameData",
    "sva",
    "variancePartition",
    "pvca",
    "illuminaio",
    "RefFreeEWAS",
    "Biobase",
    "reformulas",
    "ggplot2",
    "plotly",
    "DT",
    "data.table",
    "htmlwidgets",
    "dplyr",
    "stringr",
    "ggrepel",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
]
EPICV2_R_PACKAGES = [
    "IlluminaHumanMethylationEPICv2manifest",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
]
OPTIONAL_R_PACKAGES = [
    "FlowSorted.Blood.EPIC",
    "FlowSorted.Blood.450k",
    "wateRmelon",
    "methylclock",
    "methylclockData",
    "planet",
]
def add_conda_paths(env: dict) -> dict:
    """Ensure LD_LIBRARY_PATH/PKG_CONFIG_PATH/PATH include conda libs so xml2/xml load correctly."""
    use_conda_libs = (env.get("ILLUMETA_USE_CONDA_LIBS") == "1") or (
        os.environ.get("ILLUMETA_USE_CONDA_LIBS") == "1"
    )
    if not use_conda_libs:
        return env
    prefix = env.get("CONDA_PREFIX") or os.environ.get("CONDA_PREFIX")
    if not prefix or not os.path.isdir(prefix):
        return env
    conda_lib = os.path.join(prefix, "lib")
    conda_pkgconfig = os.path.join(conda_lib, "pkgconfig")
    conda_bin = os.path.join(prefix, "bin")
    def prepend(path_var, new_path):
        cur = env.get(path_var, "")
        parts = cur.split(os.pathsep) if cur else []
        if new_path and os.path.exists(new_path) and new_path not in parts:
            env[path_var] = os.pathsep.join([new_path] + parts if parts else [new_path])
    prepend("LD_LIBRARY_PATH", conda_lib)
    prepend("PKG_CONFIG_PATH", conda_pkgconfig)
    prepend("PATH", conda_bin)
    return env

def log(msg: str):
    """Prints a timestamped info message."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

def log_err(msg: str):
    """Prints a timestamped error message to stderr."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", file=sys.stderr)

def check_r_installation():
    """Checks if R is available in the path."""
    try:
        subprocess.run(["R", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        log_err("Error: R is not installed or not in the PATH.")
        log_err("  - Ubuntu/Debian: sudo apt-get update && sudo apt-get install r-base")
        log_err("  - macOS (Homebrew): brew install --cask r")
        log_err("  - Conda: conda install -c conda-forge r-base")
        log_err("After installation, ensure 'R' and 'Rscript' are in your PATH and re-run illumeta.py.")
        sys.exit(1)

def check_pandoc_installation():
    """Checks if pandoc is available in the path."""
    try:
        subprocess.run(["pandoc", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        log_err("Error: 'pandoc' is not installed or not in the PATH.")
        log_err("Pandoc is required for generating interactive HTML reports.")
        log_err("  - macOS: brew install pandoc")
        log_err("  - Linux: sudo apt-get install pandoc (Debian/Ubuntu) or yum install pandoc (CentOS/RHEL)")
        log_err("  - Windows: choco install pandoc")
        sys.exit(1)

def check_r_package(pkg, env):
    expr = (
        'lib <- Sys.getenv("R_LIBS_USER"); '
        'if (nzchar(lib)) .libPaths(c(lib, .Library, .Library.site)); '
        f'ok <- requireNamespace("{pkg}", quietly=TRUE); '
        'cat(if (isTRUE(ok)) "OK" else "MISSING")'
    )
    try:
        res = subprocess.run(
            ["Rscript", "-e", expr],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except FileNotFoundError:
        return False, "Rscript not found"
    if res.returncode == 0 and "OK" in res.stdout:
        return True, ""
    err = (res.stderr + res.stdout).strip()
    return False, err

def check_r_packages(pkgs, env=None):
    if not pkgs:
        return [], {}
    if env is None:
        env = ensure_r_lib_env(os.environ.copy())
        env = add_conda_paths(env)
    missing = []
    errors = {}
    for pkg in pkgs:
        ok, err = check_r_package(pkg, env)
        if not ok:
            missing.append(pkg)
            if err:
                errors[pkg] = err
    return missing, errors

def missing_r_packages(pkgs, env=None):
    """Returns a list of missing or unloadable R packages."""
    missing, _ = check_r_packages(pkgs, env=env)
    return missing

def format_r_error(err: str) -> str:
    if not err:
        return ""
    lower = err.lower()
    if "segfault" in lower:
        return "R crashed while loading this package (segfault)"
    line = err.splitlines()[0].strip()
    return line[:200]

def ensure_r_lib_env(env):
    """
    Ensure R uses a repo-local, writable library by default for reproducibility.
    To respect an externally-set R_LIBS_USER instead, set ILLUMETA_RESPECT_R_LIBS_USER=1.
    """
    respect_existing = (env.get("ILLUMETA_RESPECT_R_LIBS_USER") == "1") or (
        os.environ.get("ILLUMETA_RESPECT_R_LIBS_USER") == "1"
    )
    existing = env.get("R_LIBS_USER")
    if respect_existing and existing:
        return env

    recorded_lib = read_recorded_r_lib(BASE_DIR)
    default_lib_root = default_r_lib_base(BASE_DIR, env)
    default_lib = recorded_lib or default_r_lib(BASE_DIR, env)
    if recorded_lib:
        default_lib_root = os.path.dirname(recorded_lib.rstrip(os.sep))
    if existing and existing != default_lib and not respect_existing:
        log(f"[*] Overriding R_LIBS_USER={existing} -> {default_lib} (set ILLUMETA_RESPECT_R_LIBS_USER=1 to keep)")

    if BASE_DIR.startswith("/Volumes/") and default_lib_root != os.path.join(BASE_DIR, ".r-lib") and not respect_existing:
        log(f"[*] Project is on external volume; using local R library at {default_lib} (set ILLUMETA_ALLOW_EXTERNAL_LIB=1 to use .r-lib)")

    env["R_LIBS_USER"] = default_lib
    if not os.path.exists(default_lib):
        os.makedirs(default_lib, exist_ok=True)
        log(f"[*] Created default R library at {default_lib}")
    return env

def ensure_r_dependencies():
    """Runs the setup_env.R script to install required R packages and cache data."""
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)
    force_setup = os.environ.get("ILLUMETA_FORCE_SETUP") == "1"
    require_epicv2 = os.environ.get("ILLUMETA_REQUIRE_EPICV2") == "1"
    marker_ok = False
    marker_epicv2_required = None
    if os.path.exists(SETUP_MARKER):
        try:
            with open(SETUP_MARKER) as f:
                content = f.read()
            for line in content.splitlines():
                if line.startswith("epicv2_required="):
                    marker_epicv2_required = line.split("=", 1)[1].strip() == "1"
            if marker_epicv2_required is None:
                marker_epicv2_required = False
            marker_ok = marker_epicv2_required == require_epicv2
        except Exception:
            marker_ok = False
    core_pkgs = CORE_R_PACKAGES + (EPICV2_R_PACKAGES if require_epicv2 else [])
    optional_pkgs = OPTIONAL_R_PACKAGES + ([] if require_epicv2 else EPICV2_R_PACKAGES)
    missing_core = missing_r_packages(core_pkgs, env=env) if marker_ok else []
    missing_optional = missing_r_packages(optional_pkgs, env=env) if marker_ok else []
    if marker_ok and not force_setup and not missing_core:
        if missing_optional:
            log(f"[*] Optional R packages missing (analysis will skip related features): {', '.join(missing_optional)}")
            if any(pkg in EPICV2_R_PACKAGES for pkg in missing_optional):
                log("[*] EPIC v2 support is not installed. Use R 4.4+ and set ILLUMETA_REQUIRE_EPICV2=1 to require it.")
        log("[*] R dependencies already set up (skipping). Set ILLUMETA_FORCE_SETUP=1 to force reinstall.")
        return
    if marker_ok and missing_core and not force_setup:
        log(f"[*] R setup marker found but missing packages: {', '.join(missing_core)}. Re-running setup_env.R ...")

    log("[*] Ensuring R dependencies (this may take a few minutes on first run)...")
    try:
        subprocess.run(["Rscript", SETUP_SCRIPT], check=True, env=env)
        with open(SETUP_MARKER, "w") as f:
            f.write(f"setup completed at {datetime.now().isoformat()}\n")
            f.write(f"epicv2_required={1 if require_epicv2 else 0}\n")
        missing_optional_after = missing_r_packages(optional_pkgs, env=env)
        if missing_optional_after:
            log(f"[*] Optional R packages missing (features will be skipped): {', '.join(missing_optional_after)}")
            if any(pkg in EPICV2_R_PACKAGES for pkg in missing_optional_after):
                log("[*] EPIC v2 support is not installed. Use R 4.4+ and set ILLUMETA_REQUIRE_EPICV2=1 to require it.")
    except subprocess.CalledProcessError as e:
        log_err(f"[!] Error while installing R dependencies: {e}")
        log_err("    Hint: If you see 'library path not writable', set a user library first:")
        log_err("      export R_LIBS_USER=\"$HOME/R/library\" && mkdir -p \"$R_LIBS_USER\"")
        log_err("    Then rerun: Rscript r_scripts/setup_env.R")
        if os.path.exists(SETUP_MARKER):
            try:
                os.remove(SETUP_MARKER)
            except OSError:
                pass
        sys.exit(1)

def run_doctor(args):
    """Checks system and R dependencies without installing anything."""
    check_r_installation()
    if not args.skip_pandoc:
        check_pandoc_installation()

    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)

    log("[*] Checking required R packages...")
    require_epicv2 = os.environ.get("ILLUMETA_REQUIRE_EPICV2") == "1"
    core_pkgs = CORE_R_PACKAGES + (EPICV2_R_PACKAGES if require_epicv2 else [])
    optional_pkgs = OPTIONAL_R_PACKAGES + ([] if require_epicv2 else EPICV2_R_PACKAGES)
    missing_core, core_errors = check_r_packages(core_pkgs, env=env)
    missing_optional, optional_errors = check_r_packages(optional_pkgs, env=env)

    if missing_core:
        log_err("[!] Missing required R packages:")
        saw_segfault = False
        for pkg in missing_core:
            err = format_r_error(core_errors.get(pkg, ""))
            if "segfault" in err.lower():
                saw_segfault = True
            if err:
                log_err(f"  - {pkg}: {err}")
            else:
                log_err(f"  - {pkg}")
        if saw_segfault:
            log_err("    Hint: A package crashed R while loading. Rerun setup with `ILLUMETA_CLEAN_MISMATCHED_RLIB=1`.")
        log_err("    Fix: run `ILLUMETA_FORCE_SETUP=1 Rscript r_scripts/setup_env.R` (or run any IlluMeta command once).")
        sys.exit(1)

    log("[*] Core R packages: OK")
    if missing_optional:
        log("[*] Optional R packages missing (features will be skipped): " + ", ".join(missing_optional))
        for pkg in missing_optional:
            err = format_r_error(optional_errors.get(pkg, ""))
            if err:
                log(f"[*]   - {pkg}: {err}")
        if any(pkg in EPICV2_R_PACKAGES for pkg in missing_optional):
            log("[*] EPIC v2 support is not installed. Use R 4.4+ and set ILLUMETA_REQUIRE_EPICV2=1 to require it.")
    else:
        log("[*] Optional R packages: OK")

def run_download(args):
    """Executes the download step."""
    out_dir = args.out_dir if args.out_dir else os.path.abspath(args.gse_id)
    log(f"[*] Starting download pipeline for {args.gse_id}...")
    log(f"[*] Output directory: {out_dir}")
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    cmd = ["Rscript", DOWNLOAD_SCRIPT, "--gse", args.gse_id, "--out", out_dir]
    
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)

    try:
        subprocess.run(cmd, check=True, env=env)
        log(f"[*] Download complete. Please edit the 'primary_group' column in: {os.path.join(out_dir, 'configure.tsv')}")
    except subprocess.CalledProcessError as e:
        log_err(f"[!] Error during download step: {e}")
        sys.exit(1)

def run_analysis(args):
    """Executes the analysis step."""
    # Resolve Config Path
    config_path = None
    if args.config:
        config_path = args.config
    elif args.input_dir:
        config_path = os.path.join(args.input_dir, "configure.tsv")
    
    if not config_path:
        log_err("Error: You must provide either --input-dir or --config.")
        sys.exit(1)

    log(f"[*] Starting analysis pipeline using configuration: {config_path}...")
    
    if not os.path.exists(config_path):
        log_err(f"Error: Configuration file {config_path} not found.")
        sys.exit(1)
        
    # Determine output directory
    if args.output:
        output_dir = args.output
    else:
        # Default: test_vs_con_results (sanitize for filesystem/locale safety)
        dir_base = args.input_dir if args.input_dir else os.path.dirname(os.path.abspath(config_path))
        folder_name = f"{args.group_test}_vs_{args.group_con}_results"
        safe_folder = safe_path_component(folder_name, fallback="analysis_results")
        if safe_folder != folder_name:
            log(f"[*] Output folder normalized for filesystem safety: {folder_name} -> {safe_folder}")
        output_dir = os.path.join(dir_base, safe_folder)
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    cmd = [
        "Rscript", ANALYZE_SCRIPT, 
        "--config", config_path, 
        "--out", output_dir, 
        "--group_con", args.group_con,
        "--group_test", args.group_test,
        "--max_plots", str(args.max_plots),
        "--pval", str(args.pval),
        "--lfc", str(args.lfc),
        "--min_total_size", str(args.min_total_size),
        "--qc_intensity_threshold", str(args.qc_intensity_threshold)
    ]
    if args.force_idat:
        cmd.append("--force_idat")
    if args.disable_auto_covariates:
        cmd.append("--disable_auto_covariates")
    if args.disable_sva:
        cmd.append("--disable_sva")
    if args.include_covariates:
        cmd.extend(["--include_covariates", args.include_covariates])
    if args.include_clock_covariates:
        cmd.append("--include_clock_covariates")
    if args.tissue:
        cmd.extend(["--tissue", args.tissue])
    if args.positive_controls:
        cmd.extend(["--positive_controls", args.positive_controls])
    if args.permutations and args.permutations > 0:
        cmd.extend(["--permutations", str(args.permutations)])
    if args.vp_top:
        cmd.extend(["--vp_top", str(args.vp_top)])
    if args.id_column:
        cmd.extend(["--id_column", args.id_column])
    
    # Handle custom temp directory
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)
    if args.tmp_dir:
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir)
        env["TMPDIR"] = os.path.abspath(args.tmp_dir)
        log(f"[*] Using custom temporary directory: {env['TMPDIR']}")

    try:
        subprocess.run(cmd, check=True, env=env)
        log(f"[*] Analysis complete. Results are in: {output_dir}")
        
        # Generate Dashboard
        generate_dashboard(output_dir, args.group_test, args.group_con)
        
    except subprocess.CalledProcessError as e:
        log_err(f"[!] Error during analysis step: {e}")
        sys.exit(1)

def safe_int(val):
    try:
        return int(val)
    except (ValueError, TypeError):
        return 0

def generate_dashboard(output_dir, group_test, group_con):
    """Generates a beautiful HTML dashboard to navigate results."""
    results_folder_name = os.path.basename(os.path.normpath(output_dir))
    parent_dir = os.path.dirname(os.path.normpath(output_dir))
    dashboard_filename = f"{results_folder_name}_index.html"
    dashboard_path = os.path.join(parent_dir, dashboard_filename)
    
    # Reordered: Intersection First
    pipelines = ["Intersection", "Minfi", "Sesame"]
    
    # Load Summary Statistics
    summary_path = os.path.join(output_dir, "summary.json")
    stats = {}
    try:
        with open(summary_path, "r") as f:
            stats = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        log_err("[!] Warning: summary.json not found or invalid. Dashboard will lack stats.")
        # Default empty stats
        keys = ["n_con", "n_test", "minfi_up", "minfi_down", "sesame_up", "sesame_down", "intersect_up", "intersect_down"]
        stats = {k: 0 for k in keys}

    intersection_sections = [
        ("Consensus Highlights", [
            ("_Consensus_DMPs.html", "Consensus DMPs Table", "CpGs significant in both pipelines with consistent direction.", "TABLE"),
            ("_Consensus_DMPs.csv", "Consensus DMPs CSV", "Consensus DMP list (CSV).", "CSV"),
        ]),
        ("Agreement Checks", [
            ("_LogFC_Concordance.html", "Pipeline Concordance", "Minfi vs Sesame logFC concordance.", "PLOT"),
            ("_Significant_Overlap.html", "Pipeline Overlap", "Significant DMP counts and overlap.", "PLOT"),
        ]),
        ("Visual Summary", [
            ("_Volcano.html", "Volcano Plot", "Visualizes significant changes (LogFC vs P-value).", "PLOT"),
            ("_Manhattan.html", "Manhattan Plot", "Genomic distribution of methylation changes.", "PLOT"),
            ("_Top100_Heatmap.html", "Heatmap (Top 100)", "Clustering of the top 100 most significant CpGs.", "PLOT"),
        ]),
    ]

    pipeline_sections_base = [
        ("Primary Results", [
            ("_Volcano.html", "Volcano Plot", "Visualizes significant changes (LogFC vs P-value).", "PLOT"),
            ("_Manhattan.html", "Manhattan Plot", "Genomic distribution of methylation changes.", "PLOT"),
            ("_Top100_Heatmap.html", "Heatmap (Top 100)", "Clustering of the top 100 most significant CpGs.", "PLOT"),
            ("_Top_DMPs.html", "Top DMPs Table", "Searchable table of differentially methylated probes.", "TABLE"),
        ]),
        ("Region-Level Results", [
            ("_DMR_Volcano.html", "DMR Volcano Plot", "Visualizes significant DMRs (Est. Diff vs P-value).", "PLOT"),
            ("_DMR_Manhattan.html", "DMR Manhattan Plot", "Genomic distribution of DMRs.", "PLOT"),
            ("_Top_DMRs_Heatmap.html", "DMR Heatmap (Top 50)", "Average methylation levels of top 50 DMRs.", "PLOT"),
            ("_DMRs_Table.html", "DMRs Table", "Searchable table of differentially methylated regions.", "TABLE"),
        ]),
        ("Quality & Diagnostics", [
            ("_PCA_Before.html", "PCA (Raw)", "Principal Component Analysis before correction.", "PLOT"),
            ("_PCA_After_Correction.html", "PCA (Corrected)", "PCA after removing detected batch effects.", "PLOT"),
            ("_Sample_Clustering_Distance.html", "Sample Clustering", "Euclidean distance matrix of samples.", "PLOT"),
            ("_QQPlot.html", "Q-Q Plot", "Check for genomic inflation and systematic bias.", "PLOT"),
            ("_PVCA.html", "PVCA (Raw)", "Variance explained by group/batch/covariates.", "PLOT"),
            ("_AfterCorrection_PVCA.html", "PVCA (Corrected)", "Variance explained after correction.", "PLOT"),
        ]),
        ("Batch & Covariates", [
            ("_Batch_Method_Comparison.html", "Correction Method", "Batch residual vs group signal by method.", "PLOT"),
            ("_BatchMethodComparison.csv", "Correction Method (CSV)", "Batch correction scoring table (CSV).", "CSV"),
            ("_Batch_Evaluation_Before.html", "Batch Eval (Before)", "Covariate association with PCs (Raw).", "PLOT"),
            ("_Batch_Evaluation_After.html", "Batch Eval (After)", "Covariate association with PCs (Corrected).", "PLOT"),
            ("_AutoCovariates.csv", "Auto Covariates (CSV)", "Auto-selected covariate candidates with PC association P-values.", "CSV"),
            ("_DroppedCovariates.csv", "Dropped Covariates (CSV)", "Covariates removed due to confounding/instability.", "CSV"),
            ("_configure_with_clocks.tsv", "Clock-Merged Config (TSV)", "Config with clock covariates merged (if enabled).", "TSV"),
        ]),
        ("Clocks & Age", [
            ("_Epigenetic_Age_methylclock.html", "Epigenetic Age (methylclock)", "Clock vs chronological age (if provided).", "PLOT"),
            ("_Epigenetic_Age_methylclock.csv", "Epigenetic Age Table", "DNAm clock predictions (methylclock).", "CSV"),
            ("_Placental_Age_planet_scatter.html", "Placental Age Scatter", "RPC vs CPC gestational age (Placenta).", "PLOT"),
            ("_Placental_Age_planet.csv", "Placental Age (planet)", "RPC/CPC gestational age predictions (Placenta).", "CSV"),
        ]),
    ]

    pipeline_sections = {
        "Intersection": intersection_sections,
        "Minfi": pipeline_sections_base,
        "Sesame": pipeline_sections_base,
    }
    
    def load_metrics(pipe_name):
        path = os.path.join(output_dir, f"{pipe_name}_Metrics.csv")
        metrics = {}
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        metrics[row["metric"]] = row["value"]
            except Exception:
                pass
        return metrics

    def load_analysis_params():
        path = os.path.join(output_dir, "analysis_parameters.json")
        try:
            with open(path, "r") as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def load_qc_summary():
        path = os.path.join(output_dir, "QC_Summary.csv")
        metrics = {}
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        metrics[row["metric"]] = row["value"]
            except Exception:
                pass
        return metrics

    def load_drop_reasons(pipe_name):
        path = os.path.join(output_dir, f"{pipe_name}_DroppedCovariates.csv")
        reasons = {}
        if os.path.exists(path):
            try:
                with open(path, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        reason = (row.get("Reason") or "").strip()
                        if not reason or reason == "NA":
                            continue
                        reasons[reason] = reasons.get(reason, 0) + 1
            except Exception:
                pass
        return reasons

    def pill_link(filename, label):
        fpath = os.path.join(output_dir, filename)
        if os.path.exists(fpath):
            rel_path = f"{results_folder_name}/{filename}"
            return f'<a class="pill" href="{rel_path}" target="_blank">{label}</a>'
        return f'<span class="pill muted">{label}</span>'

    # Calculated values for Summary
    n_con = safe_int(stats.get('n_con'))
    n_test = safe_int(stats.get('n_test'))
    total_samples = n_con + n_test
    
    minfi_up = safe_int(stats.get('minfi_up'))
    minfi_down = safe_int(stats.get('minfi_down'))
    minfi_total = minfi_up + minfi_down
    
    sesame_up = safe_int(stats.get('sesame_up'))
    sesame_down = safe_int(stats.get('sesame_down'))
    sesame_total = sesame_up + sesame_down
    
    intersect_up = safe_int(stats.get('intersect_up'))
    intersect_down = safe_int(stats.get('intersect_down'))
    intersect_total = intersect_up + intersect_down
    analysis_params = load_analysis_params()
    qc_summary = load_qc_summary()

    # CSS Style Block (Using format to avoid curly brace hell)
    style_block = """
        @import url('https://fonts.googleapis.com/css2?family=Manrope:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
        :root {
            --ink: #1d2628;
            --primary: #2f6b64;
            --accent: #2f7a7b;
            --accent-2: #f2b45d;
            --rose: #e07a5f;
            --bg: #f6f1e8;
            --card: #ffffff;
            --line: #e4ddd2;
            --muted: #6f7a76;
            --shadow: 0 12px 30px rgba(27, 38, 40, 0.12);
        }
        html { scroll-behavior: smooth; }
        body {
            font-family: "Manrope", "IBM Plex Sans", "Helvetica Neue", sans-serif;
            margin: 0;
            color: var(--ink);
            background-color: var(--bg);
            background-image:
                radial-gradient(circle at 10% 10%, rgba(242, 180, 93, 0.12), transparent 45%),
                radial-gradient(circle at 90% 20%, rgba(47, 122, 123, 0.12), transparent 40%),
                radial-gradient(circle at 50% 85%, rgba(224, 122, 95, 0.08), transparent 45%);
        }
        h1, h2, h3 { font-family: "Space Grotesk", "Manrope", sans-serif; margin: 0; }
        header {
            background: linear-gradient(120deg, #1f3f3c, #2f6b64 55%, #3d8b7d);
            color: #f7f5f0;
            padding: 3.2rem 1.5rem 3.5rem;
        }
        .hero {
            max-width: 1200px;
            margin: 0 auto;
            display: grid;
            grid-template-columns: minmax(260px, 1.2fr) minmax(240px, 0.8fr);
            gap: 28px;
            align-items: center;
        }
        .hero-eyebrow {
            text-transform: uppercase;
            font-size: 0.75rem;
            letter-spacing: 0.2em;
            color: rgba(247, 245, 240, 0.65);
        }
        .hero-title { font-size: clamp(2rem, 2.8vw, 2.8rem); margin-top: 0.5rem; }
        .hero-sub { margin-top: 0.6rem; font-size: 1rem; opacity: 0.85; }
        .hero-tags { display: flex; flex-wrap: wrap; gap: 0.5rem; margin-top: 1.2rem; }
        .tag {
            background: rgba(255, 255, 255, 0.12);
            padding: 0.35rem 0.75rem;
            border-radius: 999px;
            font-size: 0.85rem;
        }
        .hero-grid { display: grid; gap: 12px; }
        .hero-card {
            background: rgba(255, 255, 255, 0.12);
            border: 1px solid rgba(255, 255, 255, 0.2);
            border-radius: 14px;
            padding: 1rem 1.2rem;
            backdrop-filter: blur(8px);
        }
        .hero-card-title { font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.12em; opacity: 0.7; }
        .hero-card-value { font-size: 1.8rem; font-weight: 700; margin-top: 0.3rem; }
        .hero-card-sub { font-size: 0.9rem; opacity: 0.8; margin-top: 0.4rem; }

        .container { max-width: 1200px; margin: -2.2rem auto 3rem; padding: 0 20px 2rem; }
        .jump-bar {
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            padding: 0.75rem 1rem;
            border-radius: 14px;
            background: var(--card);
            box-shadow: var(--shadow);
            position: sticky;
            top: 0;
            z-index: 100;
        }
        .jump-bar a {
            text-decoration: none;
            color: var(--ink);
            font-weight: 600;
            font-size: 0.9rem;
            padding: 0.35rem 0.75rem;
            border-radius: 999px;
            background: #f2eee6;
        }
        .jump-bar a:hover { background: #e7dfd2; }

        .section-title {
            font-size: 1.3rem;
            font-weight: 700;
            color: var(--primary);
            margin: 2.2rem 0 0.8rem;
        }
        .tab-content .section-title {
            font-size: 1.05rem;
            color: var(--ink);
            margin-top: 1.6rem;
        }
        .section-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 18px; margin-bottom: 20px; }
        .callout {
            background: #fdf9f2;
            border: 1px solid var(--line);
            border-radius: 14px;
            padding: 1rem 1.2rem;
            color: var(--muted);
            box-shadow: 0 8px 22px rgba(27, 38, 40, 0.06);
        }
        .callout strong { color: var(--ink); }

        .step-card {
            background: var(--card);
            border: 1px solid var(--line);
            border-radius: 16px;
            padding: 1.2rem 1.4rem;
            box-shadow: var(--shadow);
        }
        .step-card h3 { margin-bottom: 0.4rem; }
        .step-label {
            font-size: 0.75rem;
            letter-spacing: 0.18em;
            text-transform: uppercase;
            color: var(--accent);
            margin-bottom: 0.4rem;
        }
        .pill-row { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 0.8rem; }
        .pill {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 0.35rem 0.7rem;
            border-radius: 999px;
            background: #eef4f2;
            color: var(--primary);
            text-decoration: none;
            font-size: 0.85rem;
            font-weight: 600;
        }
        .pill.muted { background: #f1f0ed; color: #9aa4a1; }

        .nav-tabs { display: flex; gap: 10px; flex-wrap: wrap; margin-top: 1.5rem; }
        .nav-tab {
            padding: 0.6rem 1.4rem;
            cursor: pointer;
            font-weight: 700;
            color: var(--primary);
            border-radius: 999px;
            background: #eef4f2;
            transition: all 0.2s ease;
        }
        .nav-tab.active { background: var(--accent); color: #fff; }
        .nav-tab:hover { transform: translateY(-2px); }

        .tab-shell { margin-top: 1.2rem; }
        .tab-content { display: none; animation: riseIn 0.5s ease; }
        .tab-content.active { display: block; }

        .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 18px; }
        .card {
            background: var(--card);
            border-radius: 16px;
            border: 1px solid var(--line);
            box-shadow: 0 10px 22px rgba(27, 38, 40, 0.08);
            overflow: hidden;
            display: flex;
            flex-direction: column;
            animation: riseIn 0.55s ease both;
            animation-delay: calc(var(--i, 1) * 40ms);
        }
        .card-header {
            padding: 0.9rem 1rem;
            font-weight: 700;
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: #f7f2ea;
        }
        .card-badge {
            font-size: 0.7rem;
            padding: 2px 8px;
            border-radius: 999px;
            font-weight: 700;
            letter-spacing: 0.08em;
        }
        .badge-plot { background: rgba(47, 122, 123, 0.15); color: var(--accent); }
        .badge-table { background: rgba(242, 180, 93, 0.2); color: #a86814; }
        .badge-csv, .badge-tsv { background: rgba(224, 122, 95, 0.2); color: #b05643; }
        .badge-doc { background: rgba(31, 63, 60, 0.15); color: #1f3f3c; }
        .card-body { padding: 1.2rem 1.2rem 1.4rem; flex-grow: 1; }
        .card-desc { font-size: 0.9rem; color: var(--muted); line-height: 1.5; margin-bottom: 1.2rem; }
        .btn {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 100%;
            padding: 0.65rem 0;
            background: var(--accent);
            color: #fff;
            text-decoration: none;
            border-radius: 10px;
            font-weight: 700;
            transition: all 0.2s ease;
        }
        .btn:hover { background: #256a6b; transform: translateY(-1px); }

        .empty-msg { text-align: center; padding: 2rem; color: var(--muted); font-style: italic; }
        .metrics-card { background: var(--card); border: 1px solid var(--line); border-radius: 14px; padding: 1rem 1.2rem; box-shadow: 0 8px 18px rgba(27, 38, 40, 0.06); }
        .metrics-title { font-weight: 700; margin-bottom: 0.6rem; color: var(--primary); }
        .metrics-row { display: flex; justify-content: space-between; gap: 12px; font-size: 0.95rem; padding: 6px 0; border-bottom: 1px dashed #eee4d8; }
        .metrics-row:last-child { border-bottom: none; }

        @keyframes riseIn { from { opacity: 0; transform: translateY(12px); } to { opacity: 1; transform: translateY(0); } }
        @media (max-width: 900px) {
            header { padding: 2.6rem 1.4rem 3rem; }
            .hero { grid-template-columns: 1fr; }
            .container { margin-top: -1.6rem; }
        }
    """

    # Construct HTML using lists and joins to be cleaner
    html_parts = []
    
    # Header
    html_parts.append(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IlluMeta Analysis: {group_test} vs {group_con}</title>
    <style>{style_block}</style>
</head>
<body>

<header>
    <div class="hero">
        <div>
            <div class="hero-eyebrow">IlluMeta Analysis</div>
            <h1 class="hero-title">Analysis Dashboard</h1>
            <p class="hero-sub">Comparison: <strong>{group_test}</strong> (Test) vs <strong>{group_con}</strong> (Control)</p>
            <div class="hero-tags">
                <span class="tag">Samples: {total_samples}</span>
                <span class="tag">Consensus DMPs: {intersect_total}</span>
                <span class="tag">Minfi: {minfi_total} / Sesame: {sesame_total}</span>
            </div>
        </div>
        <div class="hero-grid">
            <div class="hero-card">
                <div class="hero-card-title">Samples</div>
                <div class="hero-card-value">{total_samples}</div>
                <div class="hero-card-sub">Control {n_con} / Test {n_test}</div>
            </div>
            <div class="hero-card">
                <div class="hero-card-title">Consensus DMPs</div>
                <div class="hero-card-value">{intersect_total}</div>
                <div class="hero-card-sub">Up {intersect_up} / Down {intersect_down}</div>
            </div>
            <div class="hero-card">
                <div class="hero-card-title">Pipeline DMPs</div>
                <div class="hero-card-value">{minfi_total} | {sesame_total}</div>
                <div class="hero-card-sub">Minfi ▲ {minfi_up} ▼ {minfi_down} · Sesame ▲ {sesame_up} ▼ {sesame_down}</div>
            </div>
        </div>
    </div>
</header>

<div class="container">
    <div class="jump-bar">
        <a href="#start">Start Here</a>
        <a href="#controls">Run Controls & QC</a>
        <a href="#docs">Run Documentation</a>
        <a href="#pipelines">Results</a>
    </div>
""")

    guide_steps = [
        ("Step 1", "Check sample quality", "Verify QC pass/fail and signal quality before interpreting results.",
         [("QC_Summary.csv", "QC Summary"), ("Sample_QC_DetectionP_FailFraction.html", "Detection P"), ("Sample_QC_Intensity_Medians.html", "Intensity")]),
        ("Step 2", "Review consensus signals", "Use the intersection set for the most conservative findings.",
         [("Intersection_Consensus_DMPs.html", "Consensus DMPs"), ("Intersection_LogFC_Concordance.html", "Concordance"), ("Intersection_Significant_Overlap.html", "Overlap")]),
        ("Step 3", "Dive into pipelines", "Explore Minfi and Sesame for method-specific depth.",
         [("Minfi_Volcano.html", "Minfi Volcano"), ("Sesame_Volcano.html", "Sesame Volcano"), ("Minfi_DMRs_Table.html", "Minfi DMRs")]),
    ]
    html_parts.append('<div id="start" class="section-title">Beginner Path</div>')
    html_parts.append('<div class="section-grid">')
    for idx, (label, title, desc, links) in enumerate(guide_steps, start=1):
        link_html = " ".join([pill_link(fname, lbl) for fname, lbl in links])
        html_parts.append(f'''<div class="step-card" style="--i:{idx};">
            <div class="step-label">{label}</div>
            <h3>{title}</h3>
            <p class="card-desc">{desc}</p>
            <div class="pill-row">{link_html}</div>
        </div>''')
    html_parts.append('</div>')
    html_parts.append('<div class="callout"><strong>Beginner tip:</strong> If results look inconsistent, check QC and batch diagnostics before focusing on DMPs/DMRs.</div>')

    if analysis_params or qc_summary:
        html_parts.append('<div id="controls" class="section-title">Run Controls & QC</div>')
        html_parts.append('<div class="section-grid">')
        if analysis_params:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">Analysis Parameters</div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>FDR / |logFC|</span><span>{analysis_params.get("pval_threshold", "N/A")} / {analysis_params.get("lfc_threshold", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Auto covariates</span><span>alpha={analysis_params.get("auto_covariate_alpha", "N/A")}, max_pcs={analysis_params.get("auto_covariate_max_pcs", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>SVA</span><span>{analysis_params.get("sva_enabled", "N/A")} ({analysis_params.get("sva_inclusion_rule", "N/A")})</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Clock covariates</span><span>{analysis_params.get("clock_covariates_enabled", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>SV group filter</span><span>P<{analysis_params.get("sv_group_p_threshold", "N/A")}, Eta^2>{analysis_params.get("sv_group_eta2_threshold", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>PVCA min n</span><span>{analysis_params.get("pvca_min_samples", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>DMR (maxgap/p)</span><span>{analysis_params.get("dmr_maxgap", "N/A")} / {analysis_params.get("dmr_p_cutoff", "N/A")}</span></div>\n')
            html_parts.append('        </div>\n')
        if qc_summary:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">QC Summary</div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Samples (pass/fail)</span><span>{qc_summary.get("Samples_passed_QC", "N/A")} / {qc_summary.get("Samples_failed_QC", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Probes final</span><span>{qc_summary.get("Probes_final", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Probes dropped (det/SNP/sex)</span><span>{qc_summary.get("Probes_failed_detection", "N/A")} / {qc_summary.get("Probes_with_SNPs", "N/A")} / {qc_summary.get("Probes_sex_chromosomes", "N/A")}</span></div>\n')
            html_parts.append('        </div>\n')
        html_parts.append('</div>')

    run_docs = [
        ("methods.md", "Methods Summary", "Auto-generated methods text for the run.", "DOC"),
        ("analysis_parameters.json", "Analysis Parameters", "Run settings and thresholds (JSON).", "DOC"),
        ("sessionInfo.txt", "R Session Info", "Full R session and package versions.", "DOC"),
        ("code_version.txt", "Code Version", "Git commit hash when available.", "DOC"),
        ("QC_Summary.csv", "QC Summary (CSV)", "Sample/probe QC counts.", "CSV"),
        ("Sample_QC_Metrics.csv", "Sample QC Metrics (CSV)", "Per-sample QC metrics.", "CSV"),
        ("cell_counts_merged.csv", "Cell Composition (CSV)", "Merged cell composition covariates (if available).", "CSV"),
        ("Cell_Group_Association.csv", "Cell vs Group (CSV)", "Cell composition vs group association.", "CSV"),
    ]
    extra_cell_files = []
    try:
        for fname in os.listdir(output_dir):
            if fname.startswith("cell_counts_") and fname.endswith(".csv") and fname not in {"cell_counts_merged.csv"}:
                extra_cell_files.append(fname)
    except Exception:
        extra_cell_files = []
    for fname in sorted(extra_cell_files):
        run_docs.append((fname, f"{fname}", "Cell composition reference output (CSV).", "CSV"))
    doc_cards = []
    for fname, title, desc, badge in run_docs:
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            rel_path = f"{results_folder_name}/{fname}"
            doc_cards.append((title, desc, rel_path, badge))
    if doc_cards:
        html_parts.append('<div id="docs" class="section-title">Run Documentation</div>')
        html_parts.append('<div class="grid">')
        card_index = 0
        for title, desc, rel_path, badge in doc_cards:
            card_index += 1
            badge_class = f"badge-{badge.lower()}"
            card_html = f"""            <div class="card" style="--i:{card_index};">
                <div class="card-header">
                    {title}
                    <span class="card-badge {badge_class}">{badge}</span>
                </div>
                <div class="card-body">
                    <p class="card-desc">{desc}</p>
                    <a href="{rel_path}" target="_blank" class="btn">Open File</a>
                </div>
            </div>
"""
            html_parts.append(card_html)
        html_parts.append('</div>')

    html_parts.append('<div id="pipelines" class="section-title">Results by Pipeline</div>')
    html_parts.append("""
<div class="nav-tabs" id="navTabs">
    <div class="nav-tab active" data-tab="Intersection" onclick="switchTab('Intersection')">Consensus</div>
    <div class="nav-tab" data-tab="Minfi" onclick="switchTab('Minfi')">Minfi (Noob)</div>
    <div class="nav-tab" data-tab="Sesame" onclick="switchTab('Sesame')">Sesame</div>
</div>
<div class="tab-shell">
""")

    # Loop for Pipelines
    for pipe in pipelines:
        active_class = "active" if pipe == "Intersection" else ""
        html_parts.append(f'    <div id="{pipe}" class="tab-content {active_class}">\n')
        
        # Metrics block (if available)
        metrics = load_metrics(pipe)
        if metrics:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">Model & Batch Summary</div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>λ (inflation)</span><span>{metrics.get("lambda", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Batch method</span><span>{metrics.get("batch_method_applied", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Samples / CpGs</span><span>{metrics.get("n_samples", "N/A")} / {metrics.get("n_cpgs", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Covariates used (n)</span><span>{metrics.get("n_covariates_used", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Covariates used</span><span>{metrics.get("covariates_used", "None") or "None"}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>SVs used</span><span>{metrics.get("n_sv_used", "0")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Batch p<0.05 (Before→After)</span><span>{metrics.get("batch_sig_p_lt_0.05_before", "N/A")} → {metrics.get("batch_sig_p_lt_0.05_after", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Min batch p (Before→After)</span><span>{metrics.get("batch_min_p_before", "N/A")} → {metrics.get("batch_min_p_after", "N/A")}</span></div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>Group min p (Before→After)</span><span>{metrics.get("group_min_p_before", "N/A")} → {metrics.get("group_min_p_after", "N/A")}</span></div>\n')
            if metrics.get("perm_mean_sig") is not None:
                html_parts.append(f'            <div class="metrics-row"><span>Perm mean/max sig (null)</span><span>{metrics.get("perm_mean_sig", "N/A")} / {metrics.get("perm_max_sig", "N/A")}</span></div>\n')
            if metrics.get("vp_primary_group_mean") is not None:
                html_parts.append(f'            <div class="metrics-row"><span>VarPart primary_group</span><span>{metrics.get("vp_primary_group_mean", "N/A")}</span></div>\n')
            if metrics.get("dropped_covariates"):
                html_parts.append(f'            <div class="metrics-row"><span>Dropped covariates</span><span>{metrics.get("dropped_covariates")}</span></div>\n')
            drop_reasons = load_drop_reasons(pipe)
            if drop_reasons:
                reasons_txt = ", ".join([f"{k}:{v}" for k, v in drop_reasons.items()])
                html_parts.append(f'            <div class="metrics-row"><span>Drop reasons</span><span>{reasons_txt}</span></div>\n')
            html_parts.append('        </div>\n')
        files_found = 0
        card_index = 0
        sections = pipeline_sections.get(pipe, [])

        for section_title, section_items in sections:
            section_cards = []
            for suffix, title, desc, badge in section_items:
                filename = f"{pipe}{suffix}"
                file_path = os.path.join(output_dir, filename)
                rel_path = f"{results_folder_name}/{filename}"
                if os.path.exists(file_path):
                    card_index += 1
                    badge_class = f"badge-{badge.lower()}"
                    card_html = f"""            <div class="card" style="--i:{card_index};">
                <div class="card-header">
                    {title}
                    <span class="card-badge {badge_class}">{badge}</span>
                </div>
                <div class="card-body">
                    <p class="card-desc">{desc}</p>
                    <a href="{rel_path}" target="_blank" class="btn">Open Result</a>
                </div>
            </div>
"""
                    section_cards.append(card_html)

            if section_cards:
                files_found += len(section_cards)
                html_parts.append(f'        <div class="section-title">{section_title}</div>\n')
                html_parts.append('        <div class="grid">\n')
                html_parts.extend(section_cards)
                html_parts.append('        </div>\n')

        if files_found == 0:
            html_parts.append('        <div class="empty-msg">No results found for this pipeline.</div>')
             
        html_parts.append('    </div>\n')

    html_parts.append('</div>\n')

    # Footer and Script
    html_parts.append("""
</div>

<script>
    function switchTab(tabName) {
        document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
        document.querySelectorAll('.nav-tab').forEach(el => el.classList.remove('active'));
        document.getElementById(tabName).classList.add('active');
        const match = Array.from(document.querySelectorAll('.nav-tab')).find(el => el.dataset.tab === tabName);
        if (match) match.classList.add('active');
    }
</script>

</body>
</html>
""")

    full_html = "".join(html_parts)

    with open(dashboard_path, "w") as f:
        f.write(full_html)
    
    log(f"[*] Dashboard generated: {dashboard_path}")

def main():
    parser = argparse.ArgumentParser(description="IlluMeta: Automated Illumina DNA Methylation Analysis Pipeline")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Download Command
    parser_download = subparsers.add_parser("download", help="Download GEO data and prepare configuration")
    parser_download.add_argument("gse_id", nargs="?", type=str, help="GEO Series ID (e.g., GSE12345)")
    parser_download.add_argument("-o", "--out-dir", type=str, default=None, help="Output directory (default: ./<GSE_ID>)")
    
    # Analysis Command
    parser_analysis = subparsers.add_parser("analysis", help="Run analysis based on configuration")
    parser_analysis.add_argument("-i", "--input-dir", type=str, help="Project input directory (containing configure.tsv)")
    parser_analysis.add_argument("-c", "--config", type=str, help="Specific path to configure.tsv (overrides --input-dir)")
    parser_analysis.add_argument("-o", "--output", type=str, help="Output directory for results (default: [input-dir]/[test]_vs_[con]_results)")
    parser_analysis.add_argument("--group_con", type=str, required=True, help="Control group label (case-insensitive)")
    parser_analysis.add_argument("--group_test", type=str, required=True, help="Test group label (case-insensitive)")
    parser_analysis.add_argument("--max_plots", type=int, default=10000, help="Max points for interactive plots (default: 10000)")
    parser_analysis.add_argument("--pval", type=float, default=0.05, help="Adjusted P-value threshold (default: 0.05)")
    parser_analysis.add_argument("--lfc", type=float, default=0.5, help="Log2 Fold Change threshold (default: 0.5)")
    parser_analysis.add_argument("--tmp-dir", type=str, help="Custom temporary directory (e.g., external drive)")
    parser_analysis.add_argument("--disable-auto-covariates", action="store_true", help="Disable automatic covariate selection via PCs")
    parser_analysis.add_argument("--disable-sva", action="store_true", help="Disable surrogate variable analysis")
    parser_analysis.add_argument("--include-covariates", type=str, help="Comma-separated covariate names to always try to include (if present in configure.tsv)")
    parser_analysis.add_argument("--include-clock-covariates", action="store_true", help="Compute epigenetic clocks and include them as candidate covariates")
    parser_analysis.add_argument("--tissue", type=str, default="Auto", help="Tissue type for cell deconvolution (default Auto = reference-free; options: Auto, Blood, CordBlood, DLPFC, Placenta)")
    parser_analysis.add_argument("--positive_controls", type=str, help="Comma-separated list of known marker genes to verify (e.g. 'AHRR,CYP1A1')")
    parser_analysis.add_argument("--permutations", type=int, default=0, help="Number of label permutations for null DMP counts (0 to skip)")
    parser_analysis.add_argument("--vp-top", type=int, default=5000, help="Top-variable CpGs for variancePartition (default: 5000)")
    parser_analysis.add_argument("--id-column", type=str, help="Column in configure.tsv to treat as sample ID (useful for non-GEO datasets)")
    parser_analysis.add_argument("--min-total-size", type=int, default=6, help="Minimum total sample size required to proceed (default: 6)")
    parser_analysis.add_argument("--qc-intensity-threshold", type=float, default=10.5, help="Median M/U intensity threshold (log2). Set <=0 to disable intensity-based sample drop (default: 10.5)")
    parser_analysis.add_argument("--force-idat", action="store_true", help="Force reading IDATs if array sizes differ but types are similar")

    # Doctor Command
    parser_doctor = subparsers.add_parser("doctor", help="Check system and R dependencies (does not install)")
    parser_doctor.add_argument("--skip-pandoc", action="store_true", help="Skip pandoc check")
    
    args = parser.parse_args()

    # Early handling for missing subcommand or required args
    if not args.command:
        parser.print_help()
        return
    if args.command == "download" and not args.gse_id:
        parser_download.print_help()
        log("Example: python illumeta.py download GSE12345 -o /path/to/project")
        return

    if args.command == "doctor":
        run_doctor(args)
        return
    
    check_r_installation()
    if args.command in ("download", "analysis"):
        ensure_r_dependencies()
    if args.command == "analysis":
        check_pandoc_installation()
    
    if args.command == "download":
        run_download(args)
    elif args.command == "analysis":
        run_analysis(args)

if __name__ == "__main__":
    main()
