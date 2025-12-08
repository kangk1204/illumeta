#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import csv
import json
from datetime import datetime

# Configuration for R script paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
R_SCRIPTS_DIR = os.path.join(BASE_DIR, "r_scripts")
DOWNLOAD_SCRIPT = os.path.join(R_SCRIPTS_DIR, "download.R")
ANALYZE_SCRIPT = os.path.join(R_SCRIPTS_DIR, "analyze.R")
SETUP_MARKER = os.path.join(BASE_DIR, ".r_setup_done")
SETUP_SCRIPT = os.path.join(R_SCRIPTS_DIR, "setup_env.R")
DEFAULT_R_LIB = os.path.join(BASE_DIR, ".r-lib")
def add_conda_paths(env: dict) -> dict:
    """Ensure LD_LIBRARY_PATH/PKG_CONFIG_PATH/PATH include conda libs so xml2/xml load correctly."""
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

def missing_r_packages(pkgs):
    """Returns a list of missing R packages using requireNamespace checks."""
    if not pkgs:
        return []
    pkg_str = '", "'.join(pkgs)
    expr = (
        f'pkgs <- c("{pkg_str}"); '
        'missing <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly=TRUE))]; '
        'if (length(missing) > 0) { cat(paste(missing, collapse=",")); quit(status=1) }'
    )
    env = ensure_r_lib_env(os.environ.copy())
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
        # R not found; handled separately in check_r_installation
        return []
    if res.returncode == 0:
        return []
    output = (res.stdout + res.stderr).strip()
    if not output:
        return pkgs
    return [p.strip() for p in output.replace("\n", " ").split(",") if p.strip()]

def ensure_r_lib_env(env):
    """Ensure R_LIBS_USER points to a repo-local, writable library."""
    if env.get("R_LIBS_USER"):
        return env
    env["R_LIBS_USER"] = DEFAULT_R_LIB
    if not os.path.exists(DEFAULT_R_LIB):
        os.makedirs(DEFAULT_R_LIB, exist_ok=True)
        log(f"[*] Created default R library at {DEFAULT_R_LIB}")
    return env

def ensure_r_dependencies():
    """Runs the setup_env.R script to install required R packages and cache data."""
    env = ensure_r_lib_env(os.environ.copy())
    env = add_conda_paths(env)
    force_setup = os.environ.get("ILLUMETA_FORCE_SETUP") == "1"
    marker_ok = False
    if os.path.exists(SETUP_MARKER):
        try:
            with open(SETUP_MARKER) as f:
                content = f.read()
            marker_ok = "epicv2_required=1" in content
        except Exception:
            marker_ok = False
    core_pkgs = ["optparse", "GEOquery", "sesame", "minfi"]
    missing_core = missing_r_packages(core_pkgs) if marker_ok else []
    if marker_ok and not force_setup and not missing_core:
        log("[*] R dependencies already set up (skipping). Set ILLUMETA_FORCE_SETUP=1 to force reinstall.")
        return
    if missing_core and not force_setup:
        log(f"[*] R setup marker found but missing packages: {', '.join(missing_core)}. Re-running setup_env.R ...")

    log("[*] Ensuring R dependencies (this may take a few minutes on first run)...")
    try:
        subprocess.run(["Rscript", SETUP_SCRIPT], check=True, env=env)
        with open(SETUP_MARKER, "w") as f:
            f.write(f"setup completed at {datetime.now().isoformat()}\n")
            f.write("epicv2_required=1\n")
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
        # Default: test_vs_con_results
        dir_base = args.input_dir if args.input_dir else os.path.dirname(os.path.abspath(config_path))
        folder_name = f"{args.group_test}_vs_{args.group_con}_results"
        output_dir = os.path.join(dir_base, folder_name)
    
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
        "--lfc", str(args.lfc)
    ]
    if args.disable_auto_covariates:
        cmd.append("--disable_auto_covariates")
    if args.disable_sva:
        cmd.append("--disable_sva")
    if args.include_covariates:
        cmd.extend(["--include_covariates", args.include_covariates])
    if args.permutations and args.permutations > 0:
        cmd.extend(["--permutations", str(args.permutations)])
    if args.vp_top:
        cmd.extend(["--vp_top", str(args.vp_top)])
    
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
    dashboard_filename = f"{group_test}_vs_{group_con}_results_index.html"
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

    # Define files always shown for all pipelines (or for Intersection)
    common_files_map = [
        ("_Volcano.html", "Volcano Plot", "Visualizes significant changes (LogFC vs P-value)."),
        ("_Manhattan.html", "Manhattan Plot", "Genomic distribution of methylation changes."),
        ("_Top100_Heatmap.html", "Heatmap (Top 100)", "Clustering of the top 100 most significant CpGs."),
        ("_Top_DMPs.html", "Top DMPs Table", "Searchable table of differentially methylated probes (sorted by P.Value)."),
        ("_DMR_Volcano.html", "DMR Volcano Plot", "Visualizes significant DMRs (Est. Diff vs P-value)."),
        ("_DMR_Manhattan.html", "DMR Manhattan Plot", "Genomic distribution of DMRs."),
        ("_Top_DMRs_Heatmap.html", "DMR Heatmap (Top 50)", "Average methylation levels of top 50 DMRs."),
        ("_DMRs_Table.html", "DMRs Table", "Searchable table of differentially methylated regions."),
        ("_PVCA.html", "PVCA (Raw)", "Variance explained by group/batch/covariates before correction."),
        ("_AfterCorrection_PVCA.html", "PVCA (Corrected)", "Variance explained after batch correction/SVA."),
        ("_Epigenetic_Age.html", "Epigenetic Age (Clocks)", "Horvath/Hannum/PhenoAge vs chronological age (if provided)."),
        ("_Epigenetic_Age.csv", "Epigenetic Age Table", "DNAm clock predictions per sample (predicted age, missing probes, residuals).")
    ]

    # Define additional files for Minfi/Sesame (full pipeline)
    full_pipeline_files_map = common_files_map + [
        ("_PCA_Before.html", "PCA (Raw)", "Principal Component Analysis before batch correction."),
        ("_PCA_After_Correction.html", "PCA (Corrected)", "PCA after removing detected batch effects."),
        ("_Sample_Clustering_Distance.html", "Sample Clustering", "Euclidean distance matrix of samples."),
        ("_QQPlot.html", "Q-Q Plot", "Check for genomic inflation and systematic bias."),
        ("_Batch_Evaluation_Before.html", "Batch Eval (Before)", "Covariate association with PCs (Raw)."),
        ("_Batch_Evaluation_After.html", "Batch Eval (After)", "Covariate association with PCs (Corrected).")
    ]
    
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

    # CSS Style Block (Using format to avoid curly brace hell)
    style_block = """
        :root { --primary: #2c3e50; --secondary: #34495e; --accent: #3498db; --bg: #f4f7f6; --text: #333; --success: #27ae60; --danger: #e74c3c; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; background-color: var(--bg); color: var(--text); }
        header { background-color: var(--primary); color: white; padding: 2rem; text-align: center; }
        header h1 { margin: 0; font-size: 2rem; }
        header p { margin: 0.5rem 0 0; opacity: 0.8; }
        
        .nav-tabs { display: flex; justify-content: center; background: white; padding: 1rem 0; box-shadow: 0 2px 5px rgba(0,0,0,0.05); position: sticky; top: 0; z-index: 100; }
        .nav-tab { padding: 0.8rem 2rem; cursor: pointer; font-weight: 600; color: var(--secondary); border-bottom: 3px solid transparent; transition: all 0.3s; margin: 0 0.5rem; }
        .nav-tab:hover { background-color: #f0f0f0; border-radius: 5px 5px 0 0; }
        .nav-tab.active { border-bottom: 3px solid var(--accent); color: var(--accent); }
        
        .container { max-width: 1200px; margin: 2rem auto; padding: 0 20px; }
        
        /* Summary Section */
        .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .stat-card { background: white; padding: 1.5rem; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); text-align: center; border-top: 4px solid var(--secondary); }
        .stat-card.highlight { border-top: 4px solid var(--success); transform: scale(1.02); }
        .stat-title { font-size: 0.9rem; color: #777; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 10px; }
        .stat-value { font-size: 2rem; font-weight: bold; color: var(--primary); }
        .stat-sub { font-size: 0.9rem; margin-top: 5px; }
        .up { color: var(--danger); } .down { color: var(--accent); } 
        
        .tab-content { display: none; animation: fadeIn 0.5s; }
        .tab-content.active { display: block; }
        
        .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 25px; }
        .card { background: white; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); overflow: hidden; transition: transform 0.2s, box-shadow 0.2s; display: flex; flex-direction: column; }
        .card:hover { transform: translateY(-5px); box-shadow: 0 8px 25px rgba(0,0,0,0.1); }
        .card-header { background: var(--primary); color: white; padding: 1rem; font-weight: bold; display: flex; justify-content: space-between; align-items: center; }
        .card-badge { font-size: 0.7rem; background: rgba(255,255,255,0.2); padding: 2px 8px; border-radius: 10px; }
        .card-body { padding: 1.5rem; flex-grow: 1; }
        .card-desc { font-size: 0.9rem; color: #666; line-height: 1.5; margin-bottom: 1.5rem; }
        .btn { display: block; width: 100%; padding: 10px 0; background-color: var(--accent); color: white; text-align: center; text-decoration: none; border-radius: 6px; font-weight: 600; transition: background 0.2s; }
        .btn:hover { background-color: #2980b9; }
        
        .empty-msg { text-align: center; padding: 3rem; color: #999; font-style: italic; }
        
        .metrics-card { background: #fdfdfd; border: 1px solid #eee; border-radius: 8px; padding: 1rem; margin-bottom: 1rem; }
        .metrics-title { font-weight: 700; margin-bottom: 0.5rem; color: var(--secondary); }
        .metrics-row { display: flex; justify-content: space-between; font-size: 0.95rem; padding: 4px 0; border-bottom: 1px dashed #eee; }
        .metrics-row:last-child { border-bottom: none; }

        @keyframes fadeIn { from { opacity: 0; transform: translateY(10px); } to { opacity: 1; transform: translateY(0); } }
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
    <h1>Analysis Dashboard</h1>
    <p>Comparison: <strong>{group_test}</strong> (Test) vs <strong>{group_con}</strong> (Control)</p>
</header>

<div class="nav-tabs" id="navTabs">
    <div class="nav-tab active" onclick="switchTab('Intersection')">Intersection (Consensus)</div>
    <div class="nav-tab" onclick="switchTab('Minfi')">Minfi (Noob)</div>
    <div class="nav-tab" onclick="switchTab('Sesame')">Sesame</div>
</div>

<div class="container">

    <!-- Summary Statistics -->
    <div class="summary-grid">
        <div class="stat-card">
            <div class="stat-title">Samples</div>
            <div class="stat-value">{total_samples}</div>
            <div class="stat-sub">Con: {n_con}, Test: {n_test}</div>
        </div>
        <div class="stat-card">
            <div class="stat-title">Minfi DMPs</div>
            <div class="stat-value">{minfi_total}</div>
            <div class="stat-sub"><span class="up">▲ {minfi_up}</span> <span class="down">▼ {minfi_down}</span></div>
        </div>
        <div class="stat-card">
            <div class="stat-title">Sesame DMPs</div>
            <div class="stat-value">{sesame_total}</div>
            <div class="stat-sub"><span class="up">▲ {sesame_up}</span> <span class="down">▼ {sesame_down}</span></div>
        </div>
        <div class="stat-card highlight">
            <div class="stat-title">Intersection DMPs</div>
            <div class="stat-value">{intersect_total}</div>
            <div class="stat-sub"><span class="up">▲ {intersect_up}</span> <span class="down">▼ {intersect_down}</span></div>
        </div>
    </div>
""")
    
    # Loop for Pipelines
    for pipe in pipelines:
        active_class = "active" if pipe == "Intersection" else ""
        html_parts.append(f'    <div id="{pipe}" class="tab-content {active_class}">\n')
        
        # Metrics block (if available)
        metrics = load_metrics(pipe)
        if metrics:
            html_parts.append('        <div class="metrics-card">\n')
            html_parts.append('            <div class="metrics-title">QC Metrics</div>\n')
            html_parts.append(f'            <div class="metrics-row"><span>λ (inflation)</span><span>{metrics.get("lambda", "N/A")}</span></div>\n')
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
            html_parts.append('        </div>\n')
        html_parts.append('        <div class="grid">\n')
        
        files_found = 0
        current_file_map = full_pipeline_files_map if pipe != "Intersection" else common_files_map

        for suffix, title, desc in current_file_map:
            filename = f"{pipe}{suffix}"
            file_path = os.path.join(output_dir, filename)
            rel_path = f"{results_folder_name}/{filename}"
            
            if os.path.exists(file_path):
                files_found += 1
                badge = "TABLE" if "Table" in title or "DMPs" in title else "PLOT"
                
                # Card HTML
                card_html = f"""            <div class="card">
                <div class="card-header">
                    {title}
                    <span class="card-badge">{badge}</span>
                </div>
                <div class="card-body">
                    <p class="card-desc">{desc}</p>
                    <a href="{rel_path}" target="_blank" class="btn">Open Result</a>
                </div>
            </div>
"""
                html_parts.append(card_html)
        
        if files_found == 0:
             html_parts.append('            <div class="empty-msg">No results found for this pipeline.</div>')
             
        html_parts.append('        </div>\n    </div>\n')

    # Footer and Script
    html_parts.append("""
</div>

<script>
    function switchTab(tabName) {
        // Hide all contents
        document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
        // Remove active class from tabs
        document.querySelectorAll('.nav-tab').forEach(el => el.classList.remove('active'));
        
        // Show selected
        document.getElementById(tabName).classList.add('active');
        // Highlight tab
        Array.from(document.querySelectorAll('.nav-tab')).find(el => el.textContent.includes(tabName)).classList.add('active');
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
    parser_analysis.add_argument("--permutations", type=int, default=0, help="Number of label permutations for null DMP counts (0 to skip)")
    parser_analysis.add_argument("--vp-top", type=int, default=5000, help="Top-variable CpGs for variancePartition (default: 5000)")
    
    args = parser.parse_args()

    # Early handling for missing subcommand or required args
    if not args.command:
        parser.print_help()
        return
    if args.command == "download" and not args.gse_id:
        parser_download.print_help()
        log("Example: python illumeta.py download GSE12345 -o /path/to/project")
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
