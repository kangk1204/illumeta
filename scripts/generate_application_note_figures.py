#!/usr/bin/env python3
"""
Generate publication-ready figures and tables for Bioinformatics Application Note.
Streamlined to 2 figures (max allowed for Application Note format).
"""
import argparse
import csv
import json
import os
import sys
from pathlib import Path

# Check for required packages
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Installing required packages...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib", "numpy", "-q"])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

# Publication-quality settings
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})


def find_results_dir(gse_dir):
    """Find the results directory inside a GSE folder (handles varied naming)."""
    gse_path = Path(gse_dir)
    # Try standard name first
    std = gse_path / 'illumeta_results'
    if std.exists() and (std / 'summary.json').exists():
        return std
    # Try any directory ending with _results that has summary.json
    for d in sorted(gse_path.iterdir()):
        if d.is_dir() and d.name.endswith('_results') and (d / 'summary.json').exists():
            return d
    return None


def load_summary_tsv(path):
    """Load benchmark summary TSV."""
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def safe_float(val, default=0.0):
    """Safely convert to float."""
    try:
        if val is None or val == '' or val == 'NA':
            return default
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=0):
    """Safely convert to int."""
    try:
        if val is None or val == '' or val == 'NA':
            return default
        return int(float(val))
    except (ValueError, TypeError):
        return default


def create_figure1_workflow(out_dir):
    """
    Figure 1: IlluMeta Workflow Diagram
    """
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis('off')

    # Box style
    box_style = dict(boxstyle='round,pad=0.3', facecolor='lightblue', edgecolor='navy', linewidth=1.5)
    arrow_style = dict(arrowstyle='->', color='navy', lw=2)

    # Input
    ax.text(1, 5, 'Input\n(IDAT/GEO)', ha='center', va='center', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f4e8', edgecolor='green', linewidth=1.5))

    # QC
    ax.annotate('', xy=(2.2, 5), xytext=(1.8, 5), arrowprops=arrow_style)
    ax.text(3, 5, 'Quality\nControl', ha='center', va='center', fontsize=10, fontweight='bold', bbox=box_style)

    # Dual Pipeline
    ax.annotate('', xy=(4.2, 5.3), xytext=(3.8, 5), arrowprops=arrow_style)
    ax.annotate('', xy=(4.2, 4.7), xytext=(3.8, 5), arrowprops=arrow_style)

    ax.text(5.5, 5.5, 'Minfi\n(Noob)', ha='center', va='center', fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff2cc', edgecolor='orange', linewidth=1.5))
    ax.text(5.5, 4.5, 'Sesame\n(pOOBAH)', ha='center', va='center', fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff2cc', edgecolor='orange', linewidth=1.5))

    # Batch Correction
    ax.annotate('', xy=(6.8, 5.5), xytext=(6.2, 5.5), arrowprops=arrow_style)
    ax.annotate('', xy=(6.8, 4.5), xytext=(6.2, 4.5), arrowprops=arrow_style)
    ax.text(7.5, 5, 'Batch\nCorrection\n(SVA/ComBat)', ha='center', va='center', fontsize=9, bbox=box_style)

    # DMP/DMR
    ax.annotate('', xy=(8.5, 5), xytext=(8.2, 5), arrowprops=arrow_style)
    ax.text(9, 5, 'DMP/DMR\nAnalysis', ha='center', va='center', fontsize=9, fontweight='bold', bbox=box_style)

    # Consensus
    ax.annotate('', xy=(7.5, 3.5), xytext=(7.5, 4.2), arrowprops=arrow_style)
    ax.text(7.5, 3, 'Consensus\nIntersection', ha='center', va='center', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffe6e6', edgecolor='red', linewidth=1.5))

    # CRF
    ax.text(5.5, 2, 'CRF\nRobustness\nAssessment', ha='center', va='center', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#e6e6ff', edgecolor='purple', linewidth=1.5))
    ax.annotate('', xy=(6.5, 2.5), xytext=(7, 2.8), arrowprops=dict(arrowstyle='->', color='purple', lw=1.5, ls='--'))

    # Output
    ax.annotate('', xy=(7.5, 1.5), xytext=(7.5, 2.2), arrowprops=arrow_style)
    ax.text(7.5, 1, 'Interactive\nDashboard', ha='center', va='center', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f4e8', edgecolor='green', linewidth=1.5))

    ax.set_title('IlluMeta: Dual-Pipeline DNA Methylation Analysis Workflow', fontsize=12, fontweight='bold', pad=10)

    plt.savefig(os.path.join(out_dir, 'Figure1_Workflow.png'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(out_dir, 'Figure1_Workflow.pdf'), facecolor='white')
    plt.close()
    print("  Created Figure 1: Workflow diagram")


def create_figure2_benchmark(summary_data, ablation_data, out_dir):
    """
    Figure 2: Compact 2x2 multi-panel benchmark summary.
    A) DMP counts by pipeline  B) Lambda (genomic inflation)
    C) Consensus overlap        D) CRF tier heatmap
    """
    if not summary_data:
        print("  Skipping Figure 2: No summary data")
        return

    # Sort datasets by platform order (450k → EPIC → EPICv2), then by sample size
    platform_order = {'450k': 0, 'EPIC': 1, 'EPICv2': 2}
    summary_data = sorted(summary_data,
                          key=lambda r: (platform_order.get(r.get('platform', ''), 9),
                                         safe_int(r.get('n_samples', 0))))

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5))

    datasets = [row.get('gse_id', row.get('label', 'Unknown')) for row in summary_data]
    n_samples = [safe_int(row.get('n_samples', 0)) for row in summary_data]
    crf_tiers = [row.get('crf_tier', '').capitalize() for row in summary_data]
    platforms = [row.get('platform', '') for row in summary_data]

    minfi_dmps = [safe_int(row.get('minfi_dmps_fdr', 0)) for row in summary_data]
    sesame_dmps = [safe_int(row.get('sesame_dmps_fdr', 0)) for row in summary_data]
    consensus_dmps = [safe_int(row.get('illumeta_intersect_dmps', 0)) for row in summary_data]

    lambda_minfi = [safe_float(row.get('lambda_minfi', 1.0), 1.0) for row in summary_data]
    lambda_sesame = [safe_float(row.get('lambda_sesame', 1.0), 1.0) for row in summary_data]

    x = np.arange(len(datasets))
    width = 0.22
    short_labels = [d.replace('GSE', '') for d in datasets]

    # Platform-colored backgrounds
    platform_colors = {'450k': '#f0f0f0', 'EPIC': '#e8f0ff', 'EPICv2': '#fff0e8'}

    def add_platform_bg(ax):
        """Add subtle platform-group shading."""
        for i, plat in enumerate(platforms):
            ax.axvspan(i - 0.45, i + 0.45, alpha=0.15,
                       color=platform_colors.get(plat, '#f0f0f0'), zorder=0)

    # ── Panel A: DMP counts by pipeline ──
    ax = axes[0, 0]
    add_platform_bg(ax)
    ax.bar(x - width, minfi_dmps, width, label='Minfi', color='#4ECDC4', edgecolor='black', linewidth=0.5, zorder=3)
    ax.bar(x, sesame_dmps, width, label='SeSAMe', color='#FF6B6B', edgecolor='black', linewidth=0.5, zorder=3)
    ax.bar(x + width, consensus_dmps, width, label='Consensus', color='#45B7D1', edgecolor='black', linewidth=0.5, zorder=3)
    ax.set_ylabel('DMPs (FDR < 0.05)')
    ax.set_title('A. DMP Detection by Pipeline', fontweight='bold', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n{p}, n={n}" for d, p, n in zip(short_labels, platforms, n_samples)],
                        fontsize=5.5, rotation=0)
    ax.legend(loc='upper left', fontsize=7, framealpha=0.9)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))
    ax.grid(axis='y', alpha=0.3, zorder=0)

    # ── Panel B: Lambda (genomic inflation) — capped at 2.5 with annotations ──
    ax = axes[0, 1]
    add_platform_bg(ax)
    LAMBDA_CAP = 2.5
    lm_plot = [min(v, LAMBDA_CAP) for v in lambda_minfi]
    ls_plot = [min(v, LAMBDA_CAP) for v in lambda_sesame]
    bars_m = ax.bar(x - width/2, lm_plot, width, label='Minfi', color='#4ECDC4', edgecolor='black', linewidth=0.5, zorder=3)
    bars_s = ax.bar(x + width/2, ls_plot, width, label='SeSAMe', color='#FF6B6B', edgecolor='black', linewidth=0.5, zorder=3)
    ax.axhline(y=1.0, color='green', linestyle='--', linewidth=1.5, label='Ideal ($\\lambda$=1.0)', zorder=2)
    ax.axhspan(0.9, 1.1, alpha=0.12, color='green', zorder=1)
    # Annotate values that exceed cap
    for i, (lm, ls) in enumerate(zip(lambda_minfi, lambda_sesame)):
        if lm > LAMBDA_CAP:
            ax.annotate(f'{lm:.1f}', xy=(i - width/2, LAMBDA_CAP), xytext=(-8, 5),
                        textcoords='offset points', ha='center', fontsize=5.5, fontweight='bold',
                        color='#2AA69A')
        if ls > LAMBDA_CAP:
            ax.annotate(f'{ls:.1f}', xy=(i + width/2, LAMBDA_CAP), xytext=(8, 5),
                        textcoords='offset points', ha='center', fontsize=5.5, fontweight='bold',
                        color='#E05555')
    ax.set_ylabel('Genomic Inflation ($\\lambda$)')
    ax.set_title('B. Lambda After Batch Correction', fontweight='bold', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n({t})" for d, t in zip(short_labels, crf_tiers)], fontsize=6)
    ax.legend(loc='upper right', fontsize=6.5, framealpha=0.9)
    ax.set_ylim(0, LAMBDA_CAP + 0.3)
    ax.grid(axis='y', alpha=0.3, zorder=0)

    # ── Panel C: Consensus overlap with pipeline concordance ──
    ax = axes[1, 0]
    add_platform_bg(ax)
    specificity = [c / min(m, s) * 100 if min(m, s) > 0 else 0
                   for m, s, c in zip(minfi_dmps, sesame_dmps, consensus_dmps)]
    colors_c = ['#45B7D1' if c > 0 else '#cccccc' for c in consensus_dmps]
    bars = ax.bar(x, [c / 1000 for c in consensus_dmps], width * 2.5,
                  color=colors_c, edgecolor='black', linewidth=0.5, zorder=3)
    for i, (bar, spec) in enumerate(zip(bars, specificity)):
        height = bar.get_height()
        if consensus_dmps[i] > 0:
            ax.annotate(f'{consensus_dmps[i]:,}\n({spec:.0f}%)',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=7, fontweight='bold')
        else:
            ax.annotate('0 DMPs',
                        xy=(bar.get_x() + bar.get_width() / 2, 0.3),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=7, color='gray')
    ax.set_ylabel('Consensus DMPs (x1,000)')
    ax.set_title('C. Dual-Pipeline Consensus', fontweight='bold', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n{p}, n={n}" for d, p, n in zip(short_labels, platforms, n_samples)],
                        fontsize=5.5, rotation=0)
    ax.grid(axis='y', alpha=0.3, zorder=0)

    # ── Panel D: CRF tier heatmap ──
    ax = axes[1, 1]
    tiers = ['Minimal\n(n<12)', 'Small\n(12\u201323)', 'Moderate\n(24\u201349)', 'Large\n(n\u226550)']
    features = ['SVA', 'ComBat', 'Full MMC', 'Full NCS', 'Full SSS']
    availability = np.array([
        [0, 0, 0.5, 0.5, 0.5],   # Minimal
        [0, 0.5, 0.5, 0.5, 0.5], # Small
        [1, 1, 1, 1, 1],         # Moderate
        [1, 1, 1, 1, 1],         # Large
    ])
    tier_colors = ['#ff6b6b', '#ffd93d', '#6bcb77']
    cmap = plt.cm.colors.ListedColormap(tier_colors)
    im = ax.imshow(availability.T, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(tiers)))
    ax.set_yticks(np.arange(len(features)))
    ax.set_xticklabels(tiers, fontsize=7)
    ax.set_yticklabels(features, fontsize=8)
    for i in range(len(features)):
        for j in range(len(tiers)):
            val = availability[j, i]
            text = 'Y' if val == 1 else ('\u2248' if val == 0.5 else 'N')
            color = 'white' if val == 0 else 'black'
            ax.text(j, i, text, ha='center', va='center', fontsize=10, color=color, fontweight='bold')
    ax.set_title('D. CRF Tier Features', fontweight='bold', fontsize=10)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#6bcb77', label='Enabled'),
        Patch(facecolor='#ffd93d', label='Limited'),
        Patch(facecolor='#ff6b6b', label='Disabled'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=6.5, ncol=3)

    plt.tight_layout(h_pad=2.0, w_pad=1.5)
    plt.savefig(os.path.join(out_dir, 'Figure2_Benchmark.png'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(out_dir, 'Figure2_Benchmark.pdf'), facecolor='white')
    plt.close()
    print("  Created Figure 2: Benchmark multi-panel")


def create_table1_features(out_dir):
    """Table 1: Feature comparison between tools."""
    headers = ['Feature', 'minfi', 'ChAMP', 'IlluMeta']
    rows = [
        ['Normalization', 'Noob/Funnorm', 'BMIQ', 'Noob + pOOBAH'],
        ['Dual Pipeline', 'No', 'No', 'Yes'],
        ['Consensus Calling', 'No', 'No', 'Yes'],
        ['Automated Batch Correction', 'No', 'ComBat only', 'SVA/ComBat/limma (auto)'],
        ['Sample Size Adaptation', 'No', 'No', 'CRF'],
        ['Interactive Dashboard', 'No', 'No', 'Yes'],
        ['GEO Integration', 'Manual', 'Manual', 'Built-in download'],
        ['Auto-generated Methods', 'No', 'No', 'Yes'],
        ['Cell Composition', 'FlowSorted', 'RefbaseEWAS', 'EpiDISH + RefFreeEWAS'],
        ['Epigenetic Clocks', 'No', 'No', 'methylclock/planet'],
        ['EPIC v2 Support', 'Manual setup', 'No', 'Auto-detected'],
    ]

    with open(os.path.join(out_dir, 'Table1_Feature_Comparison.tsv'), 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')

    print("  Created Table 1: Feature comparison")


def create_table2_benchmark_summary(summary_data, out_dir):
    """Table 2: Benchmark summary across datasets."""
    if not summary_data:
        print("  Skipping Table 2: No summary data")
        return

    # Sort by platform then sample size (matching Figure 2 order)
    platform_order = {'450k': 0, 'EPIC': 1, 'EPICv2': 2}
    summary_data = sorted(summary_data,
                          key=lambda r: (platform_order.get(r.get('platform', ''), 9),
                                         safe_int(r.get('n_samples', 0))))

    headers = ['Dataset', 'Platform', 'n', 'CRF Tier',
               'Minfi DMPs (FDR)', 'Sesame DMPs (FDR)', 'Consensus DMPs',
               'Minfi DMPs (Nominal)', 'Sesame DMPs (Nominal)',
               'Minfi DMRs', 'Sesame DMRs',
               'Lambda (Minfi)', 'Lambda (Sesame)', 'CAF Score']

    rows = []
    for row in summary_data:
        rows.append([
            row.get('gse_id', ''),
            row.get('platform', ''),
            row.get('n_samples', ''),
            row.get('crf_tier', '').capitalize(),
            f"{safe_int(row.get('minfi_dmps_fdr', 0)):,}",
            f"{safe_int(row.get('sesame_dmps_fdr', 0)):,}",
            f"{safe_int(row.get('illumeta_intersect_dmps', 0)):,}",
            f"{safe_int(row.get('minfi_dmps_nominal', 0)):,}",
            f"{safe_int(row.get('sesame_dmps_nominal', 0)):,}",
            f"{safe_int(row.get('minfi_dmrs', 0)):,}",
            f"{safe_int(row.get('sesame_dmrs', 0)):,}",
            f"{safe_float(row.get('lambda_minfi', 1.0)):.3f}",
            f"{safe_float(row.get('lambda_sesame', 1.0)):.3f}",
            f"{safe_float(row.get('caf_score', 0)):.4f}",
        ])

    with open(os.path.join(out_dir, 'Table2_Benchmark_Summary.tsv'), 'w', encoding='utf-8') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')

    print("  Created Table 2: Benchmark summary")


def collect_illumeta_results(results_dir):
    """Collect metrics from IlluMeta results directories."""
    summary = []
    results_path = Path(results_dir)

    for gse_dir in sorted(results_path.glob('GSE*')):
        if not gse_dir.is_dir():
            continue

        res_dir = find_results_dir(gse_dir)
        if res_dir is None:
            print(f"  Skipping {gse_dir.name}: no results directory found")
            continue

        gse_id = gse_dir.name
        row = {'gse_id': gse_id}

        # Load summary.json for basic info
        summary_file = res_dir / 'summary.json'
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                sj = json.load(f)
                row['n_con'] = sj.get('n_con', 0)
                row['n_test'] = sj.get('n_test', 0)
                row['caf_score'] = sj.get('primary_caf_score', 0)

        # Get sample count and tier from CRF_Sample_Tier.csv
        tier_file = res_dir / 'CRF_Sample_Tier.csv'
        if tier_file.exists():
            with open(tier_file, 'r') as f:
                reader = csv.DictReader(f)
                for r in reader:
                    row['n_samples'] = r.get('total_n', '0')
                    row['crf_tier'] = r.get('tier', '')
                    break

        # Detect platform from Probe_Filter_Summary.csv (raw probe count)
        probe_file = res_dir / 'Probe_Filter_Summary.csv'
        if probe_file.exists():
            with open(probe_file, 'r') as f:
                reader = csv.DictReader(f)
                for r in reader:
                    if r.get('step') == 'raw':
                        n_probes = safe_int(r.get('remaining', 0))
                        if n_probes > 900000:
                            row['platform'] = 'EPICv2'
                        elif n_probes > 600000:
                            row['platform'] = 'EPIC'
                        else:
                            row['platform'] = '450k'
                        break

        # Get group info
        group_file = res_dir / 'Input_Group_Distribution.csv'
        if group_file.exists():
            with open(group_file, 'r') as f:
                reader = csv.DictReader(f)
                groups = list(reader)
                if len(groups) >= 2:
                    row['group_con'] = groups[0].get('group', 'Control')[:25]
                    row['group_test'] = groups[1].get('group', 'Test')[:25]

            # Get n_samples from summary.json (post-QC count)
        if 'n_con' in row and 'n_test' in row:
            n_total = row['n_con'] + row['n_test']
            row['n_samples'] = str(n_total)

        # Get lambda from ablation summary
        for pipeline in ['Minfi', 'Sesame']:
            abl_file = res_dir / f'{pipeline}_Ablation_Summary.csv'
            if abl_file.exists():
                with open(abl_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        if r.get('metric') == 'lambda':
                            row[f'lambda_{pipeline.lower()}'] = r.get('corrected', r.get('raw', '1.0'))
                            break

        # Count DMPs (both FDR<0.05 and nominal: raw P<0.05 + |logFC|>1)
        for pipeline in ['Minfi', 'Sesame']:
            dmp_file = res_dir / f'{pipeline}_DMPs_full.csv'
            if dmp_file.exists():
                fdr_count = 0
                nominal_count = 0
                with open(dmp_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        try:
                            fdr = float(r.get('adj.P.Val', r.get('FDR', '1')))
                            if fdr < 0.05:
                                fdr_count += 1
                            pval = float(r.get('P.Value', r.get('pval', '1')))
                            logfc = abs(float(r.get('logFC', '0')))
                            if pval < 0.05 and logfc > 1:
                                nominal_count += 1
                        except (ValueError, TypeError):
                            pass
                row[f'{pipeline.lower()}_dmps_fdr'] = fdr_count
                row[f'{pipeline.lower()}_dmps_nominal'] = nominal_count

        # Consensus DMPs
        consensus_file = res_dir / 'Intersection_Consensus_DMPs.csv'
        if consensus_file.exists():
            with open(consensus_file, 'r') as f:
                row['illumeta_intersect_dmps'] = sum(1 for _ in f) - 1  # subtract header

        # Count DMRs per pipeline
        for pipeline in ['Minfi', 'Sesame']:
            dmr_file = res_dir / f'{pipeline}_DMRs.csv'
            if dmr_file.exists():
                with open(dmr_file, 'r') as f:
                    row[f'{pipeline.lower()}_dmrs'] = sum(1 for _ in f) - 1  # subtract header

        if row.get('n_samples'):
            summary.append(row)
            print(f"  Collected: {gse_id} (n={row.get('n_samples')}, tier={row.get('crf_tier', 'N/A')}, "
                  f"results_dir={res_dir.name})")

    return summary


def collect_ablation_long(results_dir):
    """Collect ablation data in long format for visualization."""
    ablation = []
    results_path = Path(results_dir)

    for gse_dir in sorted(results_path.glob('GSE*')):
        if not gse_dir.is_dir():
            continue

        res_dir = find_results_dir(gse_dir)
        if res_dir is None:
            continue

        gse_id = gse_dir.name

        for pipeline in ['Minfi', 'Sesame']:
            abl_file = res_dir / f'{pipeline}_Ablation_Summary.csv'
            if not abl_file.exists():
                continue

            with open(abl_file, 'r') as f:
                reader = csv.DictReader(f)
                for r in reader:
                    metric = r.get('metric', '')
                    raw_val = r.get('raw', '')
                    corr_val = r.get('corrected', '')

                    if raw_val:
                        ablation.append({
                            'gse_id': gse_id,
                            'pipeline': pipeline,
                            'variant': 'baseline',
                            'metric': metric,
                            'value': raw_val
                        })
                    if corr_val:
                        ablation.append({
                            'gse_id': gse_id,
                            'pipeline': pipeline,
                            'variant': 'corrected',
                            'metric': metric,
                            'value': corr_val
                        })

    return ablation


def main():
    parser = argparse.ArgumentParser(description='Generate Application Note figures and tables')
    parser.add_argument('--summary-tsv', default=None,
                        help='Path to benchmark summary TSV (optional)')
    parser.add_argument('--results-dir', default='benchmarks/application_note_rerun',
                        help='Directory containing IlluMeta results (GSE* subdirectories)')
    parser.add_argument('--out-dir', default='benchmarks/paper_figures',
                        help='Output directory for figures and tables')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"\nGenerating Application Note figures and tables...")
    print(f"Output directory: {args.out_dir}\n")

    # Load data from summary TSV if provided, otherwise collect from results
    summary_data = []
    if args.summary_tsv and os.path.exists(args.summary_tsv):
        summary_data = load_summary_tsv(args.summary_tsv)
        print(f"Loaded {len(summary_data)} datasets from summary TSV")
    elif os.path.exists(args.results_dir):
        print(f"Collecting metrics from: {args.results_dir}")
        summary_data = collect_illumeta_results(args.results_dir)
        print(f"Collected {len(summary_data)} datasets from results directories")

    ablation_data = []
    if os.path.exists(args.results_dir):
        ablation_data = collect_ablation_long(args.results_dir)
        print(f"Collected {len(ablation_data)} ablation metrics from results")

    print("\nGenerating figures (2 figures for Application Note)...")

    # Figure 1: Workflow
    create_figure1_workflow(args.out_dir)

    # Figure 2: Compact multi-panel benchmark summary
    create_figure2_benchmark(summary_data, ablation_data, args.out_dir)

    print("\nGenerating tables...")
    create_table1_features(args.out_dir)
    create_table2_benchmark_summary(summary_data, args.out_dir)

    print(f"\nDone! Files saved to: {args.out_dir}/")
    print("\nGenerated files:")
    for f in sorted(os.listdir(args.out_dir)):
        size = os.path.getsize(os.path.join(args.out_dir, f))
        print(f"  - {f} ({size:,} bytes)")


if __name__ == '__main__':
    main()
