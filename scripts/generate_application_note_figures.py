#!/usr/bin/env python3
"""
Generate publication-ready figures and tables for Bioinformatics Application Note.
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


def load_summary_tsv(path):
    """Load benchmark summary TSV."""
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def load_ablation_metrics(path):
    """Load ablation metrics TSV."""
    rows = []
    if not os.path.exists(path):
        return rows
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
    ax.text(5.5, 2, 'CRF v2.1\nRobustness\nAssessment', ha='center', va='center', fontsize=9,
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


def create_figure2_benchmark_comparison(summary_data, out_dir):
    """
    Figure 2: IlluMeta dual-pipeline comparison (Minfi vs Sesame, Lambda, Consensus)
    """
    if not summary_data:
        print("  Skipping Figure 2: No summary data")
        return

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))

    # Prepare data
    datasets = [row.get('gse_id', row.get('label', 'Unknown')) for row in summary_data]
    n_samples = [safe_int(row.get('n_samples', 0)) for row in summary_data]
    crf_tiers = [row.get('crf_tier', '').capitalize() for row in summary_data]

    minfi_dmps = [safe_int(row.get('minfi_dmps_fdr', 0)) for row in summary_data]
    sesame_dmps = [safe_int(row.get('sesame_dmps_fdr', 0)) for row in summary_data]
    consensus_dmps = [safe_int(row.get('illumeta_intersect_dmps', 0)) for row in summary_data]

    lambda_minfi = [safe_float(row.get('lambda_minfi', 1.0), 1.0) for row in summary_data]
    lambda_sesame = [safe_float(row.get('lambda_sesame', 1.0), 1.0) for row in summary_data]

    x = np.arange(len(datasets))
    width = 0.35

    # Panel A: DMP counts by pipeline
    ax = axes[0]
    bars1 = ax.bar(x - width/2, minfi_dmps, width, label='Minfi (Noob)', color='#4ECDC4', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, sesame_dmps, width, label='Sesame (pOOBAH)', color='#FF6B6B', edgecolor='black', linewidth=0.5)
    ax.set_ylabel('Significant DMPs (FDR < 0.05)')
    ax.set_title('A. DMP Detection by Pipeline', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n(n={n})" for d, n in zip(datasets, n_samples)], fontsize=8)
    ax.legend(loc='upper right', fontsize=8)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # Panel B: Lambda (genomic inflation)
    ax = axes[1]
    ax.bar(x - width/2, lambda_minfi, width, label='Minfi', color='#4ECDC4', edgecolor='black', linewidth=0.5)
    ax.bar(x + width/2, lambda_sesame, width, label='Sesame', color='#FF6B6B', edgecolor='black', linewidth=0.5)
    ax.axhline(y=1.0, color='green', linestyle='--', linewidth=1.5, label='Ideal (λ=1.0)')
    ax.axhspan(0.9, 1.1, alpha=0.1, color='green')
    ax.set_ylabel('Genomic Inflation (λ)')
    ax.set_title('B. Lambda After Correction', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n({t})" for d, t in zip(datasets, crf_tiers)], fontsize=8)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_ylim(0.5, 2.0)

    # Panel C: Consensus approach benefit
    ax = axes[2]
    union_dmps = [m + s - c for m, s, c in zip(minfi_dmps, sesame_dmps, consensus_dmps)]
    specificity = [c / (m + s) * 100 if (m + s) > 0 else 0 for m, s, c in zip(minfi_dmps, sesame_dmps, consensus_dmps)]

    bars = ax.bar(x, [c/1000 for c in consensus_dmps], width * 1.5, label='Consensus DMPs (K)',
                  color='#45B7D1', edgecolor='black', linewidth=0.5)

    # Add percentage text on bars
    for i, (bar, spec) in enumerate(zip(bars, specificity)):
        height = bar.get_height()
        if consensus_dmps[i] > 0:
            ax.annotate(f'{spec:.0f}%\noverlap',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=7)

    ax.set_ylabel('Consensus DMPs (thousands)')
    ax.set_title('C. High-Confidence Consensus', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n(n={n})" for d, n in zip(datasets, n_samples)], fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Figure2_Pipeline_Comparison.png'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(out_dir, 'Figure2_Pipeline_Comparison.pdf'), facecolor='white')
    plt.close()
    print("  Created Figure 2: Pipeline comparison")


def create_figure3_ablation(ablation_data, out_dir):
    """
    Figure 3: Ablation study - effect of batch correction on lambda
    Shows raw vs corrected lambda for each pipeline/dataset
    """
    if not ablation_data:
        print("  Skipping Figure 3: No ablation data")
        return

    # Organize data by gse_id and pipeline
    data_by_gse = {}
    for row in ablation_data:
        if row.get('metric') != 'lambda':
            continue
        gse = row.get('gse_id', 'Unknown')
        pipeline = row.get('pipeline', '')
        variant = row.get('variant', '')
        value = safe_float(row.get('value', 1.0), 1.0)

        if gse not in data_by_gse:
            data_by_gse[gse] = {}
        key = f"{pipeline}_{variant}"
        data_by_gse[gse][key] = value

    if not data_by_gse:
        print("  Skipping Figure 3: No lambda data found")
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))

    datasets = list(data_by_gse.keys())
    x = np.arange(len(datasets))
    width = 0.2

    # Get values for each category
    minfi_raw = [data_by_gse[d].get('Minfi_baseline', 1.0) for d in datasets]
    minfi_corr = [data_by_gse[d].get('Minfi_corrected', 1.0) for d in datasets]
    sesame_raw = [data_by_gse[d].get('Sesame_baseline', 1.0) for d in datasets]
    sesame_corr = [data_by_gse[d].get('Sesame_corrected', 1.0) for d in datasets]

    ax.bar(x - 1.5*width, minfi_raw, width, label='Minfi (raw)', color='#4ECDC4', alpha=0.5, edgecolor='black', linewidth=0.5)
    ax.bar(x - 0.5*width, minfi_corr, width, label='Minfi (corrected)', color='#4ECDC4', edgecolor='black', linewidth=0.5)
    ax.bar(x + 0.5*width, sesame_raw, width, label='Sesame (raw)', color='#FF6B6B', alpha=0.5, edgecolor='black', linewidth=0.5)
    ax.bar(x + 1.5*width, sesame_corr, width, label='Sesame (corrected)', color='#FF6B6B', edgecolor='black', linewidth=0.5)

    ax.axhline(y=1.0, color='green', linestyle='--', linewidth=1.5, label='Ideal (λ=1.0)')
    ax.axhspan(0.9, 1.1, alpha=0.1, color='green')

    ax.set_ylabel('Genomic Inflation (λ)')
    ax.set_xlabel('Dataset')
    ax.set_title('Effect of Batch Correction on Genomic Inflation', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=9)
    ax.legend(loc='upper right', fontsize=8, ncol=2)
    ax.set_ylim(0.5, max(max(minfi_raw), max(sesame_raw), 1.5) * 1.1)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Figure3_Ablation_Lambda.png'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(out_dir, 'Figure3_Ablation_Lambda.pdf'), facecolor='white')
    plt.close()
    print("  Created Figure 3: Ablation study")


def create_figure4_crf_tiers(out_dir):
    """
    Figure 4: CRF v2.1 Sample Size Tiers
    """
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))

    tiers = ['Minimal\n(n<12)', 'Small\n(12-23)', 'Moderate\n(24-49)', 'Large\n(n≥50)']
    features = ['SVA', 'ComBat', 'Full MMC', 'Full NCS', 'Full SSS']

    # Feature availability matrix (1=available, 0=not available, 0.5=limited)
    availability = np.array([
        [0, 0, 0.5, 0.5, 0.5],   # Minimal
        [0, 0.5, 0.5, 0.5, 1],   # Small
        [1, 1, 1, 1, 1],          # Moderate
        [1, 1, 1, 1, 1],          # Large
    ])

    colors = ['#ff6b6b', '#ffd93d', '#6bcb77']
    cmap = plt.cm.colors.ListedColormap(colors)

    im = ax.imshow(availability.T, cmap=cmap, aspect='auto', vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(tiers)))
    ax.set_yticks(np.arange(len(features)))
    ax.set_xticklabels(tiers, fontsize=10)
    ax.set_yticklabels(features, fontsize=10)

    # Add text annotations
    for i in range(len(features)):
        for j in range(len(tiers)):
            val = availability[j, i]
            text = '✓' if val == 1 else ('○' if val == 0.5 else '✗')
            color = 'white' if val == 0 else 'black'
            ax.text(j, i, text, ha='center', va='center', fontsize=14, color=color, fontweight='bold')

    ax.set_title('CRF v2.1: Sample-Size Adaptive Feature Availability', fontweight='bold', pad=10)
    ax.set_xlabel('Sample Size Tier', fontweight='bold')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#6bcb77', label='Available'),
        Patch(facecolor='#ffd93d', label='Limited'),
        Patch(facecolor='#ff6b6b', label='Disabled'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1))

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Figure4_CRF_Tiers.png'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(out_dir, 'Figure4_CRF_Tiers.pdf'), facecolor='white')
    plt.close()
    print("  Created Figure 4: CRF tiers")


def create_figure5_cai_metrics(summary_data, results_dir, out_dir):
    """
    Figure 5: Correction Adequacy Index (CAI) breakdown
    Shows calibration, preservation, and batch removal scores
    """
    results_path = Path(results_dir)
    cai_data = []

    for gse_dir in sorted(results_path.glob('GSE*')):
        cai_file = gse_dir / 'illumeta_results' / 'Correction_Adequacy_Summary.csv'
        if not cai_file.exists():
            continue

        gse_id = gse_dir.name
        metrics = {}
        with open(cai_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                metrics[row.get('metric', '')] = safe_float(row.get('value', 0))

        if metrics:
            cai_data.append({
                'gse_id': gse_id,
                'calibration': metrics.get('calibration_score', 0),
                'preservation': metrics.get('preservation_score', 0),
                'batch_removal': metrics.get('batch_removal_score', 0),
                'cai': metrics.get('cai', 0),
            })

    if not cai_data:
        print("  Skipping Figure 5: No CAI data")
        return

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    # Panel A: Component scores
    ax = axes[0]
    datasets = [d['gse_id'] for d in cai_data]
    x = np.arange(len(datasets))
    width = 0.25

    calibration = [d['calibration'] for d in cai_data]
    preservation = [d['preservation'] for d in cai_data]
    batch = [d['batch_removal'] for d in cai_data]

    ax.bar(x - width, calibration, width, label='Calibration', color='#FF6B6B', edgecolor='black', linewidth=0.5)
    ax.bar(x, preservation, width, label='Preservation', color='#4ECDC4', edgecolor='black', linewidth=0.5)
    ax.bar(x + width, batch, width, label='Batch Removal', color='#45B7D1', edgecolor='black', linewidth=0.5)

    ax.set_ylabel('Score (0-1)')
    ax.set_title('A. CAI Component Scores', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=9)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_ylim(0, 1.1)

    # Panel B: Overall CAI
    ax = axes[1]
    cai_values = [d['cai'] for d in cai_data]
    colors = ['#6bcb77' if v >= 0.7 else '#ffd93d' if v >= 0.5 else '#ff6b6b' for v in cai_values]

    bars = ax.bar(datasets, cai_values, color=colors, edgecolor='black', linewidth=0.5)

    # Add threshold lines
    ax.axhline(y=0.7, color='green', linestyle='--', linewidth=1, alpha=0.7, label='Good (≥0.7)')
    ax.axhline(y=0.5, color='orange', linestyle='--', linewidth=1, alpha=0.7, label='Acceptable (≥0.5)')

    # Add value labels on bars
    for bar, val in zip(bars, cai_values):
        ax.annotate(f'{val:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel('Correction Adequacy Index')
    ax.set_title('B. Overall CAI Score', fontweight='bold')
    ax.set_ylim(0, 1.1)
    ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Figure5_CAI_Metrics.png'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(out_dir, 'Figure5_CAI_Metrics.pdf'), facecolor='white')
    plt.close()
    print("  Created Figure 5: CAI metrics")


def create_table1_features(out_dir):
    """
    Table 1: Feature comparison between tools
    """
    headers = ['Feature', 'minfi', 'ChAMP', 'IlluMeta']
    rows = [
        ['Normalization', 'Noob/Funnorm', 'BMIQ', 'Noob + pOOBAH'],
        ['Dual Pipeline', '✗', '✗', '✓'],
        ['Consensus Calling', '✗', '✗', '✓'],
        ['Automated Batch Correction', '✗', 'ComBat only', 'SVA/ComBat/limma (auto-select)'],
        ['Sample Size Adaptation', '✗', '✗', 'CRF v2.1'],
        ['Interactive Dashboard', '✗', '✗', '✓'],
        ['GEO Integration', 'Manual', 'Manual', 'Built-in download'],
        ['Auto-generated Methods', '✗', '✗', '✓'],
        ['Cell Composition', 'FlowSorted', 'RefbaseEWAS', 'EpiDISH + RefFreeEWAS'],
        ['Epigenetic Clocks', '✗', '✗', 'methylclock/planet'],
    ]

    # Write as TSV
    with open(os.path.join(out_dir, 'Table1_Feature_Comparison.tsv'), 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')

    # Write as formatted markdown
    with open(os.path.join(out_dir, 'Table1_Feature_Comparison.md'), 'w') as f:
        f.write('# Table 1: Feature Comparison\n\n')
        f.write('| ' + ' | '.join(headers) + ' |\n')
        f.write('|' + '|'.join(['---'] * len(headers)) + '|\n')
        for row in rows:
            f.write('| ' + ' | '.join(row) + ' |\n')

    print("  Created Table 1: Feature comparison")


def create_table2_benchmark_summary(summary_data, out_dir):
    """
    Table 2: Benchmark summary across datasets
    """
    if not summary_data:
        print("  Skipping Table 2: No summary data")
        return

    headers = ['Dataset', 'n', 'CRF Tier', 'Comparison', 'Minfi DMPs', 'Sesame DMPs', 'Consensus DMPs', 'λ (Minfi)', 'λ (Sesame)']

    rows = []
    for row in summary_data:
        comparison = f"{row.get('group_test', 'Test')[:15]} vs {row.get('group_con', 'Control')[:15]}"
        rows.append([
            row.get('gse_id', ''),
            row.get('n_samples', ''),
            row.get('crf_tier', '').capitalize(),
            comparison,
            f"{safe_int(row.get('minfi_dmps_fdr', 0)):,}",
            f"{safe_int(row.get('sesame_dmps_fdr', 0)):,}",
            f"{safe_int(row.get('illumeta_intersect_dmps', 0)):,}",
            f"{safe_float(row.get('lambda_minfi', 1.0)):.3f}",
            f"{safe_float(row.get('lambda_sesame', 1.0)):.3f}",
        ])

    # Write as TSV
    with open(os.path.join(out_dir, 'Table2_Benchmark_Summary.tsv'), 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')

    # Write as markdown
    with open(os.path.join(out_dir, 'Table2_Benchmark_Summary.md'), 'w') as f:
        f.write('# Table 2: Benchmark Summary\n\n')
        f.write('| ' + ' | '.join(headers) + ' |\n')
        f.write('|' + '|'.join(['---'] * len(headers)) + '|\n')
        for row in rows:
            f.write('| ' + ' | '.join(str(x) for x in row) + ' |\n')

    print("  Created Table 2: Benchmark summary")


def create_supplementary_table_ablation(ablation_data, out_dir):
    """
    Supplementary Table: Ablation study results
    """
    if not ablation_data:
        print("  Skipping Supplementary Table: No ablation data")
        return

    # Write full ablation data
    with open(os.path.join(out_dir, 'TableS1_Ablation_Full.tsv'), 'w') as f:
        if ablation_data:
            headers = list(ablation_data[0].keys())
            f.write('\t'.join(headers) + '\n')
            for row in ablation_data:
                f.write('\t'.join(str(row.get(h, '')) for h in headers) + '\n')

    print("  Created Supplementary Table S1: Full ablation results")


def collect_illumeta_results(results_dir):
    """
    Collect metrics from IlluMeta results directories.
    Returns list of dicts with benchmark summary data.
    """
    summary = []
    results_path = Path(results_dir)

    for gse_dir in sorted(results_path.glob('GSE*')):
        illumeta_results = gse_dir / 'illumeta_results'
        if not illumeta_results.exists():
            continue

        gse_id = gse_dir.name
        row = {'gse_id': gse_id}

        # Get sample count from CRF tier
        tier_file = illumeta_results / 'CRF_Sample_Tier.csv'
        if tier_file.exists():
            with open(tier_file, 'r') as f:
                reader = csv.DictReader(f)
                for r in reader:
                    row['n_samples'] = r.get('total_n', '0')
                    row['crf_tier'] = r.get('tier', '')
                    break

        # Get group info from Input_Group_Distribution
        group_file = illumeta_results / 'Input_Group_Distribution.csv'
        if group_file.exists():
            with open(group_file, 'r') as f:
                reader = csv.DictReader(f)
                groups = list(reader)
                if len(groups) >= 2:
                    row['group_con'] = groups[0].get('group', groups[0].get('primary_group', 'Control'))[:25]
                    row['group_test'] = groups[1].get('group', groups[1].get('primary_group', 'Test'))[:25]

        # Get lambda from ablation summary (corrected value)
        for pipeline in ['Minfi', 'Sesame']:
            abl_file = illumeta_results / f'{pipeline}_Ablation_Summary.csv'
            if abl_file.exists():
                with open(abl_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        if r.get('metric') == 'lambda':
                            row[f'lambda_{pipeline.lower()}'] = r.get('corrected', r.get('raw', '1.0'))
                            break

        # Count DMPs
        for pipeline in ['Minfi', 'Sesame']:
            dmp_file = illumeta_results / f'{pipeline}_DMPs_full.csv'
            if dmp_file.exists():
                count = 0
                with open(dmp_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        try:
                            fdr = float(r.get('adj.P.Val', r.get('FDR', '1')))
                            if fdr < 0.05:
                                count += 1
                        except (ValueError, TypeError):
                            pass
                row[f'{pipeline.lower()}_dmps_fdr'] = count

        # Consensus DMPs
        consensus_file = illumeta_results / 'Intersection_Consensus_DMPs.csv'
        if consensus_file.exists():
            with open(consensus_file, 'r') as f:
                row['illumeta_intersect_dmps'] = sum(1 for _ in f) - 1  # subtract header

        if row.get('n_samples'):
            summary.append(row)
            print(f"  Collected metrics from {gse_id} (n={row.get('n_samples')}, tier={row.get('crf_tier', 'N/A')})")

    return summary


def collect_ablation_long(results_dir):
    """
    Collect ablation data in long format for visualization.
    """
    ablation = []
    results_path = Path(results_dir)

    for gse_dir in sorted(results_path.glob('GSE*')):
        illumeta_results = gse_dir / 'illumeta_results'
        gse_id = gse_dir.name

        for pipeline in ['Minfi', 'Sesame']:
            abl_file = illumeta_results / f'{pipeline}_Ablation_Summary.csv'
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
    parser.add_argument('--ablation-tsv', default=None,
                        help='Path to ablation metrics TSV (optional)')
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
    if args.ablation_tsv and os.path.exists(args.ablation_tsv):
        ablation_data = load_ablation_metrics(args.ablation_tsv)
        print(f"Loaded {len(ablation_data)} ablation metrics from TSV")
    elif os.path.exists(args.results_dir):
        ablation_data = collect_ablation_long(args.results_dir)
        print(f"Collected {len(ablation_data)} ablation metrics from results")

    print("\nGenerating figures...")

    # Generate figures
    create_figure1_workflow(args.out_dir)
    create_figure2_benchmark_comparison(summary_data, args.out_dir)
    create_figure3_ablation(ablation_data, args.out_dir)
    create_figure4_crf_tiers(args.out_dir)
    create_figure5_cai_metrics(summary_data, args.results_dir, args.out_dir)

    print("\nGenerating tables...")

    # Generate tables
    create_table1_features(args.out_dir)
    create_table2_benchmark_summary(summary_data, args.out_dir)
    create_supplementary_table_ablation(ablation_data, args.out_dir)

    print(f"\n✅ Done! Files saved to: {args.out_dir}/")
    print("\nGenerated files:")
    for f in sorted(os.listdir(args.out_dir)):
        size = os.path.getsize(os.path.join(args.out_dir, f))
        print(f"  - {f} ({size:,} bytes)")


if __name__ == '__main__':
    main()
