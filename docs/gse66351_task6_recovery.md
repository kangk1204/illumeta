# Task 6 recovery: GSE66351 occipital neuron stratified evidence

Task 6 recovers useful Lane 1b evidence after Task 5 was marked failed during the first full-run attempt. The recovery path first produced a non-overwriting Minfi-only run for the highest-value cell-type stratum and records verified configs for both feasible occipital cell-type strata.

A later worker-1 direct run also reached terminal full dual-pipeline outputs for the same occipital neuron stratum. The final interpretation therefore has two layers: the Task 6 Minfi-only recovery proves the stratum is analyzable, and the late full direct run is the preferred terminal IlluMeta evidence for manuscript-facing statements.

## Config evidence

Generated configs use the explicit schema from `configure_AD_vs_CTRL_bulk.tsv`:

`primary_group`, `geo_accession`, `cell_type`, `diagnosis`, `braak_stage`, `brain_region`, `age`, `sex`, `donor_id`, `Sentrix_ID`, `Sentrix_Position`, `source_name_ch1`, `description`.

| File | Stratum | AD | CTRL | Verification |
| --- | --- | ---: | ---: | --- |
| `projects/GSE66351/configure_AD_vs_CTRL_neuron_occipital.tsv` | Occipital cortex / Neuron | 15 | 16 | 31 rows; all `cell_type=Neuron`; all `brain_region=Occipital cortex`. |
| `projects/GSE66351/configure_AD_vs_CTRL_glia_occipital.tsv` | Occipital cortex / Glia | 15 | 16 | 31 rows; all `cell_type=Glia`; all `brain_region=Occipital cortex`. |
| `projects/GSE66351/gse66351_stratified_config_summary.tsv` | Summary | - | - | Also records frontal bulk 37/26 and temporal bulk 39/26 as documented-not-run strata. |

## Recovery run

Terminal recovery output:

`projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_worker2_minfi_vp0_limmavp0`

Command pattern used:

```bash
R_LIBS=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5:/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.5 \
R_LIBS_USER=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5 \
LC_ALL=C LANG=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Rscript r_scripts/analyze.R \
  --config projects/GSE66351/configure_AD_vs_CTRL_neuron_occipital.tsv \
  --out projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_worker2_minfi_vp0_limmavp0 \
  --group_con CTRL --group_test AD \
  --idat_dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  --cross_reactive_list /home/keunsoo/Projects/23_illumeta/projects/GSE66351/references/probe_blacklists/cross_reactive_450K.tsv \
  --skip-sesame \
  --max_plots 10000 --pval 0.05 --lfc 0.5 --delta_beta 0 \
  --permutations 0 --min_total_size 20 --qc_intensity_threshold 9.0 \
  --batch_method limma --include_covariates age,sex,Sentrix_ID,Sentrix_Position \
  --tissue DLPFC --preset conservative --vp_top 0 --tier3-on-fail skip
```

The run completed at `2026-05-21 23:08:44` and produced `methods.md`, `decision_ledger.tsv`, `QC_Summary.csv`, `Preflight_Summary.csv`, `Minfi_Metrics.csv`, `Minfi_DMPs_full.csv`, `Minfi_DMRs.csv`, plots, and session information.

## Late full direct run

Terminal full-output directory:

`projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_direct_run5_tier3skip`

This direct run completed after Task 5 had already been transitioned to failed in the team audit state. It produced `summary.json`, `analysis_parameters.json`, Minfi/Sesame/Sesame_Native metrics, DMR tables, strict/native intersection files, methods, decision ledger, QC/preflight summaries, and log evidence. The main working tree preserves the compact terminal evidence subset; the worker worktree retained the much larger matrix, plot, and full DMP artifacts.

| Output metric | Value | Source |
| --- | ---: | --- |
| Control / AD samples | 16 / 15 | `summary.json` |
| Strict consensus DMP rows | 0 | `Intersection_Consensus_DMPs.csv` |
| Native consensus DMP rows | 0 | `Intersection_Native_Consensus_DMPs.csv` |
| Minfi DMR rows | 651 | `Minfi_DMRs.csv` |
| Sesame DMR rows | 828 | `Sesame_DMRs.csv` |
| Sesame native DMR rows | 807 | `Sesame_Native_DMRs.csv` |
| Primary result mode | tier3_ineligible | `summary.json` |
| Primary no-signal flag | true | `summary.json` |

## Verified terminal evidence

| Output metric | Value | Source |
| --- | ---: | --- |
| Config samples | 31 | `Preflight_Summary.csv` |
| Control / AD samples | 16 / 15 | `Input_Group_Distribution.csv`, `Preflight_Summary.csv` |
| IDAT pairs missing | 0 | `Preflight_Summary.csv` |
| Samples passed QC | 31 | `QC_Summary.csv` |
| Detection-P probes removed | 2,374 | `Probe_Filter_Summary.csv` |
| Cross-reactive probes removed | 38,638 | `Probe_Filter_Summary.csv`, `decision_ledger.tsv` |
| Final probes | 420,309 | `Probe_Filter_Summary.csv` |
| Minfi genomic inflation lambda | 0.926590458015778 | `Minfi_Metrics.csv` |
| Lambda guard status | ok | `Minfi_Metrics.csv` |
| Primary result mode | tier3_ineligible | `Minfi_Metrics.csv` |
| DMR rows | 865 | `Minfi_DMRs.csv` |

## Caveats for manuscript integration

- Task 6 itself remains **Minfi-only recovery evidence**. Sesame was intentionally skipped in that recovery run to avoid the full-run runtime blocker.
- The late worker-1 full direct run is now the preferred terminal dual-pipeline evidence for the occipital neuron stratum, but it remains a small stratum (AD=15, CTRL=16) with `tier3_ineligible` primary mode and zero strict/native consensus DMPs.
- The result should be presented as guarded cell-type/region sensitivity evidence, not as a positive AD biomarker discovery result.
- Original Task 5 remains an acknowledged failed attempt in the OMX task audit because it was transitioned before the late terminal run finished. Task 6 plus this late-run reconciliation supplies the corrected evidence trail.
