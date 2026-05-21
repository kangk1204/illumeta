# Task 6 recovery: GSE66351 occipital neuron stratified evidence

Task 6 recovers useful Lane 1b evidence after the broader Task 5 attempt failed to produce a terminal full dual-pipeline IlluMeta result. The recovery path uses a non-overwriting Minfi-only run for the highest-value cell-type stratum and records verified configs for both feasible occipital cell-type strata.

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

- This is **recovery evidence**, not a completed full dual-pipeline IlluMeta consensus result. Sesame was intentionally skipped in the terminal recovery run to avoid the full-run runtime blocker.
- The first full worker-2 dual-pipeline attempt reached prefilter/QC and then entered Sesame, but did not reach terminal outputs before the recovery path was chosen. Its log is retained at `projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_worker2_vp0_limmavp0_run2.log` as partial evidence.
- Because Sesame is absent, strict/native intersection files and `summary.json` are not available for the recovery run; downstream text should report this as Minfi-only stratum evidence with explicit caveats.
- Original Task 5 remains an acknowledged failed attempt for terminal full IlluMeta output. Task 6 supplies verified configs plus a terminal Minfi-only neuron stratum recovery run.
