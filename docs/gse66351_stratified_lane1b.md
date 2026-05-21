# GSE66351 stratified AD-vs-control Lane 1b handoff

Task 5 creates non-overwriting, analysis-ready configs for feasible GSE66351 strata and records the run strategy used for the highest-value stratum. The intent is to treat GSE66351 as cell-type/region-sensitive AD follow-up evidence, not as a pooled universal AD biomarker cohort.

## Generated config files

| Config | Stratum | AD | CTRL | Notes |
| --- | --- | ---: | ---: | --- |
| `projects/GSE66351/configure_AD_vs_CTRL_neuron_occipital.tsv` | Occipital cortex / Neuron | 15 | 16 | Highest-value cell-type stratum; selected for worker-2 run. |
| `projects/GSE66351/configure_AD_vs_CTRL_glia_occipital.tsv` | Occipital cortex / Glia | 15 | 16 | Feasible paired cell-type stratum; config verified, not run by worker-2 unless time permits. |
| `projects/GSE66351/gse66351_stratified_config_summary.tsv` | Counts summary | - | - | Also documents frontal and temporal bulk strata. |

The generated configs use the explicit bulk-config schema:
`primary_group`, `geo_accession`, `cell_type`, `diagnosis`, `braak_stage`, `brain_region`, `age`, `sex`, `donor_id`, `Sentrix_ID`, `Sentrix_Position`, `source_name_ch1`, `description`.

## Bulk region counts documented but not run by worker-2

| Stratum | AD | CTRL | Reason |
| --- | ---: | ---: | --- |
| Frontal cortex / bulk | 37 | 26 | Already represented by prior bulk GSE66351 evidence; task priority was cell-type strata. |
| Temporal cortex / bulk | 39 | 26 | Already represented by prior bulk GSE66351 evidence; task priority was cell-type strata. |

## Run strategy

Worker-2 launched a non-overwriting neuron stratum run under:

`projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_worker2_vp0_limmavp0_run2`

Command shape:

```bash
R_LIBS=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5:/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.5 \
R_LIBS_USER=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5 \
LC_ALL=C LANG=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Rscript r_scripts/analyze.R \
  --config projects/GSE66351/configure_AD_vs_CTRL_neuron_occipital.tsv \
  --out projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_worker2_vp0_limmavp0_run2 \
  --group_con CTRL --group_test AD \
  --idat_dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  --cross_reactive_list /home/keunsoo/Projects/23_illumeta/projects/GSE66351/references/probe_blacklists/cross_reactive_450K.tsv \
  --max_plots 10000 --pval 0.05 --lfc 0.5 --delta_beta 0 \
  --permutations 0 --min_total_size 20 --qc_intensity_threshold 9.0 \
  --batch_method limma --include_covariates age,sex,Sentrix_ID,Sentrix_Position \
  --tissue DLPFC --preset conservative --vp_top 0 --tier3-on-fail skip
```

Guardrails:

- `--vp_top 0` follows the safe run5 pattern for avoiding variance-partition instability on this lane.
- `--batch_method limma` preserves the explicit/manual batch-handling pattern from the latest bulk run.
- The 450K cross-reactive blacklist is supplied explicitly from the existing GSE66351 reference bundle; the first direct run without this path stopped at the mandatory blacklist guard.
- Outputs are written under new worker-2 result directories only; existing bulk and Lane 2/3 outputs are not overwritten.
