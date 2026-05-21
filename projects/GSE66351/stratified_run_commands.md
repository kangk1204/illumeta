# GSE66351 stratified AD-vs-CTRL run commands

Generated for Task 5. Outputs are non-overwriting relative to prior bulk runs.

## Feasible strata

- Occipital Neuron: AD=15, CTRL=16 (`configure_AD_vs_CTRL_neuron.tsv`)
- Occipital Glia: AD=15, CTRL=16 (`configure_AD_vs_CTRL_glia.tsv`)
- Bulk frontal cortex: AD=37, CTRL=26 (documented but not rerun here)
- Bulk temporal cortex: AD=39, CTRL=26 (documented but not rerun here)

## Completed highest-value run

The initial Python wrapper run reached final modeling but stopped at Tier3 safety because the small sorted-cell stratum had `min_stratum_n=1`. The completed run therefore preserved the guardrail and proceeded with explicit Tier3 skip:

```bash
R_LIBS=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5:/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.5 \
R_LIBS_USER=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5 \
LC_ALL=C LANG=C Rscript projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_lcC/_illumeta_code_snapshot/analyze.R \
  --config projects/GSE66351/configure_AD_vs_CTRL_neuron.tsv \
  --out projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_direct_run5_tier3skip \
  --group_con CTRL --group_test AD \
  --max_plots 10000 --pval 0.05 --lfc 0.5 --delta_beta 0.0 \
  --permutations 0 --min_total_size 20 --qc_intensity_threshold 9.0 \
  --batch_method limma --tier3-on-fail skip \
  --idat_dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  --include_covariates age,sex,Sentrix_ID,Sentrix_Position \
  --tissue Auto --vp_top 0
```

## Optional follow-up not run

```bash
ILLUMETA_RESPECT_R_LIBS_USER=1 \
R_LIBS_USER=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5 \
LC_ALL=C LANG=C python illumeta.py analysis \
  -c projects/GSE66351/configure_AD_vs_CTRL_glia.tsv \
  --idat_dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  -o projects/GSE66351/AD_vs_CTRL_glia_occipital_results_lcC \
  --group_con CTRL --group_test AD \
  --include_covariates age,sex,Sentrix_ID,Sentrix_Position \
  --batch_method limma --tier3-on-fail skip --vp_top 0 --permutations 0 \
  --tissue Auto --min_total_size 20
```
