# GSE66351 stratified AD-vs-CTRL run commands

Generated for Task 5. Outputs are non-overwriting relative to prior bulk runs.

## Feasible strata

- Occipital Neuron: AD=15, CTRL=16 (`configure_AD_vs_CTRL_neuron.tsv`)
- Occipital Glia: AD=15, CTRL=16 (`configure_AD_vs_CTRL_glia.tsv`)
- Bulk frontal cortex: AD=37, CTRL=26 (documented but not rerun here)
- Bulk temporal cortex: AD=39, CTRL=26 (documented but not rerun here)

## Highest-value run started

```bash
ILLUMETA_RESPECT_R_LIBS_USER=1 \
R_LIBS_USER=/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.5 \
LC_ALL=C LANG=C python illumeta.py analysis \
  -c projects/GSE66351/configure_AD_vs_CTRL_neuron.tsv \
  --idat-dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  -o projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_lcC \
  --group_con CTRL --group_test AD \
  --include-covariates age,sex,Sentrix_ID,Sentrix_Position \
  --batch-method limma --vp-top 0 --permutations 0 \
  --tissue Auto --min-total-size 20
```

## Feasible follow-up command, not started until neuron result is available

```bash
ILLUMETA_RESPECT_R_LIBS_USER=1 \
R_LIBS_USER=/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.5 \
LC_ALL=C LANG=C python illumeta.py analysis \
  -c projects/GSE66351/configure_AD_vs_CTRL_glia.tsv \
  --idat-dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  -o projects/GSE66351/AD_vs_CTRL_glia_occipital_results_lcC \
  --group_con CTRL --group_test AD \
  --include-covariates age,sex,Sentrix_ID,Sentrix_Position \
  --batch-method limma --vp-top 0 --permutations 0 \
  --tissue Auto --min-total-size 20
```

## Rerun after Tier3 guard stop

The first direct run reached Minfi confounding diagnostics and stopped at Tier3 safety because a small sorted-cell stratum had `min_stratum_n=1`. The completed/continuing rerun therefore keeps this as a guardrail and proceeds with explicit skip:

```bash
R_LIBS=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5:/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.5 \
R_LIBS_USER=/home/keunsoo/R/x86_64-pc-linux-gnu-library/4.5 \
LC_ALL=C LANG=C Rscript projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_lcC/_illumeta_code_snapshot/analyze.R \
  --config projects/GSE66351/configure_AD_vs_CTRL_neuron.tsv \
  --out projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_direct_run5_tier3skip \
  --group_con CTRL --group_test AD \
  --max_plots 10000 --pval 0.05 --lfc 0.5 --delta-beta 0.0 \
  --permutations 0 --min-total-size 20 --qc-intensity-threshold 9.0 \
  --batch-method limma --tier3-on-fail skip \
  --idat-dir /home/keunsoo/Projects/23_illumeta/projects/GSE66351/idat \
  --include-covariates age,sex,Sentrix_ID,Sentrix_Position \
  --tissue Auto --vp-top 0
```
