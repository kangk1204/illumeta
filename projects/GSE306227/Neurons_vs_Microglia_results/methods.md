# IlluMeta analysis methods (auto-generated)

- Generated: 2026-05-22 20:08:35
- Config: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/configure.tsv`
- Output: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results`
- Groups: control=`Microglia` (n=18), test=`Neurons` (n=19)
- Tissue: `Auto` (source: auto)
- Array type detected: EPIC
- Genomic coordinates: standardized to minfi annotation (hg19) for cross-pipeline comparison.
- Cell reference: default references (fallback to RefFreeEWAS when unavailable)
- Sesame cell composition: minfi-derived.
- Cell adjustment guard: action=warn when Eta^2 > 0.50.
- Significance thresholds: BH FDR < 0.050 and |log2FC| > 0.50

## Overview
IlluMeta performs an end-to-end DNA methylation analysis from raw Illumina IDAT files, running minfi and sesame independently and reporting Sesame in both strict (Minfi-aligned) and native (pOOBAH-preserving) views, plus consensus (intersection) call sets.

## Data input
- Raw IDATs are read from the project `idat/` directory.
- Samples are defined by `SampleID` and `primary_group` in `configure.tsv`.
- Auto-grouping (if used) is heuristic; users must verify group labels/counts before interpretation.

## Sample-level QC
- Detection P-value: samples with > 20% probes failing (P > 0.05) are excluded.
- Signal intensity QC: median methylated/unmethylated intensity < 9.0 (log2).
- Mixed-array safeguard: samples whose IDAT array size deviates from the modal array size are excluded (unless `--force_idat`).
- Sex mismatch QC: predicted sex via minfi::getSex() vs metadata column 'sex' (action=stop; mismatches=0).

## Probe-level QC (minfi-derived probe set)
- Detection P-value filter: probes are retained only if P < 0.05 in all retained samples.
- Cross-reactive probes: removed using list (n=44254). References: Pidsley 2016; McCartney 2016; Chen 2013.
- SNP filtering: probes overlapping common SNPs (SBE/CpG; MAF >= 0.01) are removed using `minfi::dropLociWithSnps()`.
- Sex chromosome probes (chrX/chrY) are removed.

## Normalization (two pipelines)
- **minfi**: `preprocessNoob()` followed by `mapToGenome()`; beta values are extracted with `getBeta()`.
- **sesame**: `qualityMask()` + `inferInfiniumIChannel()` + dye bias correction (default: `dyeBiasL()`; optional `dyeBiasCorrTypeINorm()` when `--sesame_typeinorm` is enabled with fallback to `dyeBiasL()`), followed by `pOOBAH()` masking and `noob()` background correction; beta values are extracted with `getBetas()`.
- **sesame strict view**: probes with any masked values are removed and then intersected with the minfi QC probe set for conservative cross-validation.
- **sesame native view**: probes with <= 0.10 missingness are retained; remaining NAs are imputed via KNN (impute::impute.knn, k=10) when available, with row-mean fallback if KNN is unavailable.

## Covariates and batch control
- Automatic covariate discovery: metadata variables associated with the top PCs (base alpha=0.010) are considered; per-variable significance is Bonferroni-corrected by the number of finite PCs tested for that variable (effective alpha = alpha / n_tested_PCs), then filtered for stability/confounding.
- Auto covariate guard: variables associated with the group (p < 0.001) are excluded by default.
- Auto-selected covariates may include mediators; verify biological plausibility before interpretation.
- Small-n safeguard: if covariate count would eliminate residual degrees of freedom, covariates are capped using PC-association/variance ranking to preserve model stability.
- Overcorrection guard: total covariates + SVs are capped at 25% of sample size (excess terms are dropped to preserve power).
- Undercorrection guard: if SVA detects hidden structure, at least 1 SV(s) are retained when possible (dropping non-forced covariates first).
- Cell composition: when tissue is `Auto`, IlluMeta attempts to infer tissue from metadata/heuristics; if unresolved, reference-free deconvolution is performed via RefFreeEWAS (K=5 latent components). For `Placenta`, IlluMeta uses planet::plCellCpGsThird with minfi::projectCellType (Houseman) when available. Custom references via `--cell_reference` override defaults when provided.
- Reference-based deconvolution can be biased if disease alters methylation states; interpret cell-adjusted models cautiously and consider RefFreeEWAS outputs.
- Sesame pipelines attempt Sesame-native cell composition; if unavailable, they try planet (Placenta), EpiDISH (Blood), or RefFreeEWAS on Sesame betas, and fall back to Minfi-derived covariates if needed.
- Surrogate variable analysis (SVA): enabled unless `--disable_sva`; SVs are estimated on top-variable probes and included in the model only when selected as the best batch strategy or when no batch factor is evaluated (to avoid double correction). SVs strongly associated with the group (P < 1.0e-06 or Eta^2 > 0.50) are excluded to avoid over-correction.
- Epigenetic clock covariates: when `--include_clock_covariates` is enabled, clock outputs are merged into metadata and considered by auto covariate selection (clocks with missing/constant values are excluded).
- Batch method comparison: batch candidates are screened for VOI confounding (Cramer's V / R2 and overlap tiers). Tier 3 triggers stratified/meta-analysis fallbacks when eligible; Tier 0-2 proceed to scoring.
- Tier3 eligibility constraints: min_total_n=20, min_per_group_per_stratum=5 (on_fail=skip).
- Optimization preset: `conservative` (weights on batch removal, biology preservation, calibration, stability).
- CRF sample tier: moderate (total_n=37; min_per_group=18).
- Calibration: permutation tests shuffle labels within batch strata to assess p-value uniformity (KS) and inflation.
- Lambda guard is a heuristic inflation check; interpret with EWAS correlation structure in mind.
- Correction Adequacy Framework (CAF): combines calibration (target FPR=0.05), signal preservation, and batch removal into a single CAI; see `Correction_Adequacy_Report.txt`.
- Decision ledger: automated decisions (covariates, batch choice, method) are logged with reasons.

## Differential methylation (DMP)
- Beta values are transformed to M-values with a logit offset of 0.000100.
- Differential methylation is tested using limma (`lmFit`/`eBayes(robust=TRUE)`) with a group contrast (test - control).
- Multiple testing is controlled by Benjamini-Hochberg FDR.

## Differentially methylated regions (DMR)
- DMRs are called with `dmrff` (maxgap=500; p.cutoff=0.050; min_cpgs=2).

## Consensus (intersection) call set
- Consensus DMPs are defined as CpGs significant in **both** minfi and sesame with the **same direction** under the same thresholds.
- Intersection is intended as a high-confidence subset; pipeline-specific results may capture additional true positives and are reported as sensitivity/discovery sets.
- Consensus table ranking uses Fisher's combined probability test (chi-squared, df=4) with genome-wide Benjamini-Hochberg FDR (`P.Value`/`adj.P.Val`). Selection-rule p-values are retained as `P.Value.selection` and `adj.P.Val.selection` (also mirrored in `*.max` legacy columns).
- Consensus is computed for both the strict (Minfi-aligned) and native Sesame views.
- Consensus outputs: `Intersection_Consensus_DMPs.*` and `Intersection_Native_Consensus_DMPs.*`, plus concordance/overlap plots.
- Primary branch is selected by the optimization score (see `results/consensus/primary_branch.txt`); the other branch is reported as sensitivity.
- Primary branch override: none.
- Primary inference mode: tier3_ineligible (tier3_batch=Sentrix_ID).
WARNING: Tier3 confounding detected but eligibility failed; stratified/meta-analysis was not run.
- Tier3 meta-analysis (primary): method=fixed, I2_median=NA.
- Lambda guard (primary): status=triggered, action=warn, threshold=1.500, lambda_guard_lambda=NA.
- DMR status (primary): ok.
- Sesame dye bias normalization: dyeBiasL (typeinorm_disabled).

## Reproducibility artifacts
- `analysis_parameters.json`: run parameters and thresholds.
- `sessionInfo.txt`: full R session/package versions.
- `code_version.txt`: git commit hash (when available).
- `config_used.yaml`: resolved config and preset details.
- `decision_ledger.tsv`: auditable record of automated decisions.

## QC summary (from QC_Summary.csv)
- Total_samples_input: 37
- Samples_failed_QC: 0
- Samples_passed_QC: 37
- Samples_failed_sex_mismatch: 0
- Sex_mismatch_samples: 0
- Total_probes_raw: 865859
- Probes_cross_reactive: 44254
- Probes_final: 772594
