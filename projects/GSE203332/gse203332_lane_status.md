# GSE203332 frontal cortex AD specificity/stress-test lane status

## Raw/metadata gate
- GEO matrix metadata rows: 492 on GPL21145 EPIC.
- Raw IDAT links in metadata: 984 = 492 sample pairs.
- NCBI filelist gate (verified live in worker log): `GSE203332_RAW.tar` plus 984 IDAT files.
- QC exclusion before primary contrast: keep `is_technical_replicate == no` and `failed_genotyping_qc == no` (442 samples).

## Locked primary contrast
- Primary config: `projects/GSE203332/configure_AD_vs_Control_neuropath.tsv`.
- Contrast: neuropathological `AD` vs `CONTR`, encoded as `AD` vs `Control`.
- Counts after QC gate: AD = 61; Control = 74.
- Rationale: neuropathological AD vs CONTR is cleaner than clinical diagnosis, which includes discordant DLB/PDD/PD/mixed pathology rows.
- Planned covariates/batch fields: age, sex, pmd_min, amp_plate, slide; exclude Braak/CERAD/Thal/asyn stages from adjustment because they are pathology-linked interpretation fields.

## Secondary specificity design
- Optional secondary contrast table: `projects/GSE203332/gse203332_secondary_specificity_design.tsv`.
- AD vs DLB/PDD/PD/mixed_AD_LBD/iLBD should be interpreted only as disease-specificity/stress-test evidence, not universal AD biomarker discovery.

## Manuscript-use guard
- Use as frontal-cortex mixed neurodegeneration stress-test if the full IlluMeta run completes with branch/lambda evidence.
- If blocked, cite the raw/metadata gate and locked contrast as readiness evidence only; do not count it as a completed cohort.

## Runtime attempt 1
- Attempted full IlluMeta run after selected-IDAT download: `projects/GSE203332/logs/AD_vs_Control_neuropath_results_lcC.log`.
- Status: blocked before methylation analysis during first-run R dependency installation, not by GSE203332 data design.
- Error signature: `stringi` failed to load because conda ICU (`libicui18n.so.78`) required `CXXABI_1.3.15` unavailable from system `/usr/lib/x86_64-linux-gnu/libstdc++.so.6`.
- Action: stopped the install/run process to avoid wasting runtime; retry should sanitize conda `PATH`/`PKG_CONFIG_PATH`/compiler library discovery or use the existing working `conda run -n illumeta` execution lane.

## Runtime attempt 2
- Retried via `conda run -n illumeta` to avoid the system-R/conda-ICU ABI mismatch.
- Status: stopped during first-run setup because default setup began downloading/installing optional cell-reference packages (`FlowSorted.*`) not required for the `--tissue Auto` GSE203332 stress-test lane.
- Action: retry with `ILLUMETA_INSTALL_MINIMAL=1` to keep the run to core IlluMeta dependencies and avoid optional cell-reference overhead.

## Runtime attempt 3
- Retried with `ILLUMETA_INSTALL_MINIMAL=1`.
- Status: blocked in `r_scripts/setup_env.R` before analysis start while downloading Bioconductor core dependency `GenomeInfoDb_1.42.3.tar.gz`; the target file stayed at 0 bytes for >30 seconds and the setup process remained in socket poll.
- No `summary.json`, consensus DMP tables, DMR outputs, or lambda guard were produced for GSE203332 in this worker attempt.
- Current defensible manuscript-use decision: readiness/blocked stress-test lane only. The raw gate, selected-IDAT download, and locked neuropathological AD-vs-CONTR design are complete; do not count GSE203332 as a completed IlluMeta cohort until the R dependency setup/download blocker is resolved and the full Minfi+SeSAMe run reaches terminal artifacts.
