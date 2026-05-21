# GSE66351 stratified lane status

Scope: Task 5 Lane 1b, occipital neuron/glia stratified AD-vs-CTRL analysis.

## Verified stratified inputs

- Occipital **Neuron**: AD = 15, CTRL = 16; config `projects/GSE66351/configure_AD_vs_CTRL_neuron.tsv`.
- Occipital **Glia**: AD = 15, CTRL = 16; config `projects/GSE66351/configure_AD_vs_CTRL_glia.tsv`.

## Bulk counts preserved for reference only

- Frontal cortex bulk: AD = 37, CTRL = 26.
- Temporal cortex bulk: AD = 39, CTRL = 26.

## Completed feasible run

- Completed highest-value stratum: occipital **Neuron** AD vs CTRL.
- Result directory: `projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_direct_run5_tier3skip`.
- Samples retained: CTRL = 16, AD = 15.
- Consensus DMPs: strict = 0, native = 0.
- Branch significant DMPs: Minfi = 0, SeSAMe strict = 0, SeSAMe native = 0.
- DMR rows: Minfi = 651, SeSAMe strict = 828, SeSAMe native = 807.
- Lambda guard: Minfi ok (0.920), SeSAMe strict ok (1.017), SeSAMe native ok (1.014).
- Result mode: `tier3_ineligible`; Tier3 confounding was detected for `Sentrix_ID`, but stratum eligibility failed (`min_stratum_n=1`), so the completed run used explicit `--tier3-on-fail skip`.

## Guarded interpretation

This is negative/guarded cell-type-stratified evidence: the neuron-only occipital AD-vs-CTRL lane completes with QC/lambda guards ok but zero DMP consensus. It should be used to support tissue/cell-type/context heterogeneity and boundary conditions, not universal AD biomarker claims.
