# Cell-type-aware neuroepigenomics reframe for the dementia IlluMeta package

Generated: 2026-05-22T15:37:23.496342+00:00

## Source-grounded conclusion

The current package is strongest as a Methods/YMETH-compatible neuroepigenomics workflow paper, not as a universal Alzheimer disease biomarker paper. The completed AD/dementia public-IDAT analyses cover 771 controls and 883 AD/test samples across blood, bulk brain, EPIC brain, BS-only brain, and frontal/temporal cortex contexts. They produced 2562 strict and 3386 native consensus DMPs in aggregate, but the evidence is heterogeneous: 3/6 cohorts triggered primary lambda guard; six AD/dementia cohorts are workflow evidence with heterogeneity: 3 strict-zero and 2 native-zero consensus cohorts, plus only 4 nonzero pairwise CpG-overlap rows.

## Guarded AD evidence

- `GSE125895` remains the cleanest brain AD anchor because it has nonzero strict/native consensus and primary lambda guard `ok`.
- `GSE284764` provides recent EPIC prefrontal-cortex support, but its primary lambda guard triggered, so it remains guarded support.
- `GSE134379` and `GSE105109` completed branch-level DMP/DMR outputs but produced no consensus DMPs.
- `GSE66351` should be described as a bulk-cortex sensitivity lane: 0 strict and 10 native consensus DMPs; Bulk frontal/temporal cortex AD run supports heterogeneity/guardrail narrative rather than a universal AD CpG panel.

## Neural cell-type reference axis

`GSE306226/Neurons_vs_Microglia_results` provides the local positive-control/reference axis for cell-type-aware framing [24,25]. It compares 20 neurons with 20 microglia, passed sample QC, and produced 6398 strict plus 6269 native Minfi/SeSAMe consensus DMPs. Branch concordance was strict logFC r=0.791, Jaccard=0.302; native logFC r=0.761, Jaccard=0.304. Because primary=triggered; Minfi=triggered; SeSAMe=triggered; native=triggered, it should be used as a reference axis demonstrating that IlluMeta can recover a large neural cell-type contrast, not as an unqualified discovery claim.


## Guarded GSE306227 extension

`GSE306227/Neurons_vs_Microglia_results` is retained as a guarded independent sorted-brain extension, not as a core Figure 3 axis in the current package [24,26]. It has summary-level evidence for 324706 strict and 334911 native consensus DMPs across 18 microglia/control-axis and 19 neuron/test-axis samples, with primary result mode `tier3_ineligible` and lambda guard `triggered`. However, branch metrics present=False, cell summary present=False, and the run log records RefFreeEWAS unavailable=True. Therefore it should be cited only as guarded extension evidence until branch-level design and cell-adjustment artifacts are synced.


## Manuscript claim guardrail

Defensible central claim: IlluMeta exposes context dependence in public neurodegeneration methylation data by making raw-IDAT import, branch-specific covariate/cell/SV adjustment, dual-pipeline agreement, consensus DMP/DMR outputs, and lambda/QC guards auditable. The manuscript should explicitly avoid phrases such as "universal AD biomarker" or "validated AD signature" unless tied to a specific cohort/context and guard status.
