# GSE306226 Neurons-vs-Microglia positive-control audit

This audit treats `projects/GSE306226/Neurons_vs_Microglia_results` as a neural cell-type positive-control/reference axis for the dementia/AD reframing package. It is not an AD biomarker discovery result. The contrast is a purified neural cell-type comparison (`Neurons` vs `Microglia`) and is most useful for showing that IlluMeta can expose large, biologically coherent, cell-type-dependent methylation structure while preserving explicit QC and calibration caveats.

## Source-backed summary

| Evidence item | Value | Source file |
| --- | ---: | --- |
| Input groups | Microglia `n=20`; Neurons `n=20` | `Input_Group_Distribution.csv` |
| Samples passing QC | `40/40` | `QC_Summary.csv` |
| Raw probes | `865,859` | `QC_Summary.csv` |
| Final probes after detection / cross-reactive / SNP / sex-chromosome filters | `770,885` | `QC_Summary.csv` |
| Minfi significant DMP direction counts | up `4,513`; down `4,630` | `summary.json` |
| Sesame strict significant DMP direction counts | up `7,251`; down `12,413` | `summary.json` |
| Sesame native significant DMP direction counts | up `9,605`; down `8,908` | `summary.json` |
| Strict dual-pipeline consensus DMPs | `6,398` | `Intersection_Consensus_DMPs.csv` |
| Native dual-pipeline consensus DMPs | `6,269` | `Intersection_Native_Consensus_DMPs.csv` |
| DMR rows | Minfi `4,537`; Sesame strict `5,076`; Sesame native `5,530` | `Minfi_DMRs.csv`, `Sesame_DMRs.csv`, `Sesame_Native_DMRs.csv` |
| Primary branch / result mode | `Minfi`; `tier3_ineligible` | `summary.json` |
| CRF sample tier | `moderate` | `summary.json`, `analysis_parameters.json` |
| Tissue / array / preset | `DLPFC`; `EPIC`; `conservative` | `analysis_parameters.json` |

## QC and guardrails

The result is strong as a cell-type reference axis but should remain guarded in text and figures:

- All three tested branches triggered the lambda guard at the configured threshold `1.5`: Minfi lambda `1.5515`, Sesame strict lambda `1.5940`, and Sesame native lambda `1.8660` (`Minfi_Metrics.csv`, `Sesame_Metrics.csv`, `Sesame_Native_Metrics.csv`).
- The primary result mode is `tier3_ineligible`, with `Sentrix_ID` recorded as the tier-3 batch field (`summary.json`, branch metrics files). This makes the lane suitable for a positive-control/reference narrative, not for an unqualified discovery claim.
- Correction adequacy is mixed: batch removal score `1.0`, preservation score `0.8162`, calibration score `0`, and CAI `0.5357` (`Correction_Adequacy_Summary.csv`).
- Permutation null FPR is `0`, but the median null KS p-value is very small (`8.02e-16`), reinforcing the need to describe calibration caveats rather than treating nominal significance as sufficient.
- RefFree cell latent variables are highly associated with the cell-type contrast for latent factors 1 and 2 (eta-squared `0.7103` and `0.9053`; `Cell_vs_Group_Association.csv`). IlluMeta drops these strongly group-associated latent variables from the adjustment set and uses weaker latent factors (`Cell_Latent4;Cell_Latent3`) plus SVs, which is appropriate for preserving the biological cell-type axis.

## Methods/YMETH fit

This lane fits the Methods/YMETH-compatible reframing because it demonstrates IlluMeta as a method for transparent reanalysis rather than a universal biomarker generator:

1. **Cell-type-aware validation axis:** the purified Neurons-vs-Microglia contrast is a direct neural cell-type comparison, complementary to bulk AD cohorts where disease, region, platform, and cell composition are entangled.
2. **Dual-pipeline reproducibility:** thousands of strict and native consensus DMPs survive both Minfi and Sesame-family processing, showing that IlluMeta can expose robust cell-type methylation structure under independent preprocessing branches.
3. **Guarded inference:** lambda guard triggers, tier-3 ineligibility, and correction adequacy scores are explicitly reported, preventing overclaiming and aligning with the paper's emphasis on QC-guarded interpretation.
4. **Reference-axis role:** the result should be used as a positive control for expected neural cell-type separability and as context for why AD methylation signals can be cohort/platform/cell-type dependent. It should not be framed as an AD-specific signature.

## Manuscript-safe phrasing

Suggested guarded language:

> As a neural cell-type reference axis, GSE306226 Neurons-vs-Microglia produced thousands of dual-pipeline consensus DMPs (strict: 6,398; native: 6,269) across 20 microglia and 20 neuron samples. Because all primary branches triggered lambda guards and the run remained tier-3 ineligible, we use this result as a QC-annotated positive control for cell-type separability rather than as an unqualified biological discovery claim.

## Integration notes for Lane 3

- Use this lane to support the statement that IlluMeta exposes cell-type-dependent methylation structure with explicit branch and QC diagnostics.
- Pair the consensus counts with the lambda/tier caveat in any figure caption or manuscript paragraph.
- Do not merge this reference axis into AD cross-cohort biomarker counts; keep it visually and textually separated from AD-vs-control cohorts.
- If included in a summary figure, label it as `Neurons vs Microglia reference` or `cell-type positive control`, not as an AD cohort.
