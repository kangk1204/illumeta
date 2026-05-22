# Dementia special issue execution report

Generated: 2026-05-06 Asia/Seoul; updated: 2026-05-08 Asia/Seoul

## Scope

Direct autonomous execution is advancing the IlluMeta dementia special issue workflow from dataset triage into full public-IDAT analyses. The current completed evidence now includes a full IlluMeta dual-pipeline peripheral blood proof-of-work cohort (`GSE208623`), a full covariate-aware brain validation cohort (`GSE125895`), a large covariate-aware brain extension cohort (`GSE134379`), a recent EPIC prefrontal-cortex brain extension cohort (`GSE284764`), a BS-only brain 450K extension (`GSE105109`), and a frontal/temporal cortex bulk brain extension (`GSE66351`).

The special issue URL was rechecked through the ScienceDirect text mirror on 2026-05-07. It maps to the Methods special issue "Advancing Mental & Neural Health Assessment - YMETH", last updated 2025-12-29, with the stated scope "AI-driven and computational methodologies for digital healthcare applications in mental and neural health assessment." This is a plausible venue fit for an IlluMeta dementia methylation application note, but final article-type/formatting constraints should still be checked against the Methods author guide before submission.

## Commands and evidence

- `python illumeta.py search -k dementia -o search/dementia_idat_candidates.tsv --retmax 80 --sleep 0.1` produced 24 IDAT-backed candidates.
- `python illumeta.py search -k Alzheimer -o search/alzheimer_idat_candidates.tsv --retmax 120 --sleep 0.1` produced 16 IDAT-backed candidates.
- `./scripts/install_full.sh --preflight` passed.
- `./scripts/install_full.sh` completed and `./scripts/illumeta doctor` passed core package checks; optional saliva/cord-blood references remained unavailable.
- `./scripts/illumeta download GSE208623 -o projects/GSE208623` downloaded and extracted 120 IDAT files for 60 samples.
- `./scripts/illumeta download GSE125895 -o projects/GSE125895` found/extracted 538 IDAT files for 269 samples.
- `GSE125895` IDATs required `LC_ALL=C LANG=C` for `minfi::read.metharray`; under the default UTF-8 locale, `illuminaio` failed with `invalid UTF-8 input in readChar()`. `gzip -t` and a direct `LC_ALL=C` minfi read confirmed the raw IDATs were valid.
- `./scripts/illumeta download GSE134379 -o projects/GSE134379` downloaded/extracted 1,616 IDAT files for 808 samples using the direct RAW archive fallback.
- `GSE134379` direct `LC_ALL=C` minfi import was verified on a representative IDAT pair, and the full IlluMeta run completed with 363 result files, `summary.json`, `methods.md`, and a dashboard.
- `./scripts/illumeta download GSE284764 -o projects/GSE284764 --platform GPL21145` downloaded/extracted 484 IDAT files for 242 samples using the direct RAW archive fallback.
- `GSE284764` direct `LC_ALL=C` minfi import was verified on a representative IDAT pair; the full IlluMeta run completed in `projects/GSE284764/AD_vs_Ctl_PFC_results_lcC_run2` with `summary.json`, `methods.md`, `decision_ledger.tsv`, and a dashboard.
- `GSE105109` completed in `projects/GSE105109/AD_vs_Control_BS_only_results_lcC`; it produced DMR outputs in all branches but zero strict/native consensus DMPs.
- `GSE66351` completed in `projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare`; it produced zero strict consensus DMPs and 10 native consensus DMPs, with native-branch lambda warning requiring conservative interpretation.
- `benchmarks/dementia_special_issue/build_full_cross_cohort_summary.py` generated full-output cross-cohort tables for `GSE208623`, `GSE125895`, `GSE134379`, `GSE284764`, `GSE105109`, and `GSE66351`, including branch-level lambda guard fields.
- `benchmarks/dementia_special_issue/make_full_cross_cohort_figure.py` generated the refreshed six-cohort Figure 2 package and passed the `data-figure` validation helper.
- `benchmarks/dementia_special_issue/write_full_manuscript_draft.py` generated the refreshed six-cohort manuscript draft, evidence map, and checklist.
- `GSE197305` metadata-only triage completed; GEO metadata lacks a sample-level AD/control or neuropathology contrast, so the full 1,221-sample run was not started.
- ScienceDirect source check: `https://www.sciencedirect.com/special-issue/10GC6QLG9C4` via `https://r.jina.ai/http://www.sciencedirect.com/special-issue/10GC6QLG9C4`.

## Dataset triage

Recommended execution order:

1. `GSE208623`: first peripheral blood EPIC proof-of-work, 60 samples, balanced MCI/AD/control design.
2. `GSE125895`: preferred brain validation, 450K, 269 samples, multiple brain regions.
3. `GSE134379`: larger AD brain fallback/extension, 450K, 808 arrays.
4. `GSE284764`: recent EPIC prefrontal-cortex brain extension, 242 raw-IDAT samples, 217 retained after mixed-array-size guarding.
5. `GSE105109`: BS-only brain 450K extension, completed with zero consensus DMPs and DMR outputs.
6. `GSE66351`: frontal/temporal cortex bulk brain 450K extension, completed with native-only guarded consensus evidence.
7. `GSE197305`: large EPIC extension, metadata-blocked unless an external BDR phenotype table can be mapped to sample IDs.

Triage table: `benchmarks/dementia_special_issue/dataset_triage.tsv`

GSE197305 metadata triage note: `benchmarks/dementia_special_issue/gse197305_metadata_triage.md`

## Completed result: GSE208623

Contrast: `AD_vs_Control`, 20 AD vs 20 controls.

Primary output directory: `projects/GSE208623/AD_vs_Control_full_long_results`

Dashboard: `projects/GSE208623/AD_vs_Control_full_long_results_index.html`

Summary:

- Minfi DMPs: 684 up, 4,662 down.
- SeSAMe strict DMPs: 982 up, 1,928 down.
- SeSAMe native DMPs: 1,533 up, 2,850 down.
- Minfi/SeSAMe strict consensus: 483 up, 1,783 down; 2,266 consensus probes.
- Minfi/SeSAMe native consensus: 540 up, 2,508 down; 3,048 consensus probes.
- DMR outputs completed for all three branches: Minfi 4,685 rows, SeSAMe strict 5,784 rows, SeSAMe native 8,460 rows.
- Primary branch: Minfi; primary mode: `tier3_ineligible`; lambda guard triggered; CRF sample tier: moderate.

The earlier apparent stall after Minfi/PVCA was not a hard failure. It was a long CPU-bound late EWAS/reporting phase that produced no visible files for about 35-40 minutes and then completed.

## Completed result: GSE125895

Contrast: `AD_vs_Control`, all brain regions combined, 82 AD vs 187 controls.

Config: `projects/GSE125895/configure_AD_vs_Control_all_regions.tsv`

Command shape:

```bash
LC_ALL=C LANG=C ./scripts/illumeta analysis \
  -c projects/GSE125895/configure_AD_vs_Control_all_regions.tsv \
  --idat-dir projects/GSE125895/idat \
  --group_con Control \
  --group_test AD \
  --output projects/GSE125895/AD_vs_Control_all_regions_results_lcC \
  --tier3-on-fail skip \
  --include-covariates source_name_ch1,characteristics_ch1.6,characteristics_ch1.7,characteristics_ch1.8
```

Primary output directory: `projects/GSE125895/AD_vs_Control_all_regions_results_lcC`

Dashboard: `projects/GSE125895/AD_vs_Control_all_regions_results_lcC_index.html`

Metadata structure:

- Brain regions: ERC 69, DLPFC 68, CRB 67, HIPPO 65.
- Disease labels: Control 187, Alzheimer 82.
- Region-level disease split: CRB 24 AD / 43 Control; DLPFC 21 AD / 47 Control; ERC 20 AD / 49 Control; HIPPO 17 AD / 48 Control.
- Covariates included: brain region (`source_name_ch1`), age (`characteristics_ch1.6`), sex (`characteristics_ch1.7`), race (`characteristics_ch1.8`).

Summary:

- Minfi DMPs: 233 up, 163 down.
- SeSAMe strict DMPs: 191 up, 138 down.
- SeSAMe native DMPs: 210 up, 149 down.
- Minfi/SeSAMe strict consensus: 166 up, 118 down; 284 consensus probes.
- Minfi/SeSAMe native consensus: 183 up, 133 down; 316 consensus probes.
- DMR outputs completed for all three branches: Minfi 2,666 rows, SeSAMe strict 2,162 rows, SeSAMe native 2,542 rows.
- Primary branch: Minfi; primary mode: `tier3_ineligible`; lambda guard status: ok; CRF sample tier: large.

Full table row counts:

- `Minfi_DMPs_full.csv`: 308,879 rows.
- `Sesame_DMPs_full.csv`: 234,783 rows.
- `Sesame_Native_DMPs_full.csv`: 303,833 rows.
- `Intersection_Consensus_DMPs.csv`: 284 rows.
- `Intersection_Native_Consensus_DMPs.csv`: 316 rows.

## Completed result: GSE134379

Contrast: `AD_vs_Control`, all sampled brain regions combined, 450 AD vs 358 controls.

Config: `projects/GSE134379/configure_AD_vs_Control_all_regions.tsv`

Command shape:

```bash
ILLUMETA_TIMEOUT=172800 LC_ALL=C LANG=C ./scripts/illumeta analysis \
  -c projects/GSE134379/configure_AD_vs_Control_all_regions.tsv \
  --idat-dir projects/GSE134379/idat \
  --group_con Control \
  --group_test AD \
  --output projects/GSE134379/AD_vs_Control_all_regions_results_lcC \
  --tier3-on-fail skip \
  --include-covariates brain_region,age,sex,plate
```

Primary output directory: `projects/GSE134379/AD_vs_Control_all_regions_results_lcC`

Dashboard: `projects/GSE134379/AD_vs_Control_all_regions_results_lcC_index.html`

Metadata structure:

- Brain regions: CBL 404, MTG 404.
- Disease labels: AD 450, Control 358.
- Region-level disease split: CBL 225 AD / 179 Control; MTG 225 AD / 179 Control.
- Sex split: F 248 AD / 152 Control; M 202 AD / 206 Control.
- Covariates included: brain region, age, sex, and plate.

Summary:

- Minfi DMPs: 1 up, 0 down.
- SeSAMe strict DMPs: 0 up, 0 down.
- SeSAMe native DMPs: 0 up, 0 down.
- Minfi/SeSAMe strict consensus: 0 up, 0 down; 0 consensus probes.
- Minfi/SeSAMe native consensus: 0 up, 0 down; 0 consensus probes.
- DMR outputs completed for all three branches: Minfi 3,349 rows, SeSAMe strict 3,178 rows, SeSAMe native 3,229 rows.
- Primary branch: Minfi; primary mode: `tier3_ineligible`; lambda guard triggered; CRF sample tier: large.

Full table row counts:

- `Minfi_DMPs_full.csv`: 306,356 rows.
- `Sesame_DMPs_full.csv`: 287,419 rows.
- `Sesame_Native_DMPs_full.csv`: 289,642 rows.
- `Intersection_Consensus_DMPs.csv`: 0 rows.
- `Intersection_Native_Consensus_DMPs.csv`: 0 rows.

Interpretation note: `GSE134379` is technically complete, but the zero consensus-DMP result and triggered lambda guard mean it should not be used as a clean positive replication cohort. It is still useful as a large-cohort stress test and DMR/reference result, while the stronger manuscript-level DMP replication evidence currently comes from `GSE208623` and `GSE125895`.

## Completed result: GSE284764

Contrast: `AD_vs_Ctl`, prefrontal cortex, 119 AD vs 98 controls retained after mixed-array-size guarding.

Config: `projects/GSE284764/configure_AD_vs_Ctl_PFC.tsv`

Command shape:

```bash
R_LIBS_USER=/home/keunsoo/Projects/23_illumeta/.r-lib/R-4.4 \
ILLUMETA_TIMEOUT=172800 LC_ALL=C LANG=C ./scripts/illumeta analysis \
  -c /home/keunsoo/Projects/23_illumeta/projects/GSE284764/configure_AD_vs_Ctl_PFC.tsv \
  --idat-dir /home/keunsoo/Projects/23_illumeta/projects/GSE284764/idat \
  --group_con Ctl \
  --group_test AD \
  --output /home/keunsoo/Projects/23_illumeta/projects/GSE284764/AD_vs_Ctl_PFC_results_lcC_run2 \
  --tier3-on-fail skip \
  --include-covariates age,sex,brain_bank
```

Primary output directory: `projects/GSE284764/AD_vs_Ctl_PFC_results_lcC_run2`

Dashboard: `projects/GSE284764/AD_vs_Ctl_PFC_results_lcC_run2_index.html`

Metadata and QC structure:

- GEO metadata: 131 AD / 111 Ctl, prefrontal cortex, EPIC (`GPL21145`).
- IDAT import: 484 IDAT files, 242 paired samples.
- Mixed-array-size safeguard dropped 25 samples with non-modal IDAT array size.
- Final analysis sample count: 217 total; 119 AD and 98 Ctl.
- Covariates included: age, sex, and brain bank.
- Braak stage was not used as a covariate because it is disease/pathology-linked.

Summary:

- Minfi DMPs: 13 up, 1 down.
- SeSAMe strict DMPs: 26 up, 2 down.
- SeSAMe native DMPs: 27 up, 2 down.
- Minfi/SeSAMe strict consensus: 11 up, 1 down; 12 consensus probes.
- Minfi/SeSAMe native consensus: 11 up, 1 down; 12 consensus probes.
- DMR outputs completed for all three branches: Minfi 6,770 rows, SeSAMe strict 16,168 rows, SeSAMe native 18,327 rows.
- Primary branch: Minfi; primary mode: `tier3_ineligible`; lambda guard triggered; CRF sample tier: large.

Full table row counts:

- `Minfi_DMPs_full.csv`: 671,166 rows.
- `Sesame_DMPs_full.csv`: 391,704 rows.
- `Sesame_Native_DMPs_full.csv`: 457,109 rows.
- `Intersection_Consensus_DMPs.csv`: 12 rows.
- `Intersection_Native_Consensus_DMPs.csv`: 12 rows.

Interpretation note: `GSE284764` is a valuable recent raw-IDAT EPIC PFC extension with nonzero consensus DMPs and six strict/native consensus CpGs shared with `GSE125895`. Because the lambda guard triggered, it should be presented as guarded supporting evidence, not as an unqualified positive replication.

## Cell-type-aware neuroepigenomics reframe update

Generated: 2026-05-22T13:42:38.968915+00:00

### Source-grounded conclusion

The current package is strongest as a Methods/YMETH-compatible neuroepigenomics workflow paper, not as a universal Alzheimer disease biomarker paper. The completed AD/dementia public-IDAT analyses cover 771 controls and 883 AD/test samples across blood, bulk brain, EPIC brain, BS-only brain, and frontal/temporal cortex contexts. They produced 2562 strict and 3386 native consensus DMPs in aggregate, but the evidence is heterogeneous: 3/6 cohorts triggered primary lambda guard; six AD/dementia cohorts are workflow evidence with heterogeneity: 3 strict-zero and 2 native-zero consensus cohorts, plus only 4 nonzero pairwise CpG-overlap rows.

### Guarded AD evidence

- `GSE125895` remains the cleanest brain AD anchor because it has nonzero strict/native consensus and primary lambda guard `ok`.
- `GSE284764` provides recent EPIC prefrontal-cortex support, but its primary lambda guard triggered, so it remains guarded support.
- `GSE134379` and `GSE105109` completed branch-level DMP/DMR outputs but produced no consensus DMPs.
- `GSE66351` should be described as a bulk-cortex sensitivity lane: 0 strict and 10 native consensus DMPs; Bulk frontal/temporal cortex AD run supports heterogeneity/guardrail narrative rather than a universal AD CpG panel.

### Neural cell-type reference axis

`GSE306226/Neurons_vs_Microglia_results` provides the local positive-control/reference axis for cell-type-aware framing. It compares 20 neurons with 20 microglia, passed sample QC, and produced 6398 strict plus 6269 native Minfi/SeSAMe consensus DMPs. Branch concordance was strict logFC r=0.791, Jaccard=0.302; native logFC r=0.761, Jaccard=0.304. Because primary=triggered; Minfi=triggered; SeSAMe=triggered; native=triggered, it should be used as a reference axis demonstrating that IlluMeta can recover a large neural cell-type contrast, not as an unqualified discovery claim.


### Guarded GSE306227 extension

`GSE306227/Neurons_vs_Microglia_results` is retained as a guarded independent sorted-brain extension, not as a core Figure 3 axis in the current package. It has summary-level evidence for 324706 strict and 334911 native consensus DMPs across 18 microglia/control-axis and 19 neuron/test-axis samples, with primary result mode `tier3_ineligible` and lambda guard `triggered`. However, branch metrics present=False, cell summary present=False, and the run log records RefFreeEWAS unavailable=True. Therefore it should be cited only as guarded extension evidence until branch-level design and cell-adjustment artifacts are synced.


### Manuscript claim guardrail

Defensible central claim: IlluMeta exposes context dependence in public neurodegeneration methylation data by making raw-IDAT import, branch-specific covariate/cell/SV adjustment, dual-pipeline agreement, consensus DMP/DMR outputs, and lambda/QC guards auditable. The manuscript should explicitly avoid phrases such as "universal AD biomarker" or "validated AD signature" unless tied to a specific cohort/context and guard status.

## Figure package

The existing proof-of-work Figure 1 package predates the full dual-pipeline completion and should be replaced or demoted before manuscript submission.

Legacy proof-of-work package:

- PNG: `benchmarks/dementia_special_issue/figures/FIG1/figure1_dementia_illumeta_proof.png`
- PDF: `benchmarks/dementia_special_issue/figures/FIG1/figure1_dementia_illumeta_proof.pdf`
- SVG: `benchmarks/dementia_special_issue/figures/FIG1/figure1_dementia_illumeta_proof.svg`
- Source workbook: `benchmarks/dementia_special_issue/figures/FIG1/figure1_dementia_illumeta_proof_source_data.xlsx`
- Legend DOCX: `benchmarks/dementia_special_issue/figures/FIG1/figure1_dementia_illumeta_proof_legend.docx`
- Manifest: `benchmarks/dementia_special_issue/figures/FIG1/figure1_dementia_illumeta_proof_manifest.json`

Validation previously passed:

- `python /home/keunsoo/.codex/skills/data-figure/scripts/validate_figure_package.py ...`
- `python -m py_compile benchmarks/dementia_special_issue/make_dementia_special_issue_outputs.py`

Current full-output cross-cohort package:

- Script: `benchmarks/dementia_special_issue/make_full_cross_cohort_figure.py`
- PNG: `benchmarks/dementia_special_issue/figures/FIG2/figure2_dementia_full_cross_cohort_evidence.png`
- PDF: `benchmarks/dementia_special_issue/figures/FIG2/figure2_dementia_full_cross_cohort_evidence.pdf`
- SVG: `benchmarks/dementia_special_issue/figures/FIG2/figure2_dementia_full_cross_cohort_evidence.svg`
- Source workbook: `benchmarks/dementia_special_issue/figures/FIG2/figure2_dementia_full_cross_cohort_evidence_source_data.xlsx`
- Legend DOCX: `benchmarks/dementia_special_issue/figures/FIG2/figure2_dementia_full_cross_cohort_evidence_legend.docx`
- Manifest: `benchmarks/dementia_special_issue/figures/FIG2/figure2_dementia_full_cross_cohort_evidence_manifest.json`

Validation completed:

- `python -m py_compile benchmarks/dementia_special_issue/build_full_cross_cohort_summary.py`
- `python benchmarks/dementia_special_issue/build_full_cross_cohort_summary.py`
- `python -m py_compile benchmarks/dementia_special_issue/make_full_cross_cohort_figure.py`
- `python benchmarks/dementia_special_issue/make_full_cross_cohort_figure.py`
- `python /home/keunsoo/.codex/skills/data-figure/scripts/validate_figure_package.py ...` returned `validation passed`.
- Visual QA confirmed the regenerated combined PNG is nonblank and includes six cohorts plus the GSE125895/GSE284764 shared-CpG panel.

## Cross-cohort summary artifacts

Output directory: `benchmarks/dementia_special_issue/cross_cohort`

- `cohort_summary.tsv`: completed-cohort sample sizes, consensus counts, DMR row counts, lambda guard status, and dashboard paths.
- `branch_metrics.tsv`: Minfi/SeSAMe logFC correlations and significant-set overlaps.
- `consensus_pairwise_overlap.tsv`: strict/native CpG overlap across cohorts.
- `top_consensus_cpgs.tsv`: top consensus DMPs for cohorts with nonzero consensus sets.
- `top_dmr_regions.tsv`: top DMR rows across branches and cohorts.
- `dementia_full_cross_cohort_source_data.xlsx`: source workbook for downstream manuscript tables/figures.

Key cross-cohort interpretation:

- `GSE208623` and `GSE125895` both provide nonzero dual-pipeline consensus DMP evidence.
- `GSE208623` and `GSE125895` share four strict/native consensus CpGs: `cg00709979`, `cg04315214`, `cg14622549`, and `cg20240347`; all four have positive mean logFC in both cohorts.
- `GSE284764` adds a recent EPIC PFC brain extension with 12 strict/native consensus DMPs.
- `GSE125895` and `GSE284764` share six strict/native consensus CpGs: `cg05066959`, `cg05810363`, `cg07584855`, `cg12309456`, `cg17104258`, and `cg18102633`; all six are directionally concordant across the two brain cohorts.
- `GSE134379` has high Minfi/SeSAMe logFC correlation but no consensus DMPs, so its safest use is as a large-cohort DMR/stress-test result rather than a positive DMP-replication claim.
- `GSE105109` completed as a BS-only brain extension with zero strict/native consensus DMPs and DMR outputs in all branches.
- `GSE66351` completed as a frontal/temporal cortex extension with zero strict consensus DMPs and 10 native consensus DMPs; the native branch triggered lambda guard and should be treated as sensitivity evidence only.

## Draft

Legacy draft manuscript scaffold:

- `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_draft.md`

This draft is now stale because it describes fallback Minfi-only evidence.

Current full-output manuscript draft:

- Script: `benchmarks/dementia_special_issue/write_full_manuscript_draft.py`
- Markdown: `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_draft.md`
- DOCX: `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_draft.docx`
- Evidence map: `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_evidence_map.tsv`
- Reference map: `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_reference_map.tsv`
- Highlights file: `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_highlights.md`
- Submission checklist: `benchmarks/dementia_special_issue/manuscript/dementia_special_issue_submission_checklist.md`

Validation completed:

- `python -m py_compile benchmarks/dementia_special_issue/write_full_manuscript_draft.py`
- `python benchmarks/dementia_special_issue/write_full_manuscript_draft.py`
- `pandoc benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_draft.md -o benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_draft.docx`
- `unzip -t benchmarks/dementia_special_issue/manuscript/dementia_special_issue_full_draft.docx` reported no compressed-data errors.
- Manuscript text checks confirmed the expected key claims: six completed cohorts, `2,266 strict`, `284 strict`, `12 strict`, 10 native GSE66351 consensus DMPs, the four blood/brain shared CpGs, the six GSE125895/GSE284764 shared CpGs, zero GSE134379/GSE105109 consensus DMPs, and Figure 2 linkage.
- Citation/submission pass on 2026-05-08 confirmed the abstract is 233 words,
  the numbered reference list has 23 sequentially cited references, all four
  highlights are under 85 characters, and the regenerated DOCX passes `unzip -t`.
- PubMed/DOI metadata were refreshed for GEO/source studies and method
  references. `GSE284764` now cites the 2026 peer-reviewed Acta Neuropathologica
  article rather than relying on the 2025 bioRxiv preprint; `GSE134379` is cited
  as a GEO dataset reference because no PubMed ID is listed in the GEO summary.
- Methods author-guide check changed the draft article type from an application
  note to a full-length research article / methods application. The main
  remaining venue risk is scope: Methods states that computational, AI, and
  machine-learning approaches should uncover novel biological insights and must
  be experimentally validated. The current manuscript should stay framed as
  workflow assessment on public experimental methylation datasets with guarded
  independent-cohort evidence, not as a purely theoretical method or definitive
  biomarker discovery.

## Meta-analysis add-on candidate triage

The additional GEO search/triage pass screened AD- and dementia-relevant methylation
array cohorts for a defensible meta-analysis extension. The immediate objective is
not to pool every available dementia dataset, but to add compatible brain bulk
cohorts first and keep special assay/cell-type/liquid-biopsy cohorts in separate
lanes.

New planning and triage artifacts:

- Plan: `.omx/plans/dementia-meta-next-datasets-20260507.md`
- Context snapshot: `.omx/context/dementia-meta-gse-20260507T015458Z.md`
- Metadata triage script: `benchmarks/dementia_special_issue/triage_meta_gse_candidates.R`
- Candidate curation script: `benchmarks/dementia_special_issue/curate_meta_analysis_candidates.py`
- Candidate metadata summary: `benchmarks/dementia_special_issue/meta_candidate_triage/meta_gse_candidate_summary.tsv`
- Priority table: `benchmarks/dementia_special_issue/meta_candidate_triage/meta_analysis_candidate_priority.tsv`
- Run queue: `benchmarks/dementia_special_issue/meta_candidate_triage/meta_analysis_run_queue.tsv`
- First-run checklist: `benchmarks/dementia_special_issue/meta_candidate_triage/gse76105_prerun_checklist.md`
- `GSE76105` raw-gate status: `benchmarks/dementia_special_issue/meta_candidate_triage/gse76105_raw_gate_status.md`
- `GSE109627` raw-gate status: `benchmarks/dementia_special_issue/meta_candidate_triage/gse109627_raw_gate_status.md`
- `GSE80970` raw-gate status: `benchmarks/dementia_special_issue/meta_candidate_triage/gse80970_raw_gate_status.md`
- `GSE284764` run status: `benchmarks/dementia_special_issue/meta_candidate_triage/gse284764_run_status.md`

Recommended execution order:

1. `GSE76105` -- GPL13534, superior temporal gyrus bulk brain, 34 AD / 34 control.
   Raw gate status: `raw_blocked`; the series-level RAW tar contains 450K manifest
   files and no IDAT members, so this cohort is not analyzed in the raw-IDAT lane.
2. `GSE109627` -- GPL13534, middle temporal gyrus bulk brain, 91 AD / 71 control.
   Raw gate status: `raw_blocked`; the series-level RAW tar is also the same
   450K manifest bundle and contains no IDAT members.
3. `GSE80970` -- GPL13534, prefrontal cortex and superior temporal gyrus bulk brain,
   148 AD / 138 control; analyze with explicit region handling.
   Raw gate status: `raw_blocked`; official filelist shows only platform manifest
   members inside the RAW tar, so processed matrices are not substituted.
4. `GSE284764` -- GPL21145 EPIC prefrontal cortex bulk brain, 131 AD / 111 control;
   completed as the first cross-platform EPIC extension after the 450K brain lane.

Deferred or special-lane cohorts:

- `GSE105109`: high-value AD brain cohort, but BS/oxBS assays must be split before
  any IlluMeta run.
- `GSE66351`: sorted neuron/glia fractions; useful for stratified validation, not
  pooled bulk-brain meta-analysis.
- `GSE59685`: large multi-tissue cohort with an `Exclude` class; requires
  tissue-specific configuration before execution.
- `GSE226298`: peripheral granulocyte/monocyte EPIC cohort; keep as a peripheral
  validation lane.
- `GSE203332`, `GSE311578`, `GSE138597`, and `GSE272947`: defer because of mixed
  diagnoses, cfDNA biology, non-AD disease focus, or small/non-primary contrast.

Meta-analysis gate:

- `GSE76105`, `GSE109627`, and `GSE80970` can only enter the IlluMeta raw-IDAT lane
  after the series-level RAW archive is reachable and paired IDAT extraction is
  verified. Processed matrix files are not acceptable substitutes for this lane.
- Mixed-region cohorts require either region-specific within-cohort summaries or an
  explicit region covariate before formal pooling.
- Formal cross-cohort claims must carry a poolability status per cohort:
  `poolable`, `descriptive_only`, `region_split_required`, `raw_blocked`, or
  `qc_guarded`.
- `GSE134379` remains a completed large-cohort DMR/stress-test result with zero
  consensus DMPs and lambda-guarding; it must not be silently pooled as a positive
  DMP-replication cohort.
- `GSE284764` is completed and contributes guarded nonzero consensus evidence;
  use it as `qc_guarded` supporting evidence because the lambda guard triggered.

## Required next work before submission

1. Confirm author-specific metadata before submission: final author list,
   affiliations, funding wording, competing-interest declaration, and whether
   the included generative-AI disclosure should be edited to match author policy.
2. If an external BDR phenotype table is available, map it to `GSE197305`
   GSM/sample IDs and reconsider the full run; otherwise keep `GSE197305` excluded.
3. Confirm in Editorial Manager that the YMETH article collection remains
   selectable; no explicit deadline was visible on the mirrored ScienceDirect
   special-issue page.
4. Tighten the current full-output manuscript draft after the meta-analysis
   extension status is final.
5. Keep the GSE197305 exclusion explicit in limitations: metadata was checked, but
   a disease/neuropathology contrast was not exposed in GEO.
6. Do not substitute processed matrices for `GSE76105`, `GSE109627`, or `GSE80970`
   in the raw-IDAT IlluMeta evidence lane unless the manuscript explicitly opens a
   separate processed-data sensitivity lane.
