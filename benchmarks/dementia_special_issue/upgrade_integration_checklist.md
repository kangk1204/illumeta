# IlluMeta dementia upgrade integration checklist

Generated: 2026-05-22T08:54:01.986795+00:00

## Scope guard

- Worker-4 owns integration/verification notes and derived summary artifacts only; live result directories remain owned by their dataset workers.
- Do not update `cell_type_reframe_summary.tsv`, reports, figures, source-data workbooks, or manuscript text with GSE306227/GSE203332/GSE226298/GSE153712 until the responsible task reaches a terminal state with auditable result evidence.
- Preserve the Methods/YMETH workflow framing: reproducible, guarded, negative, blocked, and deferred evidence must remain distinct; no universal AD biomarker claim is allowed.

## Current team gate

| Task | Owner | Status | Lease | Integration action |
| --- | --- | --- | --- | --- |
| 1: GSE306227 independent CEAM neural cell-type reference analysis | worker-1 | in_progress | 2026-05-22T09:05:04.897Z | wait for terminal evidence |
| 2: GSE203332 frontal cortex AD specificity/stress-test lane | worker-2 | completed |  | integrate completed evidence |
| 3: Optional extension readiness: GSE226298 and GSE153712 | worker-3 | completed |  | record readiness/deferred evidence; do not pool with brain/neural cohorts |
| 4: Integration and verification for manuscript upgrade evidence | worker-4 | completed in task state; checklist committed as durable handoff |  | maintain this checklist until task 1 terminal evidence arrives |

## Baseline derived evidence already reproducible

| Axis | Source | n_control | n_test | Strict DMPs | Native DMPs | Lambda guard | Manuscript use |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| AD public-IDAT dementia cohorts | `benchmarks/dementia_special_issue/cross_cohort/cohort_summary.tsv` | 771 | 883 | 2562 | 3386 | 3/6 cohorts triggered primary lambda guard | reproducible but heterogeneous workflow evidence: Six AD/dementia cohorts are workflow evidence with heterogeneity: 3 strict-zero and 2 native-zero consensus cohorts, plus only 4 nonzero pairwise CpG-overlap rows. |
| GSE66351 bulk AD sensitivity lane | `projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare/summary.json` | 52 | 76 | 0 | 10 | primary=ok; native branch is treated as guarded sensitivity evidence | guarded/negative sensitivity evidence: Bulk frontal/temporal cortex AD run supports heterogeneity/guardrail narrative rather than a universal AD CpG panel. |
| GSE306226 Neurons_vs_Microglia reference axis | `projects/GSE306226/Neurons_vs_Microglia_results` | 20 | 20 | 6398 | 6269 | primary=triggered; Minfi=triggered; SeSAMe=triggered; native=triggered | guarded neural cell-type reference axis: Neural cell-type positive-control/reference axis (Microglia n=20; Neurons n=20); QC passed 40/40 samples, but lambda guard prevents unqualified biological-discovery claims. Correction preservation score=0.816. |


## Newly terminal upstream evidence captured

### GSE203332 frontal-cortex mixed neurodegeneration stress-test (task 2)

- Current terminal status: metadata/raw-gate readiness completed; no full IlluMeta result directory or `summary.json` is available yet in leader/worktree-4, so it must not be counted as a completed cohort.
- Read-only evidence source: `worktrees/worker-2/projects/GSE203332/gse203332_lane_status.md`.
- Raw/metadata gate: 492 GEO metadata rows on GPL21145 EPIC; 984 raw IDAT links = 492 sample pairs; NCBI filelist includes `GSE203332_RAW.tar` plus 984 IDAT files.
- Primary locked contrast: neuropathological `AD` vs `CONTR`, encoded as `AD` vs `Control`, after excluding technical replicates and failed genotyping QC.
- Primary counts after QC gate: AD n=61; Control n=74; usable QC-gated samples n=442.
- Planned covariates/batch fields: age, sex, post-mortem delay (`pmd_min`), amplification plate, slide; pathology-linked fields are not adjustment covariates.
- Manuscript-use decision now: readiness/contrast-design evidence only. If full IlluMeta completes later, use it as guarded frontal-cortex neurodegeneration specificity/stress-test evidence, not universal AD biomarker validation.


### GSE226298 and GSE153712 optional extension readiness (task 3)

- Current terminal status: readiness/configuration completed; no full optional IlluMeta run was launched and no optional `summary.json` exists, so neither lane should be counted in aggregate AD/dementia cohort totals.
- GSE226298 raw gate: 212 IDAT URLs = 106 sample pairs, with representative HEAD checks reported HTTP 200 on 2026-05-22.
- GSE226298 contrast designs: granulocytes AD n=39 vs Control n=26; monocytes AD n=24 vs Control n=17. Manuscript-use decision: guarded optional peripheral immune cell-type stress-test only, explicitly LOY/gender-sensitive and not neural evidence.
- GSE153712 raw gate: 1452 IDAT URLs = 726 sample pairs, with representative HEAD checks reported HTTP 200 on 2026-05-22.
- GSE153712 contrast designs: whole-blood AD n=161 vs Control n=471; MCI n=94 vs Control n=471. Manuscript-use decision: run-deferred large peripheral whole-blood robustness supplement; keep AD and MCI separate and do not pool with brain/neural cohorts.

## Artifact verification checklist for each new terminal lane

For every completed or failed dataset lane, record the dataset-specific decision before touching manuscript-facing outputs:

| Evidence class | Required checks | Manuscript handling |
| --- | --- | --- |
| Reproducible | `summary.json`; strict/native consensus CSVs; branch metrics; DMR rows; primary/native lambda guard `ok` or explicitly justified; source-data row added | Can support workflow reproducibility in a dataset-specific context |
| Guarded | Output artifacts exist but lambda guard triggered, branch concordance weak, or contrast is sensitivity-only | State as guarded support/stress-test, not biomarker validation |
| Negative | Full run completed with zero strict/native consensus or no defensible overlap | Keep as negative evidence demonstrating context dependence |
| Blocked | Metadata/raw gate/runtime/power/resources prevent a valid run | Include blocker reason and omit from quantitative pooled counts |
| Deferred | Optional lane not started because higher-priority lanes remain unresolved | Mention only as readiness/deferred supplement if documented |

## Source and dashboard consistency gates

- `cross_cohort/cohort_summary.tsv`, `branch_metrics.tsv`, and `consensus_pairwise_overlap.tsv` must be regenerated or explicitly left unchanged when new AD/dementia runs are integrated.
- `cell_type_reframe_summary.tsv` must exactly match `build_cell_type_reframe_outputs.py` output for the chosen source root before report/manuscript updates.
- FIG2 source-data/workbook/dashboard files must include any new cross-cohort rows before manuscript claims mention new aggregate counts.
- FIG3 summary/manifest/legend must retain the cell-type-reference vs AD-cohort distinction and guardrail language.
- For each result directory, verify `summary.json`, `Intersection_Consensus_DMPs.csv`, `Intersection_Native_Consensus_DMPs.csv`, `Intersection_Comparison_Metrics.csv`, `Intersection_Native_Comparison_Metrics.csv`, branch `*Metrics.csv`, and lambda/QC fields where present.

## Current file-existence audit

| Artifact | Status | Path |
| --- | --- | --- |
| cross_cohort/cohort_summary.tsv | present | `benchmarks/dementia_special_issue/cross_cohort/cohort_summary.tsv` |
| cross_cohort/branch_metrics.tsv | present | `benchmarks/dementia_special_issue/cross_cohort/branch_metrics.tsv` |
| cross_cohort/consensus_pairwise_overlap.tsv | present | `benchmarks/dementia_special_issue/cross_cohort/consensus_pairwise_overlap.tsv` |
| GSE306226 summary.json | present | `projects/GSE306226/Neurons_vs_Microglia_results/summary.json` |
| GSE306226 strict consensus | present | `projects/GSE306226/Neurons_vs_Microglia_results/Intersection_Consensus_DMPs.csv` |
| GSE306226 native consensus | present | `projects/GSE306226/Neurons_vs_Microglia_results/Intersection_Native_Consensus_DMPs.csv` |
| GSE66351 neuron-occipital summary.json | present | `projects/GSE66351/AD_vs_CTRL_neuron_occipital_results_direct_run5_tier3skip/summary.json` |
| GSE66351 bulk summary.json | present | `projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare/summary.json` |
| GSE306227 terminal result scan | project dir present | no terminal `summary.json` found yet |
| GSE203332 terminal result scan | readiness artifacts present in worker-2 worktree; no terminal run summary | `projects/GSE203332/gse203332_lane_status.md`; no `summary.json` yet |
| GSE226298 terminal result scan | readiness artifacts present in worker-3 worktree; no terminal run summary | `projects/GSE226298/README_readiness.md`; granulocytes 39/26, monocytes 24/17 |
| GSE153712 terminal result scan | readiness artifacts present in worker-3 worktree; no terminal run summary | `projects/GSE153712/README_readiness.md`; AD 161/471, MCI 94/471 |


## Native subagent findings integrated

- Subagents spawned: 2 (`019e4ee0-cdd7-7260-99ca-3b9978a49ceb` test/coverage probe; `019e4ee0-e321-7a30-9e18-0c2a21955a04` change-slice/blocker probe).
- Existing coverage to reuse: `tests/test_benchmark_table.py` for artifact consumers, `tests/test_dashboard_warnings.py`, `tests/test_integration.py`, `tests/test_preflight.py`, and `tests/test_run_smoke_pipeline.py`; downstream readers include `scripts/build_benchmark_table.py`, `scripts/build_intersection_report.py`, `scripts/build_signal_preservation_report.py`, `scripts/prepare_geo_submission.py`, and `scripts/build_supplementary_data_docx.py`.
- Regression gaps to preserve in final audit: no direct schema/golden test for `summary.json`; no automatic row-count parity check between `summary.json` and strict/native consensus CSVs; no direct lambda/status parity check between `*_Metrics.csv`, dashboards, and `summary.json`; no enforced source-data/manuscript drift test.
- Change-slice finding: final cross-cohort and cell-type-reframe edits are blocked until task 1 reaches a terminal state; task 2 contributes readiness/contrast-design evidence and task 3 contributes optional readiness/deferred evidence only at this point.

## Ready-to-run verification commands

```bash
python -m py_compile benchmarks/dementia_special_issue/build_cell_type_reframe_outputs.py
python benchmarks/dementia_special_issue/build_cell_type_reframe_outputs.py --source-root /home/keunsoo/Projects/23_illumeta --output-root /tmp/illumeta_celltype_smoke
diff -u benchmarks/dementia_special_issue/cell_type_reframe_summary.tsv /tmp/illumeta_celltype_smoke/benchmarks/dementia_special_issue/cell_type_reframe_summary.tsv
pytest tests/test_dashboard_warnings.py tests/test_integration.py
```

## Stop condition for final integration

Final manuscript-facing integration is safe only after task 1 is completed or failed with explicit evidence and the leader has merged terminal artifacts from task 2/3 worktrees. Until then, this checklist is the derived artifact of record: GSE203332, GSE226298, and GSE153712 are readiness/deferred evidence only, and GSE306227 remains pending for aggregate claims.
