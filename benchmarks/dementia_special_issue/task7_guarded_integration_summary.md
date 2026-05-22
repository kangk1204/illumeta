# Task 7 guarded dementia upgrade integration summary

Generated: 2026-05-22T11:14:23.840596+00:00

## Integration decision

- `GSE306227` is integrated as a **completed positive-control neural cell-type reference axis / cell-type-dependence evidence**. It is not Alzheimer disease case-control evidence and must not be used as universal AD biomarker validation.
- `GSE203332` is integrated as a **blocked readiness / frontal-cortex mixed-neurodegeneration stress-test lane**. It has a valid raw gate and locked neuropathological AD-vs-CONTR contrast, but no terminal IlluMeta outputs; it must not be counted as a completed cohort.
- `GSE226298` and `GSE153712` remain optional peripheral/deferred readiness evidence from task 3 unless a future task produces terminal run artifacts.

## Artifact paths used

### GSE306227 terminal outputs

- Result directory: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results`
- `summary.json`: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/summary.json`
- Strict consensus DMPs: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Intersection_Consensus_DMPs.csv`
- Native consensus DMPs: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Intersection_Native_Consensus_DMPs.csv`
- Branch concordance: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Intersection_Comparison_Metrics.csv`, `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Intersection_Native_Comparison_Metrics.csv`
- Branch/lambda metrics: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Minfi_Metrics.csv`, `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Sesame_Metrics.csv`, `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Sesame_Native_Metrics.csv`
- DMR tables: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Minfi_DMRs.csv`, `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Sesame_DMRs.csv`, `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Sesame_Native_DMRs.csv`
- QC: `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/QC_Summary.csv`, `/home/keunsoo/Projects/23_illumeta/.omx/team/execute-the-approved-203c3847/worktrees/worker-1/projects/GSE306227/Neurons_vs_Microglia_results/Input_Group_Distribution.csv`

### GSE203332 blocker/readiness outputs

- Leader-visible project directory: `/home/keunsoo/Projects/23_illumeta/projects/GSE203332`
- Locked config: `/home/keunsoo/Projects/23_illumeta/projects/GSE203332/configure_AD_vs_Control_neuropath.tsv`
- Status: `/home/keunsoo/Projects/23_illumeta/projects/GSE203332/gse203332_lane_status.md`
- Blocker summary: `/home/keunsoo/Projects/23_illumeta/projects/GSE203332/gse203332_run_blocker_summary.tsv`
- Logs directory: `/home/keunsoo/Projects/23_illumeta/projects/GSE203332/logs`

## Verified GSE306227 terminal evidence

| Field | Value |
| --- | ---: |
| Microglia control samples | 18 |
| Neurons test samples | 19 |
| QC passed / input samples | 37 / 37 |
| Strict consensus DMP rows | 324706 |
| Native consensus DMP rows | 334911 |
| Minfi DMR rows | 49744 |
| SeSAMe DMR rows | 43448 |
| SeSAMe native DMR rows | 44498 |
| Strict logFC correlation | 0.999718 |
| Strict Jaccard overlap | 0.972599 |
| Native logFC correlation | 0.999706 |
| Native Jaccard overlap | 0.977326 |

Lambda guard summary: Minfi lambda=36.8832641566773 (triggered); SeSAMe lambda=38.0124018334161 (triggered); SeSAMe native lambda=35.529493432393 (triggered). All branches are therefore **guarded positive-control/reference evidence**, not unqualified biological-discovery evidence.

## Verified GSE203332 blocker/readiness evidence

| Field | Value |
| --- | --- |
| Primary contrast | neuropathological AD vs CONTR |
| Primary counts after QC | AD=61; Control=74 |
| Raw gate | 984 IDAT links in GEO metadata/filelist; selected primary download complete: 270 files / 135 pairs / 1.99GB |
| Run status | blocked_before_analysis |
| summary.json | missing |
| strict consensus DMPs | missing |
| native consensus DMPs | missing |
| DMR rows | missing |
| lambda guard | missing |
| Manuscript use | blocked readiness/stress-test evidence only; not a completed cohort; preserve guarded interpretation |

Blocker details: system R attempt failed loading stringi because conda ICU required CXXABI_1.3.15 unavailable in system libstdc++; conda full setup spent runtime on optional FlowSorted/cell-reference downloads not needed for tissue Auto; minimal conda setup stalled downloading Bioconductor GenomeInfoDb_1.42.3.tar.gz at 0 bytes.

## Manuscript-use classification

| Lane | Classification | Use in manuscript now | Do not claim |
| --- | --- | --- | --- |
| GSE306227 | completed guarded neural cell-type reference axis | strengthens independent CEAM/sorted-brain cell-type-dependence axis and positive-control signal recovery | universal AD biomarker, AD case-control validation |
| GSE203332 | blocked readiness/stress-test lane | cite raw gate, locked contrast, selected-IDAT download, and explicit runtime blocker only | completed cohort, negative biology, lambda/DMP/DMR result |
| GSE226298 | optional peripheral readiness | peripheral immune cell-type supplement/deferred only | neural replication or pooled brain evidence |
| GSE153712 | optional large blood readiness/deferred | peripheral whole-blood robustness supplement/deferred only | brain/neural evidence or pooled universal biomarker claim |

## Builder/report consistency decision

No generated cross-cohort or cell-type-reframe TSV/figure/manuscript counts were changed in this task because the current builder does not yet model GSE306227 as an additional cell-type-reference row and GSE203332 has no terminal run counts. This task updates only guarded integration notes/checklist. A future builder update should add a schema-supported reference-axis row for GSE306227 before changing `cell_type_reframe_summary.tsv` or FIG3 source data.
