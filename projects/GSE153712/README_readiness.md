# GSE153712 optional large blood readiness

## Scope and decision

- Manuscript lane: **supplement / large whole-blood robustness**, not brain evidence and not a pooled universal AD biomarker claim.
- Current decision: **configuration-ready but run-deferred** because the raw archive is large (~10.7 GB / 1452 IDAT files) and higher-priority brain/neural lanes should not be starved.
- Primary contrast: AD vs healthy control. Secondary contrast: MCI vs healthy control.

## Raw-IDAT gate

- Local metadata source: `benchmarks/dementia_special_issue/meta_candidate_triage/GSE153712_object1_metadata.tsv`.
- Raw filelist: `raw_idat_filelist.tsv`.
- Verified from metadata: 1452 IDAT URLs = 726 Grn/Red sample pairs; sample HEAD checks returned HTTP 200 for representative Grn/Red files on 2026-05-22.
- Download command: `bash projects/GSE153712/download_idats.sh`.

## Contrast design

| contrast config | test | control | n test | n control | interpretation |
|---|---:|---:|---:|---:|---|
| `configure_AD_vs_Control.tsv` | AD | Control | 161 | 471 | Large peripheral whole-blood AD/control robustness lane. |
| `configure_MCI_vs_Control.tsv` | MCI | Control | 94 | 471 | Secondary prodromal/stage stress-test; do not mix with primary AD endpoint. |

Both contrasts include `gender`, `Sentrix_ID`, and `Sentrix_Position` for audit/covariate inspection. Whole-blood cell composition should be handled through the blood tissue/deconvolution path and interpreted separately from brain or sorted neural-cell evidence.

## Suggested run commands

```bash
# Only after resource check confirms disk/time availability and priority lanes are not starved.
bash projects/GSE153712/download_idats.sh
Rscript r_scripts/analyze.R   --config projects/GSE153712/configure_AD_vs_Control.tsv   --idat_dir projects/GSE153712/idat   --out projects/GSE153712/AD_vs_Control_results   --group_con Control --group_test AD --tissue Blood --beginner_safe
Rscript r_scripts/analyze.R   --config projects/GSE153712/configure_MCI_vs_Control.tsv   --idat_dir projects/GSE153712/idat   --out projects/GSE153712/MCI_vs_Control_results   --group_con Control --group_test MCI --tissue Blood --beginner_safe
```

## Blockers / guardrails

- Do not launch casually in this worker lane: expected download size is large and full run can compete with GSE306227/GSE203332.
- Keep as supplementary robustness unless leader explicitly promotes it after terminal brain/neural evidence exists.
- Report AD and MCI contrasts separately; do not combine AD+MCI into a single disease-positive group for manuscript claims.
