# GSE226298 optional peripheral cell-type readiness

## Scope and decision

- Manuscript lane: **supplement / peripheral immune cell-type stress-test**, not neural evidence and not a universal AD biomarker claim.
- Current decision: **ready for guarded optional run after GSE306227/GSE203332 are terminal or blocked**.
- Rationale: raw sample-level IDAT pairs and clean disease labels are present, but the biology is granulocyte/monocyte peripheral immune cell type, male/LOY-focused, and should not be pooled with brain or neural-cell cohorts.

## Raw-IDAT gate

- Local metadata source: `benchmarks/dementia_special_issue/meta_candidate_triage/GSE226298_object1_metadata.tsv`.
- Raw filelist: `raw_idat_filelist.tsv`.
- Verified from metadata: 212 IDAT URLs = 106 Grn/Red sample pairs; sample HEAD checks returned HTTP 200 for representative Grn/Red files on 2026-05-22.
- Download command: `bash projects/GSE226298/download_idats.sh`.

## Contrast design

| contrast config | test | control | n test | n control | interpretation |
|---|---:|---:|---:|---:|---|
| `configure_AD_vs_Control_granulocytes.tsv` | AD | Control | 39 | 26 | Peripheral granulocyte AD/control stress-test. |
| `configure_AD_vs_Control_monocytes.tsv` | AD | Control | 24 | 17 | Peripheral monocyte AD/control stress-test. |

Both contrasts include `gender`, `loy_status`, `donor_id`, `Sentrix_ID`, and `Sentrix_Position` for audit/covariate inspection. Because LOY is a strong study feature, LOY balance must be checked in QC and any signal interpreted as peripheral/LOY-sensitive unless branch metrics prove otherwise.

## Suggested run commands

```bash
bash projects/GSE226298/download_idats.sh
Rscript r_scripts/analyze.R   --config projects/GSE226298/configure_AD_vs_Control_granulocytes.tsv   --idat_dir projects/GSE226298/idat   --out projects/GSE226298/AD_vs_Control_granulocytes_results   --group_con Control --group_test AD --tissue Blood --beginner_safe
Rscript r_scripts/analyze.R   --config projects/GSE226298/configure_AD_vs_Control_monocytes.tsv   --idat_dir projects/GSE226298/idat   --out projects/GSE226298/AD_vs_Control_monocytes_results   --group_con Control --group_test AD --tissue Blood --beginner_safe
```

## Blockers / guardrails

- Do not start before higher-priority GSE306227/GSE203332 lanes are underway, terminal, or explicitly blocked.
- Do not call results brain/neural replication; use only as peripheral immune cell-type evidence.
- Keep LOY/gender sensitivity explicit in any report and prefer supplement/deferred placement unless results are exceptionally clean and concordant.
