# GSE306227 IlluMeta neural cell-type reference lane status

Completed: 2026-05-22T20:08:36+09:00

## Scope

- Dataset: GSE306227 (CEAM cingulate gyrus sorted brain cell types, GPL21145 EPIC).
- Primary manuscript lane: independent sorted-brain neural cell-type reference axis comparable to GSE306226.
- Primary contrast: `Neurons` vs `Microglia` (raw-IDAT only).
- Interpretation guard: cell-type reference/positive-control axis only; not a universal AD biomarker claim.

## Raw gate and metadata

- NCBI filelist: `metadata/ncbi_filelist.txt`.
- Raw gate: 150 IDAT files in NCBI filelist = 75 sample pairs; archive listed as `GSE306227_RAW.tar`.
- GEO metadata: `metadata/GSE306227_pData_GPL21145.tsv`, 75 samples.
- Cell-type counts: Astrocytes 19, Microglia 18, Neurons 19, Oligodendrocytes 19.
- Primary selected raw IDATs: 74/74 files downloaded for 37 Neurons/Microglia samples.
- Primary counts: Microglia n=18, Neurons n=19; sex counts in selected lane: Microglia Female 9/Male 9, Neurons Female 9/Male 10.

## Commands/logs

Final successful command was direct system-R analysis after wrapper setup retry/CRF guard recovery:

```bash
Rscript r_scripts/analyze.R \
  -c /abs/path/projects/GSE306227/configure.tsv \
  --idat_dir /abs/path/projects/GSE306227/idat \
  -o /abs/path/projects/GSE306227/Neurons_vs_Microglia_results \
  --group_con Microglia \
  --group_test Neurons \
  --include_covariates sex \
  --sex_check_column sex \
  --batch_column Sentrix_ID \
  --tier3-on-fail skip \
  --cross_reactive_list /home/keunsoo/Projects/23_illumeta/references/probe_blacklists/cross_reactive_EPIC.tsv \
  --force_idat
```

Full log: `logs/run_Neurons_vs_Microglia_full.log`.

## Output evidence

Output root: `Neurons_vs_Microglia_results/`.

- `summary.json` exists.
- `Intersection_Consensus_DMPs.csv`: 324,706 strict consensus DMP rows.
- `Intersection_Native_Consensus_DMPs.csv`: 334,911 native consensus DMP rows.
- `Minfi_DMRs.csv`: 49,744 DMR rows.
- `Sesame_DMRs.csv`: 43,448 DMR rows.
- `Sesame_Native_DMRs.csv`: 44,498 DMR rows.
- Branch concordance: strict logFC correlation 0.999718, strict Jaccard 0.972599; native logFC correlation 0.999706, native Jaccard 0.977326.
- QC: 37/37 samples retained; 4,824 failed detection-P probes, 44,254 cross-reactive probes, 26,679 SNP probes, 17,508 sex-chromosome probes removed.
- Tier3/batch guard: Sentrix_ID confounding detected but ineligible for stratified/meta because per-stratum group sizes are too small (min stratum n=1); primary result mode `tier3_ineligible`.
- Lambda guard: triggered/warn in all branches (Minfi lambda 36.883; Sesame lambda 38.012; Sesame_Native lambda 35.529), so interpretation must remain guarded.

## Manuscript-use decision

Use as a strong independent sorted-brain/cell-type reference replication/extension axis for IlluMeta. The large and branch-concordant Neurons-vs-Microglia signal supports cell-type/context dependence and positive-control behavior, but lambda and Sentrix_ID guards mean it should be framed as guarded reference-axis evidence, not universal AD biomarker discovery.

## Runtime notes/blockers

- The wrapper path repeatedly entered `setup_env.R`; direct `Rscript r_scripts/analyze.R` under system R 4.5 was used after verifying required packages.
- The worker worktree lacked `references/probe_blacklists`; the run used the leader repo EPIC cross-reactive list explicitly.
- Initial DMR annotation was slow; optimized DMR annotation code in `r_scripts/analyze.R` was used for completion.
- DMR top-labelled static/PDF plots were skipped by plotting code because of a missing-value plotting condition, but DMR CSV tables were written for all branches.
