# GSE66351 stratified lane status

Scope: Task 5 lane 1b, occipital neuron/glia stratified AD-vs-CTRL analysis.

## Verified stratified inputs

- Occipital **Neuron**: AD = 15, CTRL = 16
- Occipital **Glia**: AD = 15, CTRL = 16

## Bulk counts preserved for reference only

- Frontal cortex bulk: AD = 37, CTRL = 26
- Temporal cortex bulk: AD = 39, CTRL = 26

## Configs created

- `projects/GSE66351/configure_AD_vs_CTRL_neuron.tsv`
- `projects/GSE66351/configure_AD_vs_CTRL_glia.tsv`

## Run status

- Attempted occipital neuron analysis run(s) did not reach final result artifacts in this environment.
- The concrete blocker was the repeated first-run R dependency/setup path under `.r-lib/R-4.5`, which consumed the available execution window before `summary.json` or final result tables were written.
- No existing completed neuron-occipital result tree was found to reuse.

## No-edit note

- Do not overwrite the bulk lane outputs.
- Keep any further work on the neuron/glia stratified lane separate from the lane 2/3 manuscript surfaces.
