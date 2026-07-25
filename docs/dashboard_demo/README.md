# IlluMeta interactive dashboard QA fixture

Open `Case_vs_Control_results_index.html` in a browser to exercise the offline
IlluMeta dashboard. The values are synthetic interface-test data, not study
results and not a reconstruction of any cohort.

Regenerate the artifact from the repository root:

```bash
CONDA_BIN="$(command -v conda || printf '%s' "$HOME/miniforge3/bin/conda")"
"$CONDA_BIN" run -n illumeta python scripts/build_dashboard_demo.py
```

The fixture includes a keyboard/hover-responsive volcano plot plus searchable,
sortable, CSV-exportable DMP tables. `manifest.json` records deterministic hashes
for the dashboard and every synthetic input.
