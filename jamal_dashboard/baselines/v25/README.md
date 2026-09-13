# JAMAL Aerodynamic Results Dashboard v25

This package contains a lightweight local launcher and the standalone dashboard engine.
No Flask, Dash, or web framework is required.

## Files

- `jamal_dashboard_launcher_v25.py` — local browser launcher, shallow folder discovery, incremental cache, live progress, setup save/load.
- `jamal_polar_convergence_dashboard_v25.py` — validated standalone dashboard engine.

Keep both files in the same folder.

## Run

```bash
python jamal_dashboard_launcher_v25.py
```

The launcher opens in your default browser. If the browser does not open, use the local URL printed in the terminal.

## Expected JAMAL structure

```text
<BASE>/
├── 00-SUPPORT/
├── 01-GRIDS/
├── 02-RUNS/
│   └── POLAR-XXX/
│       ├── infout
│       └── FLUENT_LOG
└── 03-RESULTS/
    ├── ADF/
    │   ├── POLAR-XXX.adf
    │   └── ADF_COMP/        # ignored during discovery
    ├── DRAG-RISE/
    │   └── <optional folders>/
    └── DASHBOARD/           # created automatically
        ├── dashboard.html
        ├── dashboard.json
        └── .jamal_cache/
```

## Two-phase workflow

### Phase 1 — Fast discovery

For each configuration, select the JAMAL base directory. The launcher performs only:

- one shallow listing of `03-RESULTS/ADF`;
- one shallow listing of `03-RESULTS/DRAG-RISE`.

It does not recurse into `ADF_COMP` and does not inspect `02-RUNS` during discovery.

### Phase 2 — Generate / Update dashboard

Only selected POLARs are inspected. The launcher fingerprints the selected:

- `POLAR-XXX.adf`;
- `infout`;
- `FLUENT_LOG`.

The fingerprint uses source path, file size, and modification time.

- New POLAR: parsed and cached.
- Modified POLAR: reprocessed and cache replaced.
- Unchanged POLAR: reused from cache.
- Unchecked POLAR: omitted from the new dashboard but may remain cached for later reuse.

The standalone HTML and JSON are regenerated from the complete current selection, which is inexpensive compared with parsing large Fluent logs.

## Generation status

The top panel shows:

- generation phase;
- percentage progress;
- new, modified, cached, and forced-rebuild counts;
- current configuration/POLAR;
- parsing log;
- final parsed/reused counts and elapsed time.

## Cache controls

- **Generate / Update dashboard** — parse only new or modified sources.
- **Force full rebuild** — ignore cache for the current run.
- **Clear cache** — remove `03-RESULTS/DASHBOARD/.jamal_cache`; generated dashboard files are kept.

## Outputs

The first configuration is the output project. Files are written to:

```text
<first configuration>/03-RESULTS/DASHBOARD/dashboard.html
<first configuration>/03-RESULTS/DASHBOARD/dashboard.json
```

Automatic `convergence_comparison.csv` and static-margin CSV files are no longer generated. CSV exports remain available explicitly inside the standalone dashboard.

## Setup files

**Save setup** stores configuration labels, base directories, selected POLARs, and drag-rise directory selections in a small JSON file. It does not store CFD results.

**Load setup** restores those selections and performs fast discovery for the configured base folders.

## Dependencies

The dashboard engine uses:

- Python 3.9+
- NumPy
- pandas

Plotly is loaded by the generated standalone HTML from the configured CDN reference, as in previous versions.
