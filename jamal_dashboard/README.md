# JAMAL Aerodynamic Results Dashboard v25.4

Coefficients shows three force plots in a row, then three moment plots below,
followed by L/D. Choose Body, Stability or Wind axes, then ALPHA, BETA,
selected-axis CL or CY as the horizontal coordinate. Self-plots are omitted.
Headings show Reynolds and Mach, plus Xref for moments. The moment reference
applies to all moment plots, and L/D uses CL/CD from the chosen axes.
Delta also provides independent axis and coordinate selectors. Drag rise supports
all three axes using their source drag columns and a lowest-Mach baseline;
file groups retain their original CLS targets. Missing axis data is reported.
These controls are retained in saved views and presets. Numeric ticks and values
use an aligned monospace font; word headings retain their existing typography.
Spanwise cl/cl.c, static margin, drag rise and convergence also use two columns,
stacking on screens below 760 px. The Cp station layout is unchanged.

Each Cp station panel includes its section outline underneath, using a shared
x/c axis and an independent y/c axis with equal geometric scaling. The ordinate
comes directly from the section file divided by local chord (no sign change,
rotation or recentering). Xmin/Xmax selection applies to both Cp and geometry;
dimensional-X mode uses metres. Matching colors, dashes and legend toggles link
the pressure curve and outline. Section geometry comes from section_state1.

Cp hover labels use transparent backgrounds and curve-colored text to keep small
station plots readable. Local Git workflow instructions are in ../GIT_WORKFLOW.md;
the v25.3 tag preserves the version before this hover styling change.

## New in v25.3

All span coordinates use infout BREF, including VTAIL. The span-coordinate selector
also controls Cp station labels: 2Y/BREF and Y/BREF appear as percentages; Y [m]
appears in metres. Component span metadata remains geometry provenance only.

Separate spanwise plots show cl and cl.c (cl times local chord, in metres).
Legends display alpha and beta. Save spanwise loads (.txt) exports every station
of the displayed states and pinned overlays, with both loads, chord, BREF and
coordinate metadata. Missing values are marked NA; Cp station selection does
not limit the load export.

Cp uses a responsive grid of separate station plots with common inverted limits.
Selections are retained when switching coordinate modes. Different components
remain in separate panels. All stations are selected initially.

The synthetic dataset now includes a paired **VTAIL with a 5 m projected span and
45 degree cant**, alongside the original 10 m wing. See
`JAMAL_SYNTHETIC_CFD/README_VTAIL.md` for geometry, force-density conventions and
expected values. The original source fixtures are unchanged. Production V-tail
physics remains unvalidated; the new data is an illustrative test fixture.

Current development suite: **26 passing tests**, including numerical tail values,
span-reference changes, normalized station grouping and original protected modules.

This package contains a lightweight local launcher and the standalone dashboard engine.
No Flask, Dash, or web framework is required.

## Files

- `jamal_dashboard_launcher_v25.py` — local browser launcher, shallow folder discovery, incremental cache, live progress, setup save/load.
- `jamal_polar_convergence_dashboard_v25.py` — validated standalone dashboard engine.

Keep the two scripts and these new assets in the same folder:

- `jamal_distributions.py` — separate DISTCLCP parser, station matching, raw-file cache and integrity checks.
- `distributions.html` and `distributions.js` — embedded Distributions tab; no separate files are needed beside a generated report.

Script filenames remain `*_v25.py` for launcher compatibility. Version v25.3 is
reported inside the application. The untouched original scripts and README are
preserved in `baselines/v25`.

## New in v25.1

- **Distributions:** sectional cl and inverted-axis Cp plots; state, POLAR,
  component and configuration selection; overlays; station selection, every-nth
  selection and single-station mode; remembered selections per source project.
- Separate incremental geometry/force/Cp caches. Only selected POLARs are inspected
  during generation. Original shallow discovery remains unchanged.
- Station matching within 1e-6 m, original Cp surface order, dynamic component
  names, source provenance and explicit integrity warnings.
- Exponent-aware reference coordinates and qdin parsing. Duplicate configuration
  labels and overlapping generation/cache-clear operations are rejected.
- Reports are staged and validated before publication, with rollback if ordinary
  file replacement fails. Progress monitoring retries transient errors and can
  reconnect to active jobs after a page reload.
- Existing initialized plots use Plotly.react for updates and retain zoom when
  axis definitions remain the same. New distribution plots render when their tab opens.

The moment-transfer equations, static-margin derivative, Fluent subblock parser,
drag-rise definition, Delta interpolation, and residual thresholds remain unchanged.

## Distributions conventions

Inputs: `03-RESULTS/DISTCLCP/POLAR-XXX/<component>/` containing
`section_state1_station*`, `total_force_stateS`, and `cp_dist_stateS_station*`.
States follow the original infout case order. Force data is N/m, using the body
frame described in `CODEX_CONTEXT.md`; cl uses `L'/(qdin*chord)` with no delta-Y factor.

The new module has been checked against the supplied synthetic dataset. Production
validation is pending. Nonzero-beta cases, missing/non-positive qdin, and changing
flow conditions without per-state qdin produce warnings and omit sectional cl.
Cp remains available when valid geometry is present.

The default Cp normalization assumes Xmin is the leading edge, as explicitly
confirmed by the synthetic fixture. For other geometry, select Xmax or dimensional
X until its orientation is confirmed. The UI never sorts Cp points by X, so upper
and lower surface branches remain connected in file order. Use dimensional Y for
components where the aircraft BREF is not an appropriate span normalization.

## Validation

```bash
python -m unittest discover -s jamal_dashboard -p "test_*.py" -v
```

Run from the parent folder; from inside `jamal_dashboard`, use `-s .` instead.
Tests require NumPy, pandas and Node.js, and use disposable copies of
`JAMAL_SYNTHETIC_CFD`. They cover numerical fixture values, cache reuse, failure
handling and unchanged protected source functions. See `VALIDATION_v25_1.md`.

The launcher is intended to run as a single process per output project. Its job
lock does not coordinate separate launcher processes. Staged publication handles
ordinary writer failures, but two filesystem paths are not a crash-atomic
transaction and external readers are not covered by the HTTP read lock.

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
