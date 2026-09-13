# v25 baseline review and development priorities

Reviewed 2026-09-13 against `CODEX_CONTEXT.md` and the three existing project files.
The launcher and dashboard engine are unchanged. This pass adds a synthetic regression
suite and this development plan; it does not implement Distributions or revalidate
the aerodynamic model against external CFD results.

## Baseline verification

Eight tests pass in `test_baseline.py` using Python, NumPy, pandas, and Node.js:

- Shallow POLAR discovery ignores `ADF_COMP` and does not read source contents.
- Keyed SREF/CREF parsing and authoritative case ordering.
- Six repeated Fluent headers produce three 22-row cases, with truncation isolated to the final case.
- A known linear CM/CL relationship gives 10% static margin.
- Drag rise uses the lowest Mach as its reference, including negative CLS filenames.
- Full synthetic generation, unchanged-cache reuse without invoking source parsers, timestamp invalidation, force rebuild, and no automatic CSV files.
- Python compilation and syntax checks on generated launcher/dashboard JavaScript.
- Actual generated JavaScript: the documented 0% MAC example, body moment shifts including lateral/vertical offsets, force invariance, duplicate-CL interpolation, no extrapolation, and ALPHA versus CL grids for Delta.

Run from the parent project folder:

```powershell
python -m unittest discover -s jamal_dashboard -p test_baseline.py -v
```

Python is not on PATH in the reviewed environment. The successful run used:

```powershell
& 'C:/Users/User/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -m unittest discover -s jamal_dashboard -p test_baseline.py -v
```

No representative CFD inputs, existing test suite, or Git repository were present.
These tests establish repeatable synthetic coverage; they do not cover every axis rotation,
all Delta metrics, browser rendering, interrupted network writes, or real-project performance.
No browser interaction or visual layout validation was performed in this pass.

## Recommended order

| Priority | Finding and consequence | Targeted improvement and acceptance check |
| --- | --- | --- |
| 1: scientific correctness | `parse_infout`, engine line 626: XREF/YREF/ZREF use a number regex that splits scientific notation. Reproduced: `XREF[m]: 1.1128E+01 YREF[m]: 0 ZREF[m]: 0` becomes `(1.1128, 1, 0)`, rather than `(11.128, 0, 0)`. This supplies wrong inputs to the protected transfer function. | Patch only keyed reference metadata parsing to accept exponents. Test decimal, signed, and exponent forms. Keep all transfer equations and validated derivative logic unchanged. |
| 1: data identity | `_normalize_configurations`, launcher line 315, accepts duplicate configuration labels. Case/history and curve keys use label + POLAR (+ case); separate configurations can overwrite histories or merge curves. Duplicate labels were accepted in a synthetic reproduction. | Reject blank-after-trimming and duplicate labels before generation, and show an inline error beside the relevant field. Test two different base folders with the same label/POLAR. Longer term, introduce stable configuration IDs while preserving readable labels. |
| 1: generation safety | `start_generation_job`, launcher line 548, starts a worker for every request. Two submissions were confirmed to start two workers. The lock protects job dictionaries, not generation. Clear-cache can also run while generation writes. Cache temporary filenames and output filenames are shared; stdout redirection is process-wide. | Serialize generation and cache clearing within the launcher and report an existing active job to the UI. If multiple launcher processes must be supported, add an output-directory lock. Test simultaneous requests, clear-cache during a job, and release after failure. |
| 1: report recovery | Launcher lines 508–533 writes final JSON, HTML, and manifest directly. Failure can leave an old HTML beside new JSON, or truncate an existing report. | Build both artifacts in unique temporary files in the output directory, validate them, then publish. Use a shared generation ID and recovery strategy: two file replacements alone are not a transaction. Inject writer failures and verify that the previous complete report remains usable. |
| 1: HTML integrity | `make_html`, engine line 1410, embeds JSON without escaping `<`; a configuration label containing `</script>` was confirmed to appear verbatim in the script. Several table/preset builders also interpolate text into `innerHTML`. Ordinary quotes and angle brackets can damage interface content; untrusted setup text can introduce markup. | Escape embedded JSON for HTML script context and use `textContent`/DOM option builders for labels and preset names. Test quotes, ampersands, Unicode, and a literal script terminator. Keep formatting markup separate from data. |
| 2: interaction speed | The active `refreshAll` implementation, engine line 3251, redraws plots and tables across all tabs. Plot builders repeatedly use `Plotly.newPlot`. Several plots recompute filtered/derived data. | Render the visible tab, mark other tabs stale, and refresh them when opened. Reuse plots with `Plotly.react`, preserve user zoom, and compute shared data once per refresh. Measure initial load and filter-change timings on a representative large report before claiming a speedup. |
| 2: progress recovery | `pollGenerationJob`, launcher line 798, polls via a 450 ms interval, allowing overlapping requests. A single fetch error discards the job ID and re-enables Generate even though the worker may still be running. Once log text exists, the current work message is not shown in the status text. | Poll sequentially after each response, retry transient failures with bounded backoff, retain/recover the job ID, and show the current POLAR separately from the log. Expose the active job through state so page reloads can reconnect. Test slow requests and temporary connection loss. |
| 2: portability | The generated HTML depends on a Plotly CDN script, engine line 1428. A single HTML file is generated, but plotting is not guaranteed offline. | Offer an explicit fully offline export that embeds a pinned Plotly bundle. Test with network access disabled and explain the increased report size. |

The line numbers above refer to the unchanged v25 source files reviewed here.
Reproduced defects are observations, not fixes; the regression suite currently tests
working baseline paths rather than asserting that these known defects are resolved.

## Additional improvements

- Keep the current curve styles and mixed typography. Add persistent keyboard focus,
  accessible tab semantics, and a live progress announcement when updating controls.
- For missing sources, show a preflight list with configuration, POLAR, missing file,
  and corrective action. Preserve strict generation as the default; an optional partial
  report should clearly identify omitted data rather than silently dropping it.
- Tie scan results to the path that produced them. After editing a base path, mark
  the selection as requiring a scan; ignore stale responses from earlier scans.
- Cache entries currently deserialize the full parsed payload during planning.
  Profile network cache I/O and peak memory before introducing a lightweight index
  or separate history loading. Preserve path + size + nanosecond mtime fingerprints.
- Recheck fingerprints after parsing before committing a cache entry when source
  files may still be written by the solver. Warn and retry or defer changed inputs.
- Expose elapsed time, current phase/POLAR, and actual parsing/reuse counts together.
  Treat progress as approximate until the work estimate accounts for changed files.
- Preserve current setup compatibility. Implement New/Open/Save/Save As and
  project-specific session storage incrementally; retain the existing JSON format.
- Establish version control or versioned release archives before production patches.
  Keep independent parser/module versions and add cache invalidation tests when a
  parser changes; `engine_version == v25` alone cannot distinguish edits within v25.

## Distributions: next feature after baseline safeguards

The context explicitly identifies DISTCLCP as the next major module. Implement it
as an additive parser/data module plus a new tab; do not restructure the existing
engine to make room for it.

1. Add the generic two-column curve parser, numeric station matching, geometry
   extraction, force parsing, and targeted integrity diagnostics with synthetic fixtures.
2. Discover components only beneath selected `DISTCLCP/POLAR-XXX` directories during
   generation. Keep Phase 1 limited to the existing ADF/DRAG-RISE discovery.
3. Add separately versioned geometry and per-state caches using existing fingerprint
   conventions; do not mix unvalidated distribution data into validated polar entries.
4. Map state order from `infout [CASES]`. The current `parse_infout` does not expose
   `qdin`, rho, or velocity. Confirm the available dynamic-pressure values and their
   association with states before computing `cl = L'/(q*c)`.
5. Verify the leading-edge direction from representative geometry and Cp files.
   Do not permanently label Xmin as the leading edge without that evidence. Until
   verified, retain raw coordinates or explicitly mark normalization as unvalidated.
6. Add configuration/POLAR/component/state selection, overlays, station selection,
   inverted Cp axis, informative empty states, and consistent trace styles. Render
   distributions only when their tab is opened.
7. Validate against an independently checked state before removing the module's
   unvalidated label. Defer Delta-cl/Delta-Cp and component-specific VTAIL assumptions
   until their geometry and conventions are established.

The input needed for physical validation is a representative selected POLAR's
`infout`, `section_state1_station*`, `total_force_state*`, and `cp_dist_state*_station*`
files, with confirmation of force units/frame and leading-edge direction.

## Source baseline fingerprints

SHA-256, recorded after the read-only source review:

```text
jamal_dashboard_launcher_v25.py
0352522ABD5E6AD16DEE4FC6C0A829E053933522A2C78D292255290C935C31B0

jamal_polar_convergence_dashboard_v25.py
A657295E2AE10E4D2E6127353777D70CF26B9B4CE18AAD344DADEAC37FA2E43B
```
