# v25.1 validation — 2026-09-13

## Automated checks

- Supplied fixture validator: **1,654 checks passed**.
- All **136 packaged file hashes** match the supplied SHA256SUMS.txt; fixture inputs are unchanged.
- Project regression suite: **22 tests passed**.
- Python compilation and generated JavaScript syntax checks passed.
- Synthetic generation: two POLARs, ten 240-iteration cases, two drag-rise curves,
  ten distribution state series, 90 Cp curves and 118 separately cached input files.
- All 90 sectional cl values and geometry chords match `expected_values.json`.
  Every Cp curve retains 81 points and the correct surface order; upper Cp minima
  and span-integrated lift match the supplied expectations.
- Both static-margin curves remain at 10%; saved-v25 comparisons confirm protected
  parser, moment, derivative, drag-rise and Delta source functions are unchanged.
- Unchanged generation reused both POLARs and all 118 distribution files. Changing
  one Cp source reparsed only that file; forced rebuild reparsed all distribution inputs.
- Missing geometry, missing force states, unmapped states, invalid chords, malformed
  curves, missing dynamic pressure, and nonzero beta exercise diagnostic paths.
- Multi-configuration data retains separate series while sharing unchanged raw
  distribution cache entries.
- Concurrent generation and cache clear are rejected. Injected writer failure
  preserves previous reports; injected failure on the second report replacement
  restores the first report file.

## Browser checks

Tested through the actual local launcher and report in the Codex in-app browser:

- Shallow scan exposes exactly POLAR-001 and POLAR-002 and BASELINE drag rise.
- Generation reaches 100%, with output paths and per-module reuse counts displayed.
- Distributions opens with nine station controls, a cl plot with root cl=0.30, and
  an inverted Cp axis. It reports no distribution integrity warnings for the fixture.
- Overlaying states 1 and 5 produces two load traces and 18 Cp curves when all
  stations are selected. State 5 retains alpha=-1 degrees.
- Reload restores distribution selections, overlays and station choices.
- Switching POLAR, using single-station mode and clearing stations updates the
  visible curve count correctly. Empty station selection displays guidance.
- Legends were inspected visually at a narrow viewport and adjusted to avoid
  overlap with axis titles. Full trace identity remains available in hover text.
- Existing Body/Wind coefficient titles update correctly, static margin still
  renders at 10%, and no JavaScript errors were captured in the tested workflows.

## Limits and remaining work

These are synthetic checks, not production CFD certification. Nonzero-beta
sectional lift and component-specific VTAIL conventions remain unvalidated.
Distribution Delta-cl/Delta-Cp is deferred. Cp leading-edge direction must be
confirmed for each new geometry family.

The report still loads Plotly from the existing pinned CDN. Fully offline Plotly
embedding, broader legacy-table HTML escaping, cancellation, and complete lazy
rendering of all legacy tabs remain future work. Progress-retry and active-job
reconnection are implemented but long network outages were not fault-injected
through the browser. No large network-drive performance benchmark was run.

The original `DEVELOPMENT_REVIEW.md` remains a historical baseline audit; its
known-issue list predates this implementation. The current README describes what
has now been implemented and the remaining operational limits.
