# v25.2 validation — 2026-09-13

- All **26 project tests pass**, including the original protected-module regressions.
- The supplied original fixture validator still passes **1,654 checks**.
- Original SHA256SUMS.txt hashes remain unchanged; VTAIL inputs have a separate hash manifest.
- Combined fixture: 20 distribution state series (10 WING + 10 VTAIL), 236 raw
  geometry/force/Cp files, and two component-reference JSON files.
- Every VTAIL cl, chord and Cp minimum agrees with expected_vtail_values.json.
  Independent trapezoidal Cp integration matches the expected normal loading;
  left/right side forces cancel. State 5 retains alpha=-1 degree.
- Tail projected span is 5 m while the shared aircraft BREF remains 10 m.
  Metadata-only reference changes leave unchanged raw distribution caches reusable.
- Actual JavaScript helper tests check full-span percentages, semispan percentages,
  grouping equivalent normalized stations across references, separation of components,
  and dimensional fallback for invalid/missing span references.
- Browser checks show a responsive two-column grid at the inspected narrow desktop
  viewport, one panel per station. Nine stations with two state overlays produce
  nine plots and eighteen curves. All plots share inverted Cp limits.
- Full-span/semispan selection, single-station mode, clearing stations and every-nth
  selection update panel counts correctly. No browser errors were captured in those checks.

The paired V-tail is an illustrative 45-degree-cant synthetic fixture, not production
validation. Its force densities use projected Y metres. The unchanged whole-aircraft
ADF fixture does not include the new tail loads. See README_VTAIL.md in the dataset
for the construction and coordinate assumptions.

The operational limits documented for v25.1 remain: beta=0 sectional lift, verified
dynamic-pressure mapping, confirmed leading-edge direction, one launcher per output
project, and the existing Plotly CDN dependency.
