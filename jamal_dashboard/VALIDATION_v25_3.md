# v25.3 validation - 2026-09-13

All 26 regression tests pass, including protected scientific calculations,
parser behavior, cache reuse, original fixture hashes and generated JavaScript syntax.
Actual JavaScript tests cover BREF-only normalization, all three station-label
modes, component grouping, cl.c including zero/missing values, and text-export rows.

Browser checks verified nine VTAIL station panels across coordinate changes,
two-state overlays (18 Cp curves), separate cl and cl.c plots, alpha/beta legends,
and no console errors. Tail tip labels use 50% (2Y/BREF), 25% (Y/BREF), or 2.500 m.
The export button executed and reported 18 rows for two displayed states.
The in-app browser did not expose a download completion event; native save
completion remains to be checked in the target browser. Export contents are
covered by the actual JavaScript tests.

Validated aerodynamic calculations and all force/geometry/Cp inputs are unchanged.
Component metadata remains provenance; the UI always uses infout BREF.
Existing beta=0, leading-edge and synthetic-only validation limits remain.
