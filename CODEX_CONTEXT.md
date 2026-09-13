# Current update: v25.3 (2026-09-13)

Post-v25.3 interface update: Cp hover labels have a transparent background and
border, with text matching the curve color. Numerical calculations are unchanged.
Local Git history now tracks this project, with v25.3 marking the validated
baseline. See GIT_WORKFLOW.md for commands and the Windows ownership note.

User correction supersedes earlier component-normalization descriptions: all plot
coordinates use infout BREF, including VTAIL. Physical tail projected span stays
5 m; BREF stays 10 m. Cp labels follow the selected 2Y/BREF, Y/BREF or Y [m]
coordinate. Normalized station labels display percentages.
Separate cl and cl.c plots use alpha/beta legends. The text-export button includes
both loads for all displayed series/stations, with chord and reference metadata.
Validated numerical modules and original synthetic inputs remain unchanged.

# CODEX_CONTEXT.md

## Latest update — v25.2, 2026-09-13

Cp is now displayed in a responsive grid, one panel per selected station. Default
labels use 100*Y/component span reference; optional semispan labels use 200*Y/reference.
Overlays share panels by normalized station within each component. Components are
kept separate, Cp axes share limits, and all stations are selected initially.

Distributions v0.2 reads optional `component_reference.json` with `span_reference_m`;
otherwise it uses infout BREF. The user confirmed a **paired V-tail**, with 5 m
projected span (half the 10 m wing). The new additive synthetic VTAIL fixture uses
an assumed 45-degree cant, body forces per projected-Y metre, and a separate tail
reference. Existing WING, infout, logs, ADF and drag-rise files are unchanged.
Tail loads are not added to the old whole-aircraft ADF. See the new
`jamal_dashboard/JAMAL_SYNTHETIC_CFD/README_VTAIL.md` for the exact assumptions.

The current suite has 26 passing tests. Protected v25 calculations remain unchanged;
production V-tail validation is still pending. The v25.1 notes below are historical.

## Implementation update — 2026-09-13

The working package is now **v25.1**, retaining the `*_v25.py` launcher/engine
filenames for compatibility. The untouched original v25 package is saved under
`jamal_dashboard/baselines/v25`.

Distributions v0.1 has been implemented as the additive `jamal_distributions.py`
module with `distributions.html` and `distributions.js`, embedded during report
generation. It supports state/configuration/POLAR overlays, geometry/force/Cp
caching, station matching, cl and Cp plots, integrity diagnostics and source provenance.
It is verified against `JAMAL_SYNTHETIC_CFD`; production validation remains pending.
The fixture confirms Xmin is its leading edge, qdin is dynamic pressure, forces
are N/m, and numeric station matching needs 1e-6 m tolerance. The UI offers Xmax
and dimensional-X alternatives for other geometry. Sectional cl currently requires
beta=0, positive qdin and constant flow within a POLAR.

Targeted safeguards fix exponent reference parsing, duplicate labels, concurrent
jobs/cache clearing, staged report publication, and progress retries. Existing
plots reuse Plotly.react; distribution rendering is deferred until its tab opens.
Protected moment transfer, static-margin derivative, Fluent subblock parsing,
drag-rise and Delta calculations and residual thresholds remain unchanged.

Read `jamal_dashboard/README.md` and `jamal_dashboard/VALIDATION_v25_1.md` for current
capabilities and limits. Run `python -m unittest discover -s jamal_dashboard -p "test_*.py" -v`
from this folder. The specifications below remain the historical
requirements; sections describing Distributions as entirely unimplemented are
superseded by this update. Distribution deltas remain future work.

## Purpose

This file is the development context for the JAMAL CFD post-processing dashboard.

Read this file **before modifying the codebase**.

The current project is not a prototype anymore. It contains several validated and stable behaviors that must not be casually rewritten. Future development should preserve working functionality unless a change is explicitly required.

The dashboard has evolved from a simple aerodynamic polar plotter into a broader CFD post-processing application covering:

- aerodynamic coefficients;
- convergence assessment;
- numerical quality;
- aerodynamic final-window stability;
- neighbor/outlier consistency;
- mesh quality;
- drag-rise analysis;
- static margin;
- moment-reference transformation;
- multi-configuration comparison;
- incremental loading/cache;
- portable standalone HTML reporting;
- planned spanwise load and chordwise Cp distributions.

---

# 1. Current Baseline and Version History

## Current launcher/package baseline

**v25** is the current working launcher/package baseline.

It includes:

- local browser-based launcher;
- base-folder selection;
- shallow ADF and DRAG-RISE discovery;
- POLAR checkbox selection;
- incremental parsing;
- persistent cache;
- live generation progress;
- output under `03-RESULTS/DASHBOARD`;
- no automatic CSV generation.

## Current dashboard behavior baseline

The aerodynamic dashboard evolved through v15-v24.

### v15
Moment-reference transformation was validated against an independent source.

Treat the validated moment-transfer behavior as a protected module.

### v17+
Added:
- Convergence Classification 2.0;
- mesh-quality integration;
- automatic outlier detection;
- module versioning;
- navigation tabs.

### v18
Mixed typography:
- Claude-inspired serif style for words/headings/labels;
- original sans-serif style retained for numerical values, ticks, legends, hover values, table numbers.

### v20
Configuration-difference logic:
- `ΔCL` vs `ALPHA`;
- `ΔCD`, `ΔCM`, `ΔL/D`, `ΔStatic Margin` vs `CLS`;
- comparison curve interpolated onto the reference curve’s CL grid.

### v21
Dedicated tabs:
- Coefficients;
- Static Margin;
- Delta;
- Drag Rise;
- Convergence;
- Quality;
- Tables.

Also:
- Xref title formatting;
- Mach/Re formatting;
- multiple comparison pairs.

### v24
Curve-style control:
- Auto;
- Colors only;
- Color by POLAR;
- Color by configuration;
- Custom.

### v25
Incremental cache + progress + new output location.

---

# 2. JAMAL Project Folder Structure

The user now supplies the JAMAL **base folder**, not the ADF path.

Expected base layout:

```text
BASE_FOLDER/
├── 00-SUPPORT/
├── 01-GRIDS/
├── 02-RUNS/
└── 03-RESULTS/
    ├── ADF/
    │   ├── POLAR-001.adf
    │   ├── POLAR-002.adf
    │   ├── ...
    │   └── ADF_COMP/
    │       └── many component ADF files
    ├── CONVERGENCE/
    ├── DISTCLCP/
    ├── DRAG-RISE/
    ├── FLOWVIS/
    ├── PROBE/
    └── DASHBOARD/
```

## Important discovery rule

The launcher must perform a **shallow scan only** of:

```text
03-RESULTS/ADF
03-RESULTS/DRAG-RISE
```

Do **not** recursively traverse `ADF`.

`ADF_COMP` must be ignored during POLAR discovery.

Available POLARs are only root-level files matching:

```text
POLAR-XXX.adf
```

---

# 3. Two-Phase Launcher Workflow

## Phase 1 — Fast Discovery

Initial scan must be cheap, especially because user paths may be on a mapped network drive such as `Z:`.

Only do:

1. shallow list of `03-RESULTS/ADF`;
2. collect root `POLAR-XXX.adf`;
3. shallow list of first-level directories under `03-RESULTS/DRAG-RISE`.

Do **not** inspect during Phase 1:

- `02-RUNS`;
- `infout`;
- `FLUENT_LOG`;
- ADF file contents;
- drag-rise file contents.

The POLAR checkbox list should remain compact and should **not** display repeated text like:

```text
ADF discovered · run files deferred
```

## Phase 2 — Generate / Update Dashboard

Only selected POLARs are processed.

For each selected POLAR:

```text
02-RUNS/POLAR-XXX/infout
02-RUNS/POLAR-XXX/FLUENT_LOG
03-RESULTS/ADF/POLAR-XXX.adf
```

Only selected drag-rise directories are processed.

---

# 4. Incremental Cache Architecture

Persistent cache location:

```text
03-RESULTS/DASHBOARD/.jamal_cache/
```

The output dashboard belongs under the **first configuration's base folder** by default.

## Cache fingerprint

Use at least:

- full path;
- file size;
- modification timestamp.

Do not hash large Fluent logs by default because network hashing can eliminate performance gains.

## Cache behavior

Each selected POLAR is classified as:

- NEW;
- MODIFIED;
- CACHED / UNCHANGED;
- REMOVED from current selection;
- FORCED REBUILD.

### Normal Generate / Update

Reuse unchanged POLAR parsed data.

Parse only:
- new POLARs;
- modified POLARs;
- new/modified drag-rise files.

### Force full rebuild

Ignore valid cache entries for the current generation.

Re-read every selected source file and then refresh cached entries.

It does **not** delete the cache first.

### Clear cache

Delete:

```text
03-RESULTS/DASHBOARD/.jamal_cache/
```

Do not delete:
- dashboard HTML;
- dashboard JSON;
- ADF;
- infout;
- FLUENT_LOG;
- drag-rise source files.

The next generation behaves as a fresh parse.

---

# 5. Generation Progress

The launcher must show live progress at the **top of the page**.

Desired display:

```text
Generating dashboard                                  63%

███████████████████████████░░░░░░░░░░

Parsing Config_2 · POLAR-018
New: 2 · Modified: 1 · Cached: 15
```

Suggested phase weighting:

```text
5%   Validate configuration and paths
10%  Determine cache/update plan
55%  Parse new and modified POLARs
10%  Process drag-rise files
10%  Build classifications and derived data
7%   Generate dashboard HTML/JSON
3%   Finalize and verify outputs
```

Progress should be based on actual work required in an incremental update.

---

# 6. Dashboard Output

Output directory:

```text
03-RESULTS/DASHBOARD/
```

Create it if missing.

Generate:

```text
dashboard.html
dashboard.json
```

Do **not** automatically generate:

```text
convergence_comparison.csv
convergence_comparison_static_margin.csv
```

CSV should only be generated through explicit dashboard export actions.

The dashboard header should be:

```text
JAMAL · CFD post-processing
Aerodynamic results dashboard
Interactive aerodynamic, convergence, drag-rise, and reference-moment analysis in a portable standalone report.
```

---

# 7. Main Dashboard Tabs

Current UI organization should remain close to:

- Overview
- Coefficients
- Static Margin
- Delta
- Drag Rise
- Convergence
- Quality
- Tables

A new tab is planned:

- Distributions

---

# 8. Typography and Visual Style

The user prefers the general v15/v18 visual style.

Use:

- Claude-inspired serif style for words/headings/labels/buttons;
- original sans-serif style for numbers.

Keep numerical fonts unchanged for:

- plot ticks;
- legend numerical content;
- hover numerical values;
- input values;
- table numerical values.

Avoid the full Claude-style color theme used in v16; it did not fit well with the plots.

---

# 9. Curve Style System

Supported modes:

## Auto

Default.

### One configuration, multiple POLARs

- different color per POLAR;
- all solid lines.

### Multiple configurations with matching POLARs

- same color for the same POLAR;
- different line styles for configurations.

Example:

```text
POLAR-014:
Baseline = solid blue
Config_2 = dashed blue

POLAR-015:
Baseline = solid orange
Config_2 = dashed orange
```

## Colors only

- unique colors;
- all solid lines.

## Color by POLAR

- color maps to POLAR;
- line style can map to configuration.

## Color by configuration

- color maps to configuration;
- line style/marker differentiates POLARs.

## Custom

Allow independent mapping of:

- color;
- line style;
- marker.

The style selection must remain consistent across dashboard plots.

---

# 10. ADF Data

ADF stands for **Aerodynamic Data File**.

The ADF contains coefficient columns such as:

```text
Body:
CDB, CYB, CLB, CRB25, CMB25, CNB25

Stability:
CDS, CYS, CLS, CRS25, CMS25, CNS25

Wind:
CDW, CYW, CLW, CRW25, CMW25, CNW25
```

The dashboard Coefficients tab has:

## Axis selector

- Body (B)
- Stability (S)
- Wind (W)

## Coefficient buttons

- CD
- CY
- CL
- CR
- CM
- CN

Dynamic mapping:

```text
Body:
CD -> CDB
CY -> CYB
CL -> CLB
CR -> CRB
CM -> CMB
CN -> CNB

Stability:
CD -> CDS
CY -> CYS
CL -> CLS
CR -> CRS
CM -> CMS
CN -> CNS

Wind:
CD -> CDW
CY -> CYW
CL -> CLW
CR -> CRW
CM -> CMW
CN -> CNW
```

The plot title must reflect the selected axis system and coefficient.

Do not use misleading fixed title text such as:

```text
Stability coefficient vs ALPHA
```

when Body or Wind axis is selected.

---

# 11. Aerodynamic Plot Titles

Do not put Xref, Mach, or Reynolds in trace names.

They belong in titles only.

Formatting:

```text
Xref = integer % MAC
Mach = 2 decimals
Re = scientific notation with no decimal mantissa
```

Examples:

```text
Xref = 25% MAC
M = 0.20
Re = 2E+07
```

For moments, do not expose the historical `25` suffix in plot labels.

Use:

```text
CMB
CMS
CMW
```

instead of:

```text
CMB25
CMS25
CMW25
```

Static-margin titles should also show the active X reference in %MAC.

---

# 12. Moment Reference Transfer — VALIDATED MODULE

This is one of the most important protected behaviors.

The validated implementation originated in v15 and was later exposed across Body/Stability/Wind views.

Do **not** rewrite casually.

## infout reference frame

`XREF`, `YREF`, `ZREF` in `infout` use a different coordinate frame:

```text
X_info positive backward
Y_info same lateral direction
Z_info positive upward
```

## aerodynamic/body calculation frame

```text
X_body positive forward
Y_body same
Z_body positive downward
```

Therefore displacement conversion is:

```text
dx_body = -dx_info
dy_body =  dy_info
dz_body = -dz_info
```

## Original reference assumption

The user assumes:

```text
XREF = 25% MAC
CREF = MAC
```

Correct %MAC mapping:

```text
X_LE_MAC = XREF - 0.25*CREF

XNEW = XREF + ((MAC_PERCENT - 25.0)/100.0)*CREF
```

For example:

```text
CREF = 4.629
XREF = 11.128

0% MAC:
XNEW = 9.97075
dx_info = -1.15725
```

## Body force signs

ADF reports intuitive positive drag/lift:

```text
CDB > 0 drag
CLB > 0 lift
```

Physical body-axis force coefficients are:

```text
CX = -CDB
CY =  CYB
CZ = -CLB
```

## Moment transfer

Apply in body axes:

```text
M_new = M_old - r × F
```

Coefficient form:

```text
CRB_new = CRB_old - (dy*CZ - dz*CY)/BREF

CMB_new = CMB_old - (dz*CX - dx*CZ)/CREF

CNB_new = CNB_old - (dx*CY - dy*CX)/BREF
```

Then rotate corrected body moments into Stability and Wind axes using the established validated rotation implementation.

For beta = 0, body/stability Y axes coincide.

## UI

Moment reference selection includes:

```text
Original
User
```

The active reference selection should affect moment plots consistently for:

- Body;
- Stability;
- Wind.

Forces are reference-point independent.

## Validation

The user explicitly validated the resulting moment values against an independent source.

Treat the module as:

```text
Moment Transfer Module v1.0 — VALIDATED
```

---

# 13. Static Margin

Definition:

```text
SM[%] = -(dCM/dCL)*100
```

Historically the dashboard used:

```text
CMS vs CLS
```

with central-difference / second-order behavior.

The stable branch should be preserved.

Do not casually rewrite the derivative logic.

Static margin is displayed in a dedicated tab.

Moment-reference changes must propagate consistently to the static-margin calculation only through the validated path.

---

# 14. Delta / Configuration Difference Tab

The tab name is:

```text
Delta
```

not `Configuration difference`.

Supports multiple comparison pairs.

For each comparison pair:

- Reference curve;
- Comparison curve.

Rules:

## ΔCL

Plot:

```text
ΔCL vs ALPHA
```

because CL difference is naturally compared at the same angle.

## Other aerodynamic deltas

Use the reference curve's CL grid.

Plot:

```text
ΔCD vs CLS
ΔCM vs CLS
ΔL/D vs CLS
ΔStatic Margin vs CLS
```

The comparison curve must be linearly interpolated onto the reference curve’s CL values.

Duplicate CL values must be consolidated safely before interpolation.

This logic was validated operationally in v20+.

---

# 15. Fluent Convergence Parser

Log filename:

```text
FLUENT_LOG
```

Typical iteration header:

```text
iter  continuity  x-velocity  y-velocity  z-velocity
      energy  nut  tstep-ave  cp-max
      cnzb  cmyb  crxb  clzb  cyyb  cdxb  time/iter
```

Meaning:

```text
continuity, x-velocity, y-velocity, z-velocity, energy, nut
    = equation residuals

tstep-ave
    = average time step

cp-max
    = maximum Cp in the fluid

cnzb, cmyb, crxb, clzb, cyyb, cdxb
    = body-axis aerodynamic coefficients

time/iter
    = time per iteration
```

## Critical Fluent subblock behavior

Fluent prints the monitor table in repeated screen subblocks, typically 11 iteration rows followed by a repeated header.

Do **not** interpret every repeated header as a new aerodynamic case.

The parser must concatenate repeated subblocks.

The actual ALPHA/BETA cases are mapped using the `[CASES]` table from `infout`.

This bug was previously observed as:

```text
actual_iters = 11
```

for every case.

That is incorrect.

The fixed parser is a protected behavior.

---

# 16. infout File

The corresponding `infout` lives at:

```text
02-RUNS/POLAR-XXX/infout
```

Important sections include:

```text
[REFERENCE_DATA]
[FLOW_CONDITION]
[RUDDER]
[ELEVON]
[AILERON]
[FLAP]
[CASES]
```

Typical reference values:

```text
SREF[m2]
CREF[m]
BREF[m]

XREF[m]
YREF[m]
ZREF[m]
```

Typical flow values:

```text
Mach
Reynolds
qdin
rho
V
```

The `[CASES]` table maps state sequence:

```text
CASE MACH REYNOLDS ALPHA BETA NAME ITERS
```

The state sequence is authoritative.

For example:

```text
CASE 0001 -> state1
CASE 0002 -> state2
CASE 0003 -> state3
```

This is important for the planned DISTCLCP module.

---

# 17. Residual Thresholds

Current desired thresholds:

```text
RESIDUAL_SUCCESSFUL = 8.0e-4
RESIDUAL_ACCEPTABLE = 8.0e-3
RESIDUAL_BAD        = 1.0e-2
```

Interpretation:

```text
< 8e-4           successful
8e-4 to 8e-3     acceptable
8e-3 to 1e-2     suspicious
>= 1e-2          bad/diverged
```

Early residual or cp-max spikes can be acceptable.

Do not mark a case diverged solely because:

```text
cp-max spike ratio > some threshold
```

Final behavior is more important.

---

# 18. Convergence Classification 2.0

Do not collapse all diagnostics into a single black-box label.

Maintain separate fields:

```text
Numerical convergence
Aerodynamic convergence
Neighbor consistency
Mesh quality
Overall
```

Example:

```text
Numerical:      ACCEPTABLE
Aerodynamic:    WARNING
Neighbor check: PASS
Mesh:           PASS
Overall:        SUSPICIOUS
```

Reasons should be visible.

---

# 19. Final-Window Aerodynamic Statistics

For final-window analysis, compute per relevant quantity:

```text
mean
minimum
maximum
range
standard deviation
linear drift per 100 iterations
```

for at least:

```text
CL
CD
CM
cp-max
```

The dashboard previously used roughly the final 200 iterations.

The final-window behavior should remain inspectable and configurable.

A message such as:

```text
CLZB final-window range
```

means:

```text
max(CLZB in final window) - min(CLZB in final window)
```

This is aerodynamic stability, not an equation residual.

---

# 20. Residual Plots

Residual axes should use explicit scientific notation:

```text
1e-2
1e-3
1e-4
1e-5
```

Do not use SI-prefix notation such as `μ`.

Residual selection buttons should include:

- continuity;
- x-velocity;
- y-velocity;
- z-velocity;
- energy;
- nut;
- All.

---

# 21. Mesh Quality Integration

Extract when present from Fluent logs:

- cell count;
- node count;
- minimum orthogonal quality;
- maximum aspect ratio;
- negative volume warnings;
- mesh-related warnings.

Expose mesh status independently in the Quality tab.

---

# 22. Automatic Outlier Detection

Use local neighbor consistency for aerodynamic data.

Primary quantities:

```text
CLS
CDS
CMS
```

Use neighboring sweep points and local interpolation.

Flag suspicious points rather than automatically declaring them unusable.

The dashboard should keep numerical convergence and aerodynamic/neighbor consistency as separate concepts.

---

# 23. Control-Surface Metadata

Parse deflections from `infout`:

```text
RUD1-RUD4
ELV1-ELV4
AIL1-AIL4
FLP1-FLP4
```

Display only deflection values in the normal dashboard.

Do not clutter the UI with hinge vectors unless explicitly needed.

---

# 24. Drag-Rise Module

Each configuration may optionally define one drag-rise directory.

Directory:

```text
03-RESULTS/DRAG-RISE/<dir_name>/
```

Files:

```text
drag_rise_cl0p20.dat
drag_rise_cl0p25.dat
drag_rise_clm0p25.dat
```

Filename mapping:

```text
cl0p20   -> +0.20
cl0p25   -> +0.25
clm0p25  -> -0.25
```

File columns include:

```text
POLAR
MACH
REYNOLDS
ALPHA
BETA
...
CDS
CYS
CLS
...
```

For each file:

```text
Delta_CDS(M) = CDS(M) - CDS(M_lowest)
```

Plot:

```text
Delta_CDS vs Mach
```

Create one plot per target CLS.

Overlay configurations on that same CLS plot.

Drag-rise module is optional.

If no drag-rise directory is selected, omit the drag-rise plots.

---

# 25. Detailed Convergence Histories

Selectable case history should include:

## Residuals

- continuity;
- x-velocity;
- y-velocity;
- z-velocity;
- energy;
- nut.

## Aero histories

- CLZB;
- CDXB;
- CMYB.

## Flow histories

- cp-max;
- tstep-ave;
- optionally time/iter.

Linked plot selection should update the selected case when possible.

---

# 26. Data Provenance / Integrity

The dashboard should expose:

- report generation time;
- source base folders;
- script version;
- module versions;
- source-file modification times;
- threshold set;
- validated module status.

Integrity checks should include:

- missing ADF points;
- duplicate ALPHA/BETA;
- non-monotonic sweep order;
- mismatched case counts;
- missing FLUENT_LOG;
- incomplete runs;
- unexpected Mach/Re changes;
- drag-rise filename CLS mismatch;
- missing DISTCLCP station mappings.

---

# 27. Browser Presets / Session Memory

Use browser `localStorage` to retain:

- filters;
- curve style;
- axis selection;
- coefficient;
- moment reference mode;
- thresholds;
- display density;
- comparison pairs;
- saved presets.

Presets may include:

- Longitudinal review;
- Lateral-directional review;
- Convergence audit;
- Drag-rise comparison.

---

# 28. Future Project Setup UX

Current launcher supports saved setup/load setup.

Future recommendation:

```text
New
Open Setup
Save
Save As...
```

Potential future project extension:

```text
*.jamal
```

Internally JSON is acceptable.

A JAMAL project file could eventually store:

- configurations;
- base folders;
- selected POLARs;
- drag-rise dirs;
- moment reference;
- thresholds;
- comparisons;
- display options;
- notes.

This is not implemented yet and can wait.

---

# 29. Planned New Module: DISTCLCP

This is the next major requested feature.

Add a new:

```text
Distributions
```

tab.

Directory:

```text
03-RESULTS/DISTCLCP/POLAR-XXX/<COMPONENT>/
```

Examples:

```text
WING
VTAIL
...
```

The component name must be discovered dynamically.

Do not hard-code `WING`.

---

# 30. DISTCLCP File Types

There are three file families.

## 30.1 Section geometry

Example filename:

```text
section_state1_stationY.YYY
```

The geometry files are used to determine local chord.

Sample format:

```text
*KEYWORD
*DEFINE_CURVE_TITLE
DEFAULT_PLANE_ZX section curve
$     LCID      SIDR       SFA       SFO      OFFA      OFFO    DATTYP
        22         0  1.000000  1.000000  0.000000  0.000000         0
$ (Color) ...
$           ABSCISSA            ORDINATE
        2.327996E+00       -1.375462E-01
        ...
*END
```

Parsing rule:

1. find the line containing:

```text
ABSCISSA            ORDINATE
```

2. data starts on the next line;
3. read two numeric columns;
4. stop at `*END`.

Interpret as:

```text
X, section-coordinate
```

Compute:

```text
Xmin = min(X)
Xmax = max(X)
chord = Xmax - Xmin
```

Do this for every station.

Geometry is based on `section_state1_*` and can be treated as the reference geometry unless later evidence shows state-varying geometry.

---

# 31. total_force_stateS

Filename:

```text
total_force_stateS
```

where:

```text
S = 1, 2, 3, ...
```

No header.

Read rows directly.

Columns:

```text
Y   FX   FY   FZ
```

The values are **force per unit span**.

This is confirmed.

State mapping comes from the sequence in `infout [CASES]`.

For example:

```text
state1 -> first case
state2 -> second case
...
```

Each state therefore gets:

```text
ALPHA
BETA
Mach
Reynolds
q
BREF
CREF
```

from `infout`.

---

# 32. Sectional cl Calculation

Because force is per unit span:

```text
cl(y) = L'(y) / (q_inf * c(y))
```

No `Δy` term.

For beta = 0 and body-axis forces with:

```text
X forward
Z downward
```

positive lift per unit span should be computed as:

```text
L' = FX*sin(alpha) - FZ*cos(alpha)
```

Sanity check at alpha = 0:

```text
L' = -FZ
```

So an upward body force (`FZ < 0`) gives positive lift.

If beta is later required, generalize using the validated body-to-stability force transformation rather than inventing a separate formula.

---

# 33. Spanwise Coordinate

Default non-dimensional coordinate:

```text
eta = 2Y / BREF
```

For a full wing this gives approximately:

```text
-1 <= eta <= +1
```

Recommended UI options:

```text
2Y/BREF   (default)
Y/BREF
dimensional Y
```

For future VTAIL support, coordinate normalization may need component-specific generalization, but do not guess before seeing actual VTAIL data.

---

# 34. Cp Distribution Files

Example filename:

```text
cp_dist_stateS_stationY.YYY
```

Format is the same LS-DYNA-style two-column curve format as the section file:

```text
$           ABSCISSA            ORDINATE
X                              Cp
...
*END
```

Use the same generic two-column curve parser.

Interpret:

```text
ABSCISSA = X
ORDINATE = Cp
```

For the corresponding station, use geometry:

```text
X_LE_candidate = Xmin
c = Xmax - Xmin
```

Initial normalization:

```text
x/c = (X - Xmin)/c
```

But verify leading-edge direction with actual geometry before permanently assuming `Xmin` is the LE.

Final display convention must always be:

```text
x/c = 0 -> leading edge
x/c = 1 -> trailing edge
```

Cp plot:

```text
Cp vs x/c
```

Use inverted Cp y-axis:

```text
more negative Cp upward
```

---

# 35. DISTCLCP Station Matching

Station values appear in filenames:

```text
stationY.YYY
```

Force file also has Y values.

Do not use exact string matching.

Use a small numeric station tolerance.

Integrity checks should flag:

- section exists but no force station;
- force station has no section;
- Cp station has no section geometry;
- duplicate stations;
- non-positive chord;
- missing state;
- state count inconsistent with infout;
- mismatched station values.

---

# 36. Proposed Distributions Tab

Controls:

```text
Configuration: [Baseline ▼]
POLAR:         [POLAR-014 ▼]
Component:     [WING ▼]
State:         [State 5 · alpha=4.00° ▼]
```

## Spanwise load plot

```text
cl vs 2Y/BREF
```

Allow overlay of:
- multiple states;
- multiple configurations;
- optionally multiple POLARs.

## Chordwise Cp plot

Station selection:

```text
☑ Y=-5.000
☐ Y=-4.500
☑ Y=-4.000
...
```

Controls:

```text
Select all
Clear
Every nth station
Single-station mode
```

Hover should include:

```text
Configuration
POLAR
Component
State
ALPHA/BETA
Station Y
eta
x/c
Cp
```

---

# 37. Planned Distribution Delta Functions

Future comparison support:

```text
Δcl(y) = comparison - reference
ΔCp(x/c) = comparison - reference
```

For Cp comparison, interpolate the comparison curve onto the reference `x/c` grid, following the same philosophy used in the Delta tab.

---

# 38. DISTCLCP Caching

Cache separately:

## Geometry cache

```text
section_state1_station*.*
```

## Per-state cache

```text
total_force_stateS
cp_dist_stateS_station*.*
```

Fingerprint by:

```text
path + size + modification timestamp
```

Only new/modified files should be reparsed.

---

# 39. Protected / Validated Behaviors

Do not casually rewrite these:

1. Fluent 11-line subblock parser fix.
2. v15-derived validated moment-reference transfer.
3. Correct %MAC -> XNEW conversion.
4. infout frame to body-frame sign conversion.
5. B/S/W coefficient selection behavior.
6. Drag-rise delta-CD definition.
7. CL-based Delta interpolation for CD/CM/L/D/SM.
8. incremental cache fingerprint logic.
9. shallow ADF discovery that ignores `ADF_COMP`.
10. final residual threshold set unless user changes it.

---

# 40. Known Historical Pitfalls

These bugs happened before and must not be reintroduced:

## Wrong XNEW at 0% MAC

Incorrect behavior once produced:

```text
XNEW = -1.1255
dx = -12.2535
```

Correct example:

```text
XREF = 11.128
CREF = 4.629

XNEW(0% MAC) = 9.97075
dx_info = -1.15725
```

## Wrong CREF parsing

At one point the parser accidentally read the `2` from:

```text
SREF[m2]
```

as a number, shifting positional parsing.

Use robust keyed parsing.

## Fluent subblock bug

Repeated monitor headers caused:

```text
actual_iters = 11
```

for every case.

Must concatenate subblocks until actual case boundary / expected count mapping.

## cp-max early spike misclassification

Early cp-max spikes should not automatically cause `DIVERGED`.

## Moment reference direct-stability approach

A direct stability-axis transfer path caused sign/rotation confusion.

The validated implementation uses the protected body-first method.

## Static margin regressions

Several attempts to rewrite the static-margin derivative caused incorrect results.

Preserve the stable implementation unless explicitly revalidated.

---

# 41. Recommended Module Versioning

Maintain logical module versions independently.

Example:

```text
ADF Parser................v1.x
infout Parser.............v1.x
Fluent Parser.............v1.x
Moment Transfer...........v1.0 VALIDATED
Static Margin.............v1.x
Convergence Engine........v2.x
Mesh Quality..............v1.x
Outlier Detection.........v1.x
Drag Rise.................v1.x
DISTCLCP..................v0.x until validated
HTML Dashboard............vX
Launcher / Cache..........vX
```

Expose module versions in the provenance/Quality panel.

---

# 42. Development Guidance for Codex

Before editing:

1. inspect the current project files;
2. locate the v25 launcher;
3. locate the current dashboard engine used by v25;
4. identify the validated moment-transfer function;
5. identify the cache manifest format;
6. run existing tests/smoke checks if present;
7. preserve backward compatibility where practical.

Do not rebuild the project from scratch.

Prefer small, testable patches.

For every significant change:

- Python compile check;
- embedded JavaScript syntax check;
- minimal synthetic parser test;
- regression check on moment transfer;
- regression check on Fluent subblock parsing;
- regression check on drag rise;
- regression check on cache reuse;
- regression check on Delta interpolation.

---

# 43. Immediate Next Development Task

Implement the new **Distributions** module using the specifications above.

Recommended order:

1. generic LS-DYNA-style two-column curve parser;
2. shallow component discovery under `DISTCLCP/POLAR-XXX`;
3. station filename parser;
4. section geometry/chord extraction;
5. `total_force_stateS` parser;
6. state mapping through `infout [CASES]`;
7. sectional `cl(y)`;
8. Cp parser + x/c normalization;
9. cache integration;
10. Distributions tab UI;
11. integrity checks;
12. multi-state/config overlay;
13. later Delta-cl and Delta-Cp comparisons.

Do not change the validated moment-reference module while implementing DISTCLCP.

---

# 44. Summary of User Preferences

The user prefers:

- engineering-focused UI;
- compact but readable controls;
- portable standalone HTML;
- local Python launcher for filesystem access;
- fast shallow discovery;
- incremental update rather than full reprocessing;
- explicit progress indication;
- clear traceability;
- multiple configurations;
- interactive Plotly plots;
- scientifically meaningful axes and conventions;
- consistent curve styling;
- minimal file clutter;
- separate numerical and aerodynamic convergence concepts;
- validated physics over clever but unverified transformations.

The project is intended for serious CFD engineering use, not as a demo.
