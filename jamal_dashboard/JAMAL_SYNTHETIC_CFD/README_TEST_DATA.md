# Synthetic JAMAL CFD test dataset

Entirely synthetic parser/dashboard fixtures, patterned on the infout, Fluent monitor, drag-rise and LS-DYNA curve examples in the referenced CFD plotter conversation. No production CFD data, proprietary geometry, mesh, solver restart, or real simulation results are included.

## Getting started

Extract the ZIP and select **JAMAL_SYNTHETIC_CFD** as the JAMAL base folder. Select POLAR-001 and POLAR-002. For drag rise, select **03-RESULTS/DRAG-RISE/BASELINE**. For distributions, select either POLAR, component WING, and states 1 through 5.

The supplied `validate_test_data.py` uses only the Python 3 standard library:

```text
python validate_test_data.py
```

Run it from the extracted folder, or pass the base-folder path as its first argument. It checks the raw files and compares them with `expected_values.json`. The supplied validation report records the successful run. The dashboard's v25 source/package was not available in the retrieved attachments, so an actual v25 dashboard generation was not run. Fixture validation does not guarantee every parser version accepts the files unchanged.

## Folder structure

```text
JAMAL_SYNTHETIC_CFD/
  README_TEST_DATA.md
  expected_values.json
  validate_test_data.py
  VALIDATION_REPORT.txt
  SHA256SUMS.txt
  00-SUPPORT/.keep
  01-GRIDS/.keep
  02-RUNS/
    POLAR-001/{infout,FLUENT_LOG}
    POLAR-002/{infout,FLUENT_LOG}
  03-RESULTS/
    ADF/
      POLAR-001.adf
      POLAR-002.adf
      ADF_COMP/WING/POLAR-999.adf
    DRAG-RISE/BASELINE/
      drag_rise_cl0p20.dat
      drag_rise_cl0p50.dat
    DISTCLCP/
      POLAR-001/WING/
        section_state1_stationY.YYY       (9 files)
        total_force_stateS                (5 files)
        cp_dist_stateS_stationY.YYY       (45 files)
      POLAR-002/WING/                     (same inventory)
    CONVERGENCE/.keep
    FLOWVIS/.keep
    PROBE/.keep
```

The letters S and Y.YYY above are placeholders; actual filenames include `total_force_state3`, `section_state1_station-2.500`, and `cp_dist_state3_station2.500`. Distribution files have no extension. DASHBOARD is initially absent so output-folder creation can be tested. `.keep` files only preserve empty project folders; the grid path in infout is illustrative and has no corresponding mesh.

There are 127 CFD input files: 2 infout, 2 FLUENT_LOG, 3 ADF including the ignored decoy, 2 drag-rise tables, 18 geometry curves, 10 force tables, and 90 Cp curves.

## States, reference data and flow

State S is the S-th row in infout `[CASES]`, starting at 1. Do not sort cases by angle before mapping distributions or histories. Angles are degrees; both sweeps have BETA = 0.

| POLAR | Mach | Reynolds | qdin, Pa | State 1 to 5 ALPHA, degrees |
|---|---:|---:|---:|---|
| POLAR-001 | 0.20 | 7.0518518519E+06 | 2774.4 | 0, 2, 4, 6, -1 |
| POLAR-002 | 0.70 | 2.4681481481E+07 | 33986.4 | 0, 1, 3, 5, -1 |

Common references: SREF = 15 m2, BREF = 10 m, CREF = 14/9 m (1.5555555556), XREF = 3 m, YREF = 0 m, ZREF = +0.1 m. Density = 1.2 kg/m3, sound speed = 340 m/s, viscosity = 1.8E-05 Pa.s. V = Mach * 340; qdin = 0.5 * rho * V^2; Reynolds = rho * V * CREF / mu. Use **qdin**, not qimp, as dynamic pressure. Pressure/temperature and total conditions use gamma = 1.4 and R = 287.05.

infout coordinates follow the earlier geometry/reference convention: X increases aft, Z increases upward. Physical body-axis forces in the force tables use X forward, Z downward. Do not apply the reference-coordinate signs directly to the force columns. ADF moments are supplied at the listed reference; no extra moment-transfer fixture is implied by the coefficient suffix `25`.

## ADF and coefficient sanity values

Each root ADF is a whitespace table with five rows and these 22 columns:

```text
MACH REYNOLDS ALPHA BETA
CDB CYB CLB CRB25 CMB25 CNB25
CDS CYS CLS CRS25 CMS25 CNS25
CDW CYW CLW CRW25 CMW25 CNW25
```

CDB is positive backward and CLB positive upward. For alpha in radians:

```text
CDB = CDS*cos(alpha) - CLS*sin(alpha)
CLB = CDS*sin(alpha) + CLS*cos(alpha)
CMS25 = CMB25 = CMW25 = -0.02 - 0.10*CLS
```

Stability and wind sets coincide because BETA is zero. Side force, roll and yaw are zero by symmetry. A negative CDB at positive alpha is intentional; stability drag remains positive.

| POLAR | State-order CLS values |
|---|---|
| POLAR-001 | 0.26625, 0.44375, 0.62125, 0.79875, 0.17750 |
| POLAR-002 | 0.26625, 0.3683125, 0.5724375, 0.7765625, 0.1641875 |

POLAR-001 state 1: CDS = 0.0228355625, CMS25 = -0.046625. POLAR-001 uses CDS = 0.020 + 0.04*CLS^2; POLAR-002 uses CDS = 0.023 + 0.04*CLS^2. Static margin `-dCMS25/dCLS` is 0.10, or 10% CREF, after sorting appropriately for differentiation. Machine-readable values for every state and station are in `expected_values.json`.

Shallow ADF discovery must return exactly POLAR-001 and POLAR-002. `ADF_COMP/WING/POLAR-999.adf` is a valid synthetic component table, deliberately without a matching run folder. Recursive discovery incorrectly exposes a third POLAR. Do not process that decoy as a whole-aircraft run.

## Fluent history checks

Each FLUENT_LOG has five ALPHA/BETA banners, 1,200 numeric monitor rows, 60 repeated monitor headers and 10 solve subblocks. Each state has 240 rows split into two 120-iteration solve commands; the header repeats every 20 rows. Concatenate all 12 headers' data into one state history. Headers and solve commands do not create states. Iterations remain global across the POLAR: 1-240, 241-480, 481-720, 721-960, 961-1200.

The monitor names and time/iteration trailer follow the supplied Fluent sample:

```text
iter continuity x-velocity y-velocity z-velocity energy nut
tstep-ave cp-max cnzb cmyb crxb clzb cyyb cdxb time/iter
```

Each row ends with a time token and remaining-iterations integer. Ignore both as aerodynamic fields. `cdxb` maps to CDB, `clzb` to CLB, `cmyb` to CMB25, `crxb` to CRB25, and `cnzb` to CNB25. Final monitor coefficients agree with the ADF row within text-rounding tolerance.

Expected final continuity residual: 2E-04 for every POLAR-001 state; 2E-03 for POLAR-002 states 1-4; 9E-03 for POLAR-002 state 5. All other residuals finish below continuity. These exercise the previously discussed 8E-04/8E-03/1E-02 residual boundaries; exact overall dashboard quality labels can also depend on other criteria and are not asserted here.

Every state has an early cp-max spike of 3.6 at its second row, followed by recovery to exactly 1.05 at its final row. A classification based solely on the global maximum would incorrectly penalize recovered states. The plateauing coefficient histories test convergence trends separately from residual thresholds. POLAR-001 uses CRLF in infout/log/Cp files and POLAR-002 uses LF; both use CRLF geometry curves.

## Geometry, force and Cp checks

Nine full-span stations per POLAR: -5, -3.75, -2.5, -1.25, 0, 1.25, 2.5, 3.75, 5 m. `eta = 2Y/BREF` runs from -1 to +1 in steps of 0.25.

```text
chord(Y) = 2 - 0.2*abs(Y) m
Xmin(Y) = 2 + 0.15*abs(Y) m
Xmax(Y) = Xmin(Y) + chord(Y)
cl(Y) = cl_root * (1 - 0.4*eta^2)
POLAR-001 cl_root = 0.30 + 0.100*ALPHA
POLAR-002 cl_root = 0.30 + 0.115*ALPHA
```

The geometry has a symmetric 12%-thick section, sweep and mild dihedral. Each geometry/Cp file has exactly **81 numeric curve points after ABSCISSA/ORDINATE**, ending at `*END`. The earlier numeric LCID/scaling row is metadata and must never enter the point array. Geometry order is upper trailing edge to leading edge (41 points), then lower surface to trailing edge (40 additional points). Geometry state1 files are reused for all five states.

`total_force_stateS` is headerless: Y, FX, FY, FZ, with forces **per unit span in N/m**. No delta-Y factor belongs in the sectional coefficient:

```text
Lprime = FX*sin(alpha) - FZ*cos(alpha)
Dprime = -FX*cos(alpha) - FZ*sin(alpha)
cl = Lprime / (qdin*chord)
```

At POLAR-001 state1 / Y = 0: chord = 2 m, eta = 0, FX = -126.7099692 N/m, FY = 0, FZ = -1664.64 N/m, Lprime = 1664.64 N/m, cl = 0.30. At Y = +/-2.5: chord = 1.5 m and cl = 0.27. At Y = +/-5: chord = 1 m and cl = 0.18. All local lift values are positive, including the negative-alpha states of this cambered aerodynamic model.

The force row for positive Y = 2.5 stores **2.4999998** deliberately. Match to filename station2.500 with absolute tolerance 1E-06 m, then use the canonical geometry stations for integration. Trapezoidal span integration of chord gives SREF = 15 m2. Integrating Lprime or Dprime and dividing by qdin*SREF reproduces ADF CLS/CDS. The span-weighted shape factor on these nine stations is 0.8875, so CLS = 0.8875*cl_root. Finite positive loads at the last station are an intentional coarse-grid approximation; this is not an exact tip solution.

Cp abscissae are dimensional geometry X, not x/c. Normalize as `(X-Xmin)/chord`. Cp uses the same point order as geometry and has negative upper-surface suction, positive lower-surface loading, leading-edge Cp = 1, and pressure recovery toward zero at the trailing edge. Plot with the Cp vertical axis inverted. At POLAR-001 state1 / Y = 0, minimum upper Cp is approximately **-0.412893343**.

For an additional cross-file check, trapezoidal integration of `(Cp_lower-Cp_upper)` over x/c equals the local body normal coefficient `-FZ/(qdin*chord)`. At alpha = 0 this is cl; at other angles it includes the drag projection. Cp curves are analytically constructed test curves, not a compressible CFD solution; full pressure integration around the finite-thickness surface and viscous drag are not modeled. Preserving upper/lower branches is essential: sorting all Cp points by X would mix both surfaces.

## Drag rise

Both files have seven rows at Mach = 0.30, 0.50, 0.60, 0.70, 0.75, 0.80, 0.85. The column layout is the ADF layout with an initial POLAR column. POLAR identifiers 101-107 are provenance labels for synthetic Mach samples and do not imply additional run folders. CLS is exactly 0.20 or 0.50 in its respective file.

```text
CDS = 0.018 + 0.04*CLS^2 + 0.002*(Mach-0.30)
      + 0.8*max(Mach-0.68,0)^2
Delta_CD = CDS(Mach) - CDS(lowest Mach)
```

| Mach | Delta_CD, both files |
|---:|---:|
| 0.30 | 0 |
| 0.50 | 0.00040 |
| 0.60 | 0.00060 |
| 0.70 | 0.00112 |
| 0.75 | 0.00482 |
| 0.80 | 0.01252 |
| 0.85 | 0.02422 |

The baseline CDS values are 0.0196 and 0.0280; final values are 0.04382 and 0.05222. Each constant-CLS file should form its own drag-rise plot. Their Delta_CD curves intentionally coincide, which makes baseline subtraction easy to verify.

## Suggested dashboard checks

1. Shallow discovery exposes two POLARs and the BASELINE drag-rise directory without opening run files.
2. Select one POLAR, generate, then add the other. With unchanged inputs, incremental caching should reuse the first and parse the second.
3. Repeated monitor subblocks concatenate to 240 samples per state; angle banners preserve state5 = -1 degree.
4. Coefficient, static-margin, drag-rise and distribution plots agree with the values above. State labels carry the matching ALPHA.
5. Geometry discovery finds WING under either POLAR, and only state1 supplies geometry. Station selectors expose nine stations including negative values.
6. Use a disposable copy for missing-force, missing-geometry, altered-chord or missing-state fault injection. The shipped fixture is complete; no missing-file warning is expected for selected root POLARs.

There is one configuration, two ALPHA sweeps, and one component. Nonzero-BETA frame rotation, multiple-configuration Delta plots, VTAIL geometry, reversed geometric X direction, and malformed-file recovery are outside this fixture's coverage. SHA256SUMS.txt covers all packaged files except itself.
