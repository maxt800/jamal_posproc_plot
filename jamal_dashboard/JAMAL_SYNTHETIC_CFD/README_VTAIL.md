# Paired V-tail synthetic extension

This additive fixture implements the requested **paired V-tail**, with a **5 m
projected span**, half the existing 10 m wing span. The assumed cant is **45 degrees**
on each side. These are illustrative distributions, not a validated aerodynamic
design or production CFD result.

The original wing, infout, Fluent logs, whole-aircraft ADF, drag rise, expected
values, and SHA256SUMS.txt are unchanged. The VTAIL loads have **not** been added
to the existing whole-aircraft ADF coefficients. The original README_TEST_DATA.md
and validator describe the original WING-only dataset.

## Inventory and coordinates

Each of POLAR-001 and POLAR-002 now contains a `VTAIL` folder under DISTCLCP:

- Nine section_state1 geometry curves, with 81 points each.
- Five total_force_stateS tables, with nine stations each.
- 45 cp_dist_stateS_station* curves, with 81 points each.
- component_reference.json: projected span reference 5 m and documented conventions.

The nine projected Y stations run from -2.5 m to +2.5 m in 0.625 m steps.
The panel center height is Z = abs(Y)*tan(45 degrees). Each actual panel span is
2.5/cos(45 degrees) = 3.5355 m. The chord tapers from 1.2 m at the centerline to
0.6 m at either tip; projected area is 4.5 m2. The section geometry stores X and
local section-normal coordinates, rather than a full 3D surface.

All plotted coordinates use aircraft infout BREF = 10 m. In v25.3, VTAIL tips
are -50%/+50% for 2Y/BREF, -25%/+25% for Y/BREF, or -2.5/+2.5 m for Y.
The component reference records the physical 5 m projected span only; it does
not override plotting normalization.

## Force and Cp construction

State order and qdin are taken from the original infout fixture. Forces use body
X forward and Z downward, expressed **per projected Y metre**. Let eta = 2Y/5:

```text
cn = (0.18 + 0.035*alpha_deg) * (1 + 0.1*(Mach-0.2)) * (1-0.4*eta^2)
cd = 0.012 + 0.025*cn^2
N_per_projected_Y = qdin*chord*cn/cos(45 degrees)
FX = -qdin*chord*cd/cos(45 degrees)
FY = -sign(Y)*N_per_projected_Y*sin(45 degrees)
FZ = -N_per_projected_Y*cos(45 degrees)
cl_vertical = (FX*sin(alpha)-FZ*cos(alpha))/(qdin*chord)
```

FY at the centerline is zero, representing the average of the two joining panels.
The left/right side forces cancel in span integration. The 1/cos(cant) conversion
accounts for the difference between panel-span and projected-Y force density.
The displayed cl is the vertical lift coefficient per projected-Y metre; it is
not the local panel-normal coefficient at nonzero alpha.

Cp keeps upper-then-lower surface order. Its pressure difference is normalized so
that trapezoidal integration over x/c equals cn, also equal to -FZ/(qdin*chord)
for this projected-span convention. LE Cp is 1 and TE Cp is zero. The drag term
is prescribed separately; full 3D pressure-force integration is not modeled.

POLAR-001 state 1 at the root: qdin=2774.4 Pa, chord=1.2 m, cl=0.18,
FZ=-599.2704 N/m and FY=0. Tips have cl=0.108. State 5 remains alpha=-1 degree.
The +1.25 m force station is stored as 1.2499998 to exercise tolerant matching.

## Reproduction and checks

`expected_vtail_values.json` contains per-state, per-station expected values.
`VTAIL_SHA256SUMS.txt` covers the generated tail inputs and expected-value file.
From the application folder, regenerate only the extension with:

```text
python generate_vtail_fixture.py JAMAL_SYNTHETIC_CFD
```

The development tests compare all 90 tail sectional coefficients, pressure minima,
Cp integrals, side-force cancellation, station percentages and cache reuse. The
original supplied validator still passes all 1,654 WING/ADF/Fluent/drag-rise checks.
