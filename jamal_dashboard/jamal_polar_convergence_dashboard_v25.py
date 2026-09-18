#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JAMAL Aerodynamic Results Dashboard
==========================================

Spyder-friendly script that combines:

1) ADF aerodynamic polar plotting
   - Reads POLAR-XXX.adf files from the configured 03-RESULTS/ADF folders.
   - Plots stability-axis coefficients: CDS, CYS, CLS, CRS25, CMS25, CNS25.
   - Computes static margin from finite difference:
       Static Margin = -dCMS25/dCLS

2) FLUENT_LOG convergence analysis
   - Reads infout and FLUENT_LOG from the corresponding 02-RUNS/POLAR-XXX folders.
   - Uses infout CASE table as the master list of cases and expected iterations.
   - Fixes Fluent's repeated 11-line monitor subblock issue by reading all numeric
     monitor rows continuously and splitting them according to infout ITERS.
   - Classifies each ALPHA/BETA case as CONVERGED, SUSPICIOUS, or DIVERGED.

Output:
    CFD/03-RESULTS/DASHBOARD/dashboard.html
    CFD/03-RESULTS/DASHBOARD/dashboard.json

How to use in Spyder:
    1. Edit USER CONFIGURATION below.
    2. Press Run.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import jamal_distributions


# =============================================================================
# USER CONFIGURATION - EDIT THIS SECTION IN SPYDER
# =============================================================================

CASES = [
    {
        "label": "Baseline",
        "directory": r"Z:\cae\jeniah2\cfd\Wing_optim\DBank\MU11_209_663462\CFD\03-RESULTS\ADF",
        "polars": [14, 15, 16, 18, 19, 20, 21],
        # Optional drag-rise folder under 03-RESULTS/DRAG-RISE/<name>
        # "drag_rise_dir": "my_drag_rise_folder",
    },
    {
        "label": "Config_2",
        "directory": r"Z:\cae\jeniah2\cfd\Wing_optim\DBank\MU11_209_663462_PLAIN_FLAP\CFD\03-RESULTS\ADF",
        "polars": [14, 15, 16, 18, 19, 20, 21],
        # Optional drag-rise folder under 03-RESULTS/DRAG-RISE/<name>
        # "drag_rise_dir": "my_drag_rise_folder",
    },
]

OUTPUT_NAME = "dashboard"

# Intended maximum number of configurations to keep the dashboard readable.
# The script accepts 1 to MAX_CONFIGURATIONS entries in CASES.
MAX_CONFIGURATIONS = 5

# If True, also save a compact static-margin CSV.
SAVE_STATIC_MARGIN_CSV = False

# Plot static margin in percent of CREF.
STATIC_MARGIN_AS_PERCENT = True

# Optional Y-axis limits for static-margin plots in HTML.
# Use percent values when STATIC_MARGIN_AS_PERCENT = True.
# Set to None for automatic scale.
STATIC_MARGIN_YLIM = (-10.0, 20.0)

# Optional hard filtering of static-margin points before plotting.
FILTER_STATIC_MARGIN_OUTLIERS_FOR_PLOT = False

# If True, each POLAR is plotted as a separate curve.
# If False, all selected polars of one configuration are combined into one curve.
PLOT_EACH_POLAR_SEPARATELY = True

# Residual thresholds
RESIDUAL_SUCCESSFUL = 8.0e-4
RESIDUAL_ACCEPTABLE = 8.0e-3
RESIDUAL_BAD = 1.0e-2

# Final-window stability thresholds
CL_WINDOW_RANGE_MAX = 0.003
CD_WINDOW_RANGE_MAX = 0.0005
CM_WINDOW_RANGE_MAX = 0.001

# cp-max thresholds
CPMAX_WINDOW_RANGE_MAX = 0.25
CPMAX_SPIKE_RATIO_MAX = 2.0

# Neighbor consistency thresholds
NEIGHBOR_CL_JUMP_MAX = 0.12
NEIGHBOR_CD_JUMP_MAX = 0.03
NEIGHBOR_CM_JUMP_MAX = 0.05

# The HTML embeds selected-case histories. Limit points per case to keep file size reasonable.
MAX_HISTORY_POINTS_PER_CASE = 1200

# Possible Fluent log file names
FLUENT_LOG_NAMES = [
    "FLUENT_LOG",
    "fluent.log",
    "FLUENT.LOG",
    "FLUENT_LOG.txt",
]

# Module versions shown in the dashboard and JSON output.
SCRIPT_VERSION = "v25.6"
MODULE_VERSIONS = {
    "infout parser": "1.3",
    "Distributions": jamal_distributions.VERSION,
    "ADF parser": "1.2",
    "Fluent history parser": "1.5",
    "Convergence engine": "2.1",
    "Aerodynamic stability assessment": "1.1",
    "Neighbor outlier detector": "1.2",
    "Mesh-quality integration": "1.0",
    "Drag-rise parser": "1.0",
    "Moment-reference transfer": "1.0 validated",
    "Configuration comparison": "3.0",
    "Curve-style engine": "1.0",
    "HTML dashboard": "5.9",
}

# Fluent's final-window metrics use at most the last 200 printed iterations.
FINAL_WINDOW_MAX_ITERS = 200

# Mesh-quality advisory thresholds. Aspect ratio is reported but not used as an
# automatic failure because high boundary-layer aspect ratios can be intentional.
MESH_MIN_ORTHOGONAL_WARNING = 0.01


# =============================================================================
# CONSTANTS
# =============================================================================

HISTORY_COLUMNS = [
    "iter",
    "continuity",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "energy",
    "nut",
    "tstep-ave",
    "cp-max",
    "cnzb",
    "cmyb",
    "crxb",
    "clzb",
    "cyyb",
    "cdxb",
    "time/iter",
]

RESIDUAL_COLUMNS = [
    "continuity",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "energy",
    "nut",
]

AERO_HISTORY_COLUMNS = [
    "cnzb",
    "cmyb",
    "crxb",
    "clzb",
    "cyyb",
    "cdxb",
]

STABILITY_COEFF_CANDIDATES = ["CDS", "CYS", "CLS", "CRS25", "CMS25", "CNS25"]
ADF_COEFF_CANDIDATES = [
    "CDB", "CYB", "CLB", "CRB25", "CMB25", "CNB25",
    "CDS", "CYS", "CLS", "CRS25", "CMS25", "CNS25",
    "CDW", "CYW", "CLW", "CRW25", "CMW25", "CNW25",
]


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class CaseConfig:
    label: str
    directory: Path
    polars: List[int]
    # Optional. If provided, the script reads drag-rise files from:
    #   <CFD>/03-RESULTS/DRAG-RISE/<drag_rise_dir>/drag_rise_cl*.dat
    drag_rise_dir: Optional[str] = None


@dataclass
class PolarData:
    case_label: str
    polar_number: int
    path: Path
    data: pd.DataFrame
    sweep_var: str

    @property
    def polar_name(self) -> str:
        return f"POLAR-{self.polar_number:03d}"

    @property
    def curve_label(self) -> str:
        return f"{self.case_label} | {self.polar_name}"


# =============================================================================
# BASIC UTILITIES
# =============================================================================

def safe_float(text):
    try:
        value = float(str(text).replace(",", ""))
        return value if math.isfinite(value) else None
    except Exception:
        return None


def percentile(values, p):
    if not values:
        return None
    values = sorted(values)
    k = (len(values) - 1) * p / 100.0
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return values[int(k)]
    return values[f] * (c - k) + values[c] * (k - f)


def json_safe(obj):
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, (np.floating, np.integer)):
        val = obj.item()
        if isinstance(val, float):
            return val if math.isfinite(val) else None
        return val
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    return obj


def polar_name_from_number(polar):
    return f"POLAR-{int(polar):03d}"


def polar_filename(polar_number: int) -> str:
    return f"POLAR-{polar_number:03d}.adf"


def get_paths_from_adf_dir(adf_dir):
    """
    Input:
        .../CFD/03-RESULTS/ADF

    Output:
        results_dir     = .../CFD/03-RESULTS
        runs_dir        = .../CFD/02-RUNS
        dashboard_dir = .../CFD/03-RESULTS/DASHBOARD
    """
    adf_dir = Path(adf_dir)
    results_dir = adf_dir.parent
    cfd_dir = results_dir.parent
    runs_dir = cfd_dir / "02-RUNS"
    dashboard_dir = results_dir / "DASHBOARD"
    return results_dir, runs_dir, dashboard_dir


def find_fluent_log(polar_path):
    for name in FLUENT_LOG_NAMES:
        candidate = polar_path / name
        if candidate.is_file():
            return candidate
    transcripts = [p for p in polar_path.glob('*') if p.is_file() and p.suffix.lower() == '.trn']
    if transcripts:
        return max(transcripts, key=lambda p: (p.stat().st_mtime_ns, p.name.lower()))
    matches = [p for p in list(polar_path.glob("*LOG*")) + list(polar_path.glob("*log*")) if p.is_file()]
    if matches:
        return matches[0]
    return polar_path / "FLUENT_LOG"


def downsample_history(rows, max_points=MAX_HISTORY_POINTS_PER_CASE):
    if len(rows) <= max_points:
        return rows
    step = max(1, int(math.ceil(len(rows) / max_points)))
    sampled = rows[::step]
    if sampled[-1] is not rows[-1]:
        sampled.append(rows[-1])
    return sampled


def normalize_cases(raw_cases: Sequence[Dict]) -> List[CaseConfig]:
    return [
        CaseConfig(
            label=str(item["label"]),
            directory=Path(item["directory"]).expanduser(),
            polars=[int(p) for p in item["polars"]],
            drag_rise_dir=item.get("drag_rise_dir"),
        )
        for item in raw_cases
    ]


# =============================================================================
# ADF READING AND STATIC MARGIN
# =============================================================================

def is_float_token(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def find_table_header(lines: Sequence[str]) -> Tuple[int, List[str]]:
    required = {"MACH", "REYNOLDS", "ALPHA", "BETA"}
    for idx, line in enumerate(lines):
        cols = line.strip().split()
        if required.issubset(set(cols)):
            return idx, cols
    raise ValueError("Could not find aerodynamic coefficient table header.")


def read_polar_file(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    header_idx, columns = find_table_header(lines)
    ncols = len(columns)
    rows = []
    for line in lines[header_idx + 1:]:
        tokens = line.strip().split()
        if len(tokens) < ncols:
            continue
        if not all(is_float_token(tok) for tok in tokens[:ncols]):
            continue
        rows.append([float(tok) for tok in tokens[:ncols]])
    if not rows:
        raise ValueError(f"No numeric coefficient data found in {path}")
    return pd.DataFrame(rows, columns=columns)


def detect_sweep_variable(df: pd.DataFrame, tol: float = 1.0e-10) -> str:
    alpha_range = df["ALPHA"].max() - df["ALPHA"].min()
    beta_range = df["BETA"].max() - df["BETA"].min()
    if alpha_range > tol and beta_range <= tol:
        return "ALPHA"
    if beta_range > tol and alpha_range <= tol:
        return "BETA"
    return "ALPHA" if alpha_range >= beta_range else "BETA"


def load_all_adf_polars(cases: Sequence[CaseConfig]) -> List[PolarData]:
    polar_data: List[PolarData] = []
    for case in cases:
        for polar_number in case.polars:
            path = case.directory / polar_filename(polar_number)
            try:
                df = read_polar_file(path)
                sweep_var = detect_sweep_variable(df)
                polar_data.append(
                    PolarData(
                        case_label=case.label,
                        polar_number=polar_number,
                        path=path,
                        data=df,
                        sweep_var=sweep_var,
                    )
                )
                print(f"  ADF loaded: {case.label} | POLAR-{polar_number:03d} | {len(df)} points")
            except Exception as exc:
                print(f"  WARNING: ADF skipped: {path} | {exc}")
    return polar_data


def sort_by_x(df: pd.DataFrame, xcol: str) -> pd.DataFrame:
    return df.sort_values(by=xcol).reset_index(drop=True)


def stability_coefficient_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in STABILITY_COEFF_CANDIDATES if c in df.columns]


def compute_static_margin_dataframe(polar: PolarData, cl_col="CLS", cm_col="CMS25") -> pd.DataFrame:
    required = [cl_col, cm_col, polar.sweep_var, "MACH", "REYNOLDS"]
    missing = [c for c in required if c not in polar.data.columns]
    if missing:
        raise ValueError(f"Missing columns for static margin: {missing}")

    df = polar.data.copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=[cl_col, cm_col, polar.sweep_var])
    df = sort_by_x(df, cl_col)

    if len(df) < 3:
        raise ValueError("Need at least 3 points for second-order finite difference.")

    if df[cl_col].duplicated().any():
        agg_cols = {cm_col: "mean", polar.sweep_var: "mean", "MACH": "mean", "REYNOLDS": "mean"}
        df = df.groupby(cl_col, as_index=False).agg(agg_cols)
        df = sort_by_x(df, cl_col)

    if len(df) < 3:
        raise ValueError("Need at least 3 unique CLS points for static margin.")

    cl = df[cl_col].to_numpy(dtype=float)
    cm = df[cm_col].to_numpy(dtype=float)

    if np.any(np.diff(cl) <= 0.0):
        raise ValueError("CLS values must be strictly increasing after duplicate handling.")

    dcm_dcl = np.gradient(cm, cl, edge_order=2)
    sm = -dcm_dcl

    return pd.DataFrame({
        "case_label": polar.case_label,
        "polar": polar.polar_name,
        "polar_number": polar.polar_number,
        "curve_label": polar.curve_label,
        "sweep_var": polar.sweep_var,
        polar.sweep_var: df[polar.sweep_var].to_numpy(dtype=float),
        "MACH": df["MACH"].to_numpy(dtype=float),
        "REYNOLDS": df["REYNOLDS"].to_numpy(dtype=float),
        "CLS": cl,
        "CMS25": cm,
        "dCMS25_dCLS": dcm_dcl,
        "STATIC_MARGIN": sm,
        "STATIC_MARGIN_PERCENT": 100.0 * sm,
    })


def compute_all_static_margin(polars: Sequence[PolarData]) -> pd.DataFrame:
    frames = []
    for polar in polars:
        try:
            frames.append(compute_static_margin_dataframe(polar))
        except Exception as exc:
            print(f"  WARNING: Static margin skipped for {polar.curve_label}: {exc}")
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def make_adf_plot_rows(polars: Sequence[PolarData], sm_df: pd.DataFrame) -> Dict:
    """Prepare compact JSON-friendly ADF plot data for HTML."""
    curves = []
    for polar in polars:
        df = polar.data.copy().replace([np.inf, -np.inf], np.nan)
        available = [c for c in ADF_COEFF_CANDIDATES if c in df.columns]
        cols = ["MACH", "REYNOLDS", "ALPHA", "BETA"] + available
        rows = df[cols].dropna(how="all").to_dict(orient="records")
        curves.append({
            "case_label": polar.case_label,
            "polar": polar.polar_name,
            "curve_label": polar.curve_label,
            "sweep_var": polar.sweep_var,
            "path": str(polar.path),
            "modified": datetime.fromtimestamp(polar.path.stat().st_mtime, tz=timezone.utc).isoformat() if polar.path.exists() else None,
            "rows": rows,
        })

    sm_rows = []
    if not sm_df.empty:
        sm_rows = sm_df.replace([np.inf, -np.inf], np.nan).to_dict(orient="records")

    return {"curves": curves, "static_margin": sm_rows}


# =============================================================================
# OPTIONAL DRAG-RISE READING
# =============================================================================

def parse_drag_rise_cls_from_filename(path: Path) -> Tuple[str, Optional[float]]:
    """Parse names like drag_rise_cl0p20.dat and drag_rise_clm0p25.dat."""
    m = re.search(r"drag[_-]rise[_-]cl(m?)([0-9]+(?:p[0-9]+)?)", path.stem, re.I)
    if not m:
        return path.stem, None
    sign = -1.0 if m.group(1).lower() == "m" else 1.0
    value = sign * float(m.group(2).replace("p", "."))
    return f"CLS = {value:+.2f}", value


def read_drag_rise_file(path: Path) -> pd.DataFrame:
    """Read one drag-rise .dat table and compute DELTA_CDS from the lowest Mach point."""
    df = read_polar_file(path)
    missing = [c for c in ["MACH", "CDS"] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required drag-rise columns {missing} in {path}")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["MACH", "CDS"])
    df = df.sort_values("MACH").reset_index(drop=True)
    if df.empty:
        raise ValueError(f"No numeric drag-rise data found in {path}")
    cd0 = float(df.iloc[0]["CDS"])
    df["DELTA_CDS"] = df["CDS"] - cd0
    return df


def load_drag_rise_data(cases: Sequence[CaseConfig]) -> Dict:
    """Load optional drag-rise folders for each configuration.

    Each CASES entry may define:
        "drag_rise_dir": "folder_name"

    The script reads:
        <CFD>/03-RESULTS/DRAG-RISE/<folder_name>/drag_rise_cl*.dat
    """
    curves = []
    for case in cases:
        if not case.drag_rise_dir:
            continue
        results_dir, _, _ = get_paths_from_adf_dir(case.directory)
        drag_dir = results_dir / "DRAG-RISE" / str(case.drag_rise_dir)
        if not drag_dir.exists():
            print(f"  WARNING: drag-rise directory not found for {case.label}: {drag_dir}")
            continue
        files = sorted(list(drag_dir.glob("drag_rise_cl*.dat")) + list(drag_dir.glob("drag-rise-cl*.dat")))
        if not files:
            print(f"  WARNING: no drag-rise files found in {drag_dir}")
            continue
        for path in files:
            try:
                cls_label, cls_value = parse_drag_rise_cls_from_filename(path)
                df = read_drag_rise_file(path)
                cols = [c for c in ["POLAR", "MACH", "REYNOLDS", "ALPHA", "BETA", "CDB", "CDW", "CDS", "CLB", "CLW", "CLS", "DELTA_CDS"] if c in df.columns]
                curves.append({
                    "case_label": case.label,
                    "drag_rise_dir": str(case.drag_rise_dir),
                    "path": str(path),
                    "modified": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat() if path.exists() else None,
                    "file_name": path.name,
                    "cls_label": cls_label,
                    "cls_value": cls_value,
                    "rows": df[cols].to_dict(orient="records"),
                })
                print(f"  Drag-rise loaded: {case.label} | {path.name} | {len(df)} points")
            except Exception as exc:
                print(f"  WARNING: drag-rise skipped: {path} | {exc}")
    return {"curves": curves}


# =============================================================================
# INFOUT PARSER
# =============================================================================

def parse_infout(infout_path: Path) -> Dict:
    text = infout_path.read_text(errors="ignore").splitlines()

    meta = {
        "polar": None,
        "aircraft": None,
        "configuration": None,
        "dimension_ref": None,
        "numeric_set": None,
        "mach": None,
        "reynolds": None,
        "sref": None,
        "cref": None,
        "bref": None,
        "xref": None,
        "yref": None,
        "zref": None,
    }

    cases = []
    deflections = {}
    in_cases = False

    for line in text:
        if "POLAR:" in line and meta["polar"] is None:
            m = re.search(r"POLAR:\s*(\S+)", line)
            if m:
                meta["polar"] = m.group(1)
        if "AIRCRAFT:" in line:
            m = re.search(r"AIRCRAFT:\s*(\S+)", line)
            if m:
                meta["aircraft"] = m.group(1)
        if "CONFIGURATION:" in line:
            m = re.search(r"CONFIGURATION:\s*(\S+)", line)
            if m:
                meta["configuration"] = m.group(1)
        if "DIMENSION_REF:" in line:
            m = re.search(r"DIMENSION_REF:\s*(\S+)", line)
            if m:
                meta["dimension_ref"] = m.group(1)
        if "NUMERIC_SET:" in line:
            m = re.search(r"NUMERIC_SET:\s*(\S+)", line)
            if m:
                meta["numeric_set"] = m.group(1)
        if "SREF[m2]:" in line:
            # Important: do not use a generic number regex on the whole line,
            # because it would also capture the "2" in "SREF[m2]".
            # Parse the labeled values explicitly.
            m = re.search(
                r"SREF\[m2\]:\s*([0-9Ee+\-.]+)\s+"
                r"CREF\[m\]:\s*([0-9Ee+\-.]+)\s+"
                r"BREF\[m\]:\s*([0-9Ee+\-.]+)",
                line,
            )
            if m:
                meta["sref"] = float(m.group(1))
                meta["cref"] = float(m.group(2))
                meta["bref"] = float(m.group(3))
        for key, unit in (("XREF", "m"), ("YREF", "m"), ("ZREF", "m"),
                          ("qdin", "Pa"), ("rho", "kg/m3"), ("V", "m/s")):
            match = re.search(rf"\b{key}\[{re.escape(unit)}\]:\s*({jamal_distributions.NUMBER})", line, re.I)
            if match:
                meta[key.lower()] = jamal_distributions.finite_number(match.group(1))
        if "Mach:" in line and "ptot" in line:
            m = re.search(r"Mach:\s*([-+]?\d+\.\d+|[-+]?\d+)", line)
            if m:
                meta["mach"] = float(m.group(1))
        if "Reynolds:" in line:
            m = re.search(r"Reynolds:\s*([-+]?\d+\.\d+E[+-]?\d+|[-+]?\d+\.\d+|[-+]?\d+)", line, re.I)
            if m:
                meta["reynolds"] = float(m.group(1))

        # Control-surface/flap deflection values.
        # Capture only lines such as RUD1:, ELV1:, AIL1:, FLP1:.
        # Ignore geometry lines such as X0R1, DXR1, etc.
        mdef = re.match(r"\s*(RUD|ELV|AIL|FLP)([1-4]):\s*(\S+)", line)
        if mdef:
            key = f"{mdef.group(1)}{mdef.group(2)}"
            val_text = mdef.group(3)
            val = safe_float(val_text)
            deflections[key] = val if val is not None else val_text

        if "[CASES]" in line:
            in_cases = True
            continue

        if in_cases:
            stripped = line.strip()
            if stripped.startswith("CASE"):
                continue
            if not stripped:
                continue
            if stripped.startswith("Grid_files"):
                break
            parts = stripped.split()
            if len(parts) >= 7 and parts[0].isdigit():
                cases.append({
                    "case": parts[0],
                    "mach": safe_float(parts[1]),
                    "reynolds": safe_float(parts[2]),
                    "alpha": safe_float(parts[3]),
                    "beta": safe_float(parts[4]),
                    "name": parts[5],
                    "expected_iters": int(float(parts[6])),
                })

    return {"meta": meta, "cases": cases, "deflections": deflections}


# =============================================================================
# FLUENT LOG PARSER - FIXED FOR REPEATED 11-LINE SUBBLOCK HEADERS
# =============================================================================

def is_history_header(line: str) -> bool:
    low = line.lower()
    tokens = ["iter", "continuity", "x-velocity", "y-velocity", "z-velocity", "energy", "nut"]
    return all(t in low for t in tokens)


def numeric_tokens(line: str) -> List[str]:
    return re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", line)


def parse_history_row(line: str):
    toks = numeric_tokens(line.strip())
    if len(toks) < len(HISTORY_COLUMNS):
        return None

    first = safe_float(toks[0])
    if first is None:
        return None
    if abs(first - round(first)) > 1.0e-9:
        return None

    values = []
    for token in toks[:len(HISTORY_COLUMNS)]:
        val = safe_float(token)
        if val is None:
            return None
        values.append(val)

    return dict(zip(HISTORY_COLUMNS, values))


def parse_fluent_log_all_rows(log_path: Path) -> Dict:
    """
    Read all Fluent monitor rows as one continuous stream.

    Fluent may reprint the monitor header every ~11 rows. Older code treated each
    header as a new block, producing actual_iters = 11. This parser ignores repeated
    headers and collects all numeric monitor rows continuously. Splitting into
    ALPHA/BETA cases is done later using infout expected ITERS.
    """
    lines = log_path.read_text(errors="ignore").splitlines()

    all_rows = []
    headers_found = 0
    first_header_line = None
    first_history_row = None

    mesh_info = {
        "cell_count": None,
        "node_count": None,
        "min_orthogonal_quality": None,
        "max_aspect_ratio": None,
        "orthogonal_quality_history": [],
        "aspect_ratio_history": [],
        "negative_volume_warnings": 0,
        "warnings": [],
        "status": "UNKNOWN",
    }

    for line in lines:
        if mesh_info["cell_count"] is None:
            m = re.match(r"\s*(\d+)\s+cells,", line)
            if m:
                mesh_info["cell_count"] = int(m.group(1))
        if mesh_info["node_count"] is None:
            m = re.match(r"\s*(\d+)\s+nodes,", line)
            if m:
                mesh_info["node_count"] = int(m.group(1))

        if "Minimum Orthogonal Quality" in line:
            m = re.search(r"Minimum Orthogonal Quality\s*=\s*([0-9Ee+\-.]+)", line)
            if m:
                value = safe_float(m.group(1))
                if value is not None:
                    mesh_info["min_orthogonal_quality"] = value
                    mesh_info["orthogonal_quality_history"].append(value)
        if "Maximum Aspect Ratio" in line:
            m = re.search(r"Maximum Aspect Ratio\s*=\s*([0-9Ee+\-.]+)", line)
            if m:
                value = safe_float(m.group(1))
                if value is not None:
                    mesh_info["max_aspect_ratio"] = value
                    mesh_info["aspect_ratio_history"].append(value)

        low_line = line.lower()
        if "negative volume" in low_line or "negative cell volume" in low_line:
            mesh_info["negative_volume_warnings"] += 1
        if "warning" in low_line:
            warning = line.strip()
            if warning and warning not in mesh_info["warnings"]:
                mesh_info["warnings"].append(warning)

        if is_history_header(line):
            headers_found += 1
            if first_header_line is None:
                first_header_line = line.strip()
            continue

        row = parse_history_row(line)
        if row is not None:
            all_rows.append(row)
            if first_history_row is None:
                first_history_row = line.strip()

    min_oq = mesh_info.get("min_orthogonal_quality")
    if mesh_info["negative_volume_warnings"] > 0:
        mesh_info["status"] = "CRITICAL"
    elif min_oq is not None and min_oq < MESH_MIN_ORTHOGONAL_WARNING:
        mesh_info["status"] = "WARNING"
    elif min_oq is not None:
        mesh_info["status"] = "PASS"

    diagnostics = {
        "log_path": str(log_path),
        "history_headers_found": headers_found,
        "history_rows_found": len(all_rows),
        "first_header_line": first_header_line,
        "first_history_row": first_history_row,
    }

    return {"mesh": mesh_info, "all_rows": all_rows, "diagnostics": diagnostics}


def split_history_rows_by_infout(all_rows: List[Dict], cases: List[Dict]) -> Tuple[List[List[Dict]], Dict]:
    blocks = []
    idx = 0
    for case in cases:
        n = int(case.get("expected_iters") or 0)
        if n <= 0:
            blocks.append([])
            continue
        blocks.append(all_rows[idx:idx + n])
        idx += n

    split_diag = {
        "rows_consumed_by_infout": idx,
        "rows_remaining_after_split": max(0, len(all_rows) - idx),
        "total_expected_iters": sum(int(c.get("expected_iters") or 0) for c in cases),
        "blocks_found": len(blocks),
    }
    return blocks, split_diag


# =============================================================================
# METRICS AND CLASSIFICATION
# =============================================================================

def column_values(rows, column):
    return [row[column] for row in rows if column in row and row[column] is not None and math.isfinite(row[column])]


def final_window_statistics(rows, column):
    """Return mean/min/max/std and linear drift per 100 printed iterations."""
    values = column_values(rows, column)
    if not values:
        return {
            "mean": None,
            "min": None,
            "max": None,
            "std": None,
            "drift_per_100": None,
        }

    arr = np.asarray(values, dtype=float)
    drift = 0.0
    if len(arr) >= 2:
        x = np.arange(len(arr), dtype=float)
        drift = float(np.polyfit(x, arr, 1)[0] * 100.0)

    return {
        "mean": float(np.mean(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "std": float(np.std(arr, ddof=0)),
        "drift_per_100": drift,
    }


def block_metrics(rows):
    n = len(rows)
    if n == 0:
        base = {"has_history": False, "actual_iters": 0}
        for col in RESIDUAL_COLUMNS:
            base[f"{col}_final"] = None
            base[f"{col}_p95"] = None
        return base

    window_size = min(FINAL_WINDOW_MAX_ITERS, max(20, int(0.2 * n)))
    final_window = rows[-window_size:]
    final = rows[-1]

    residual_final_values = {col: abs(final[col]) for col in RESIDUAL_COLUMNS if col in final and final[col] is not None}
    residual_window_p95 = {}
    for col in RESIDUAL_COLUMNS:
        vals = [abs(v) for v in column_values(final_window, col)]
        residual_window_p95[col] = percentile(vals, 95)

    worst_residual_eq = max(residual_final_values, key=residual_final_values.get)
    worst_residual_final = residual_final_values[worst_residual_eq]
    worst_residual_p95_eq = max(
        residual_window_p95,
        key=lambda k: residual_window_p95[k] if residual_window_p95[k] is not None else -1.0,
    )
    worst_residual_p95 = residual_window_p95[worst_residual_p95_eq]

    force_ranges = {}
    for col in AERO_HISTORY_COLUMNS:
        vals = column_values(final_window, col)
        force_ranges[col] = max(vals) - min(vals) if vals else None

    cp_vals = column_values(rows, "cp-max")
    cp_window_vals = column_values(final_window, "cp-max")
    cp_final = final.get("cp-max")
    cp_window_range = max(cp_window_vals) - min(cp_window_vals) if cp_window_vals else None
    cp_spike_ratio = None
    if cp_vals:
        cp_med = median(cp_vals)
        if abs(cp_med) > 1.0e-12:
            cp_spike_ratio = max(cp_vals) / abs(cp_med)

    residual_trend_ratio = None
    if n >= 40:
        early = rows[-40:-20]
        late = rows[-20:]
        early_max = max(max(abs(r[col]) for col in RESIDUAL_COLUMNS) for r in early)
        late_max = max(max(abs(r[col]) for col in RESIDUAL_COLUMNS) for r in late)
        if early_max > 0:
            residual_trend_ratio = late_max / early_max

    final_window_stats = {
        "clzb": final_window_statistics(final_window, "clzb"),
        "cdxb": final_window_statistics(final_window, "cdxb"),
        "cmyb": final_window_statistics(final_window, "cmyb"),
        "cpmax": final_window_statistics(final_window, "cp-max"),
    }

    metrics = {
        "has_history": True,
        "actual_iters": n,
        "final_iter": int(final["iter"]),
        "window_size": window_size,
        "worst_residual_final": worst_residual_final,
        "worst_residual_eq": worst_residual_eq,
        "worst_residual_p95": worst_residual_p95,
        "worst_residual_p95_eq": worst_residual_p95_eq,
        "residual_trend_ratio": residual_trend_ratio,
        "clzb_final": final.get("clzb"),
        "cdxb_final": final.get("cdxb"),
        "cmyb_final": final.get("cmyb"),
        "cnzb_final": final.get("cnzb"),
        "cyyb_final": final.get("cyyb"),
        "crxb_final": final.get("crxb"),
        "clzb_range": force_ranges.get("clzb"),
        "cdxb_range": force_ranges.get("cdxb"),
        "cmyb_range": force_ranges.get("cmyb"),
        "cnzb_range": force_ranges.get("cnzb"),
        "cyyb_range": force_ranges.get("cyyb"),
        "crxb_range": force_ranges.get("crxb"),
        "cpmax_final": cp_final,
        "cpmax_window_range": cp_window_range,
        "cpmax_spike_ratio": cp_spike_ratio,
        "time_per_iter_final": final.get("time/iter"),
        "clzb_mean": final_window_stats["clzb"]["mean"],
        "clzb_min": final_window_stats["clzb"]["min"],
        "clzb_max": final_window_stats["clzb"]["max"],
        "clzb_std": final_window_stats["clzb"]["std"],
        "clzb_drift_per_100": final_window_stats["clzb"]["drift_per_100"],
        "cdxb_mean": final_window_stats["cdxb"]["mean"],
        "cdxb_min": final_window_stats["cdxb"]["min"],
        "cdxb_max": final_window_stats["cdxb"]["max"],
        "cdxb_std": final_window_stats["cdxb"]["std"],
        "cdxb_drift_per_100": final_window_stats["cdxb"]["drift_per_100"],
        "cmyb_mean": final_window_stats["cmyb"]["mean"],
        "cmyb_min": final_window_stats["cmyb"]["min"],
        "cmyb_max": final_window_stats["cmyb"]["max"],
        "cmyb_std": final_window_stats["cmyb"]["std"],
        "cmyb_drift_per_100": final_window_stats["cmyb"]["drift_per_100"],
        "cpmax_mean": final_window_stats["cpmax"]["mean"],
        "cpmax_min": final_window_stats["cpmax"]["min"],
        "cpmax_max": final_window_stats["cpmax"]["max"],
        "cpmax_std": final_window_stats["cpmax"]["std"],
        "cpmax_drift_per_100": final_window_stats["cpmax"]["drift_per_100"],
    }

    for col in RESIDUAL_COLUMNS:
        metrics[f"{col}_final"] = abs(final.get(col)) if final.get(col) is not None else None
        metrics[f"{col}_p95"] = residual_window_p95.get(col)

    return metrics


def assess_numerical_convergence(case, metrics):
    """Assess solver convergence independently of force/moment stability."""
    if not metrics.get("has_history"):
        return "DIVERGED", ["No iteration history found"]

    reasons = []
    worst = metrics["worst_residual_final"]
    worst_p95 = metrics["worst_residual_p95"]

    if worst >= RESIDUAL_BAD:
        reasons.append(f"final residual {worst:.3e} >= {RESIDUAL_BAD:.1e}")
    if worst_p95 is not None and worst_p95 >= RESIDUAL_BAD:
        reasons.append(f"p95 residual {worst_p95:.3e} >= {RESIDUAL_BAD:.1e}")
    if metrics.get("residual_trend_ratio") is not None and metrics["residual_trend_ratio"] > 10.0:
        reasons.append(f"residual grew by {metrics['residual_trend_ratio']:.1f}x near the end")
    if reasons:
        return "DIVERGED", reasons

    warning_reasons = []
    if worst >= RESIDUAL_ACCEPTABLE:
        warning_reasons.append(f"final residual {worst:.3e} between acceptable and bad")
    expected = case.get("expected_iters")
    actual = metrics.get("actual_iters")
    if expected is not None and actual is not None and actual < 0.8 * expected:
        warning_reasons.append(f"actual iterations {actual} < 80% of expected {expected}")
    if warning_reasons:
        return "WARNING", warning_reasons

    if worst < RESIDUAL_SUCCESSFUL:
        return "SUCCESSFUL", ["residual below successful threshold"]
    return "ACCEPTABLE", ["residual below acceptable threshold"]


def assess_aerodynamic_convergence(metrics):
    """Assess stability of monitored aero coefficients in the final window."""
    if not metrics.get("has_history"):
        return "NOT_AVAILABLE", ["No history available for aerodynamic assessment"]

    reasons = []
    if metrics.get("clzb_range") is not None and metrics["clzb_range"] > CL_WINDOW_RANGE_MAX:
        reasons.append(f"CLZB final-window range {metrics['clzb_range']:.4f}")
    if metrics.get("cdxb_range") is not None and metrics["cdxb_range"] > CD_WINDOW_RANGE_MAX:
        reasons.append(f"CDXB final-window range {metrics['cdxb_range']:.5f}")
    if metrics.get("cmyb_range") is not None and metrics["cmyb_range"] > CM_WINDOW_RANGE_MAX:
        reasons.append(f"CMYB final-window range {metrics['cmyb_range']:.5f}")
    if metrics.get("cpmax_window_range") is not None and metrics["cpmax_window_range"] > CPMAX_WINDOW_RANGE_MAX:
        reasons.append(f"cp-max final-window range {metrics['cpmax_window_range']:.3f}")
    return ("WARNING", reasons) if reasons else ("PASS", ["final-window aero monitors are stable"])


def update_overall_classification(row):
    numerical = row.get("numerical_status", "WARNING")
    aerodynamic = row.get("aerodynamic_status", "WARNING")
    neighbor = row.get("neighbor_status", "NOT_CHECKED")

    if numerical == "DIVERGED":
        overall, quality = "DIVERGED", "bad"
    elif numerical == "WARNING" or aerodynamic == "WARNING" or neighbor == "OUTLIER":
        overall, quality = "SUSPICIOUS", "warning"
    elif numerical == "ACCEPTABLE":
        overall, quality = "ACCEPTABLE", "acceptable"
    else:
        overall, quality = "CONVERGED", "successful"

    row["status"] = overall
    row["overall_status"] = overall
    row["quality"] = quality
    combined = []
    combined.extend(row.get("numerical_reasons", []))
    combined.extend(row.get("aerodynamic_reasons", []))
    combined.extend(row.get("neighbor_reasons", []))
    row["reasons"] = combined
    return row


def classify_case(case, metrics):
    """Compatibility wrapper returning the Convergence 2.0 assessments."""
    numerical_status, numerical_reasons = assess_numerical_convergence(case, metrics)
    aerodynamic_status, aerodynamic_reasons = assess_aerodynamic_convergence(metrics)
    row = {
        "numerical_status": numerical_status,
        "numerical_reasons": numerical_reasons,
        "aerodynamic_status": aerodynamic_status,
        "aerodynamic_reasons": aerodynamic_reasons,
        "neighbor_status": "NOT_CHECKED",
        "neighbor_reasons": [],
    }
    update_overall_classification(row)
    return row


def attach_adf_coefficients_to_convergence(rows, adf_data):
    """Attach stability-axis ADF coefficients to matching convergence cases.

    The outlier detector prefers CLS/CDS/CMS25 from the aerodynamic data file,
    falling back to Fluent body-axis monitor values when no exact ADF match exists.
    """
    curve_map = {
        (curve.get("case_label"), curve.get("polar")): curve.get("rows", [])
        for curve in adf_data.get("curves", [])
    }
    for row in rows:
        candidates = curve_map.get((row.get("case_label"), row.get("polar")), [])
        alpha = row.get("alpha")
        beta = row.get("beta")
        best = None
        best_distance = float("inf")
        for candidate in candidates:
            ca = candidate.get("ALPHA")
            cb = candidate.get("BETA")
            if None in (alpha, beta, ca, cb):
                continue
            distance = abs(float(alpha) - float(ca)) + abs(float(beta) - float(cb))
            if distance < best_distance:
                best_distance = distance
                best = candidate
        if best is not None and best_distance <= 1.0e-5:
            row["adf_cls"] = best.get("CLS")
            row["adf_cds"] = best.get("CDS")
            row["adf_cms25"] = best.get("CMS25")
        else:
            row["adf_cls"] = None
            row["adf_cds"] = None
            row["adf_cms25"] = None
    return rows


def _interpolate_neighbor_value(prev_r, curr_r, next_r, x_key, y_key):
    x0, x1, x2 = prev_r.get(x_key), curr_r.get(x_key), next_r.get(x_key)
    y0, y2 = prev_r.get(y_key), next_r.get(y_key)
    if None in (x0, x1, x2, y0, y2) or abs(x2 - x0) < 1.0e-12:
        return None
    return y0 + (x1 - x0) * (y2 - y0) / (x2 - x0)


def add_neighbor_consistency(rows):
    """Flag local polar outliers using linear interpolation between adjacent points."""
    grouped = {}
    for row in rows:
        key = (row["case_label"], row["polar"])
        grouped.setdefault(key, []).append(row)

    for group in grouped.values():
        alpha_values = [r.get("alpha") for r in group if r.get("alpha") is not None]
        beta_values = [r.get("beta") for r in group if r.get("beta") is not None]
        alpha_span = max(alpha_values) - min(alpha_values) if alpha_values else 0.0
        beta_span = max(beta_values) - min(beta_values) if beta_values else 0.0
        sweep_key = "alpha" if alpha_span >= beta_span else "beta"

        valid = [
            r for r in group
            if r.get("numerical_status") != "DIVERGED"
            and r.get(sweep_key) is not None
        ]
        valid = sorted(valid, key=lambda r: r[sweep_key])

        use_adf = sum(
            1 for r in valid
            if r.get("adf_cls") is not None and r.get("adf_cds") is not None and r.get("adf_cms25") is not None
        ) >= 3
        cl_key, cd_key, cm_key = (
            ("adf_cls", "adf_cds", "adf_cms25")
            if use_adf else
            ("clzb_final", "cdxb_final", "cmyb_final")
        )
        source_name = "ADF stability axes (CLS/CDS/CMS25)" if use_adf else "Fluent monitor body axes (CLZB/CDXB/CMYB)"

        for row in group:
            row["neighbor_status"] = "NOT_CHECKED"
            row["neighbor_reasons"] = []
            row["outlier_score"] = None
            row["neighbor_dcl"] = None
            row["neighbor_dcd"] = None
            row["neighbor_dcm"] = None
            row["neighbor_sweep"] = sweep_key.upper()
            row["outlier_source"] = source_name

        for i in range(1, len(valid) - 1):
            prev_r, curr, next_r = valid[i - 1], valid[i], valid[i + 1]
            expected_cl = _interpolate_neighbor_value(prev_r, curr, next_r, sweep_key, cl_key)
            expected_cd = _interpolate_neighbor_value(prev_r, curr, next_r, sweep_key, cd_key)
            expected_cm = _interpolate_neighbor_value(prev_r, curr, next_r, sweep_key, cm_key)

            dcl = abs(curr[cl_key] - expected_cl) if expected_cl is not None and curr.get(cl_key) is not None else None
            dcd = abs(curr[cd_key] - expected_cd) if expected_cd is not None and curr.get(cd_key) is not None else None
            dcm = abs(curr[cm_key] - expected_cm) if expected_cm is not None and curr.get(cm_key) is not None else None
            curr["neighbor_dcl"], curr["neighbor_dcd"], curr["neighbor_dcm"] = dcl, dcd, dcm

            scores = []
            reasons = []
            if dcl is not None:
                scores.append(dcl / NEIGHBOR_CL_JUMP_MAX)
                if dcl > NEIGHBOR_CL_JUMP_MAX:
                    reasons.append(f"neighbor CL interpolation error {dcl:.4f}")
            if dcd is not None:
                scores.append(dcd / NEIGHBOR_CD_JUMP_MAX)
                if dcd > NEIGHBOR_CD_JUMP_MAX:
                    reasons.append(f"neighbor CD interpolation error {dcd:.4f}")
            if dcm is not None:
                scores.append(dcm / NEIGHBOR_CM_JUMP_MAX)
                if dcm > NEIGHBOR_CM_JUMP_MAX:
                    reasons.append(f"neighbor CM interpolation error {dcm:.4f}")

            curr["outlier_score"] = max(scores) if scores else None
            curr["neighbor_status"] = "OUTLIER" if reasons else "PASS"
            curr["neighbor_reasons"] = reasons or ["consistent with adjacent sweep points"]

    for row in rows:
        update_overall_classification(row)
    return rows


# =============================================================================
# POLAR CONVERGENCE PROCESSING
# =============================================================================

def process_polar_convergence(case_label, polar_name, polar_path):
    infout_path = polar_path / "infout"
    log_path = find_fluent_log(polar_path)

    if not infout_path.exists():
        raise FileNotFoundError(f"Missing infout: {infout_path}")
    if not log_path.exists():
        raise FileNotFoundError(f"Missing FLUENT_LOG: {log_path}")

    infout = parse_infout(infout_path)
    log = parse_fluent_log_all_rows(log_path)
    cases = infout["cases"]
    blocks, split_diag = split_history_rows_by_infout(log["all_rows"], cases)

    diag = dict(log["diagnostics"])
    diag.update(split_diag)

    print(f"  FLUENT_LOG: {log_path.name}")
    print(f"  Monitor headers found      : {diag['history_headers_found']}")
    print(f"  Monitor rows found         : {diag['history_rows_found']}")
    print(f"  Total expected iters infout: {diag['total_expected_iters']}")
    print(f"  Rows consumed by infout    : {diag['rows_consumed_by_infout']}")
    print(f"  Rows remaining after split : {diag['rows_remaining_after_split']}")

    rows = []
    history = {}

    for i, case in enumerate(cases):
        hist_rows = blocks[i] if i < len(blocks) else []
        metrics = block_metrics(hist_rows)
        assessment = classify_case(case, metrics)

        row = {}
        row.update(case)
        row.update(metrics)
        row.update(assessment)
        row["case_label"] = case_label
        row["polar"] = polar_name
        row["polar_path"] = str(polar_path)
        row["case_key"] = f"{case_label}|{polar_name}|{case['case']}"
        rows.append(row)

        history[row["case_key"]] = downsample_history(hist_rows)

    summary = {
        "case_label": case_label,
        "polar": polar_name,
        "polar_path": str(polar_path),
        "meta": infout["meta"],
        "deflections": infout.get("deflections", {}),
        "mesh": log["mesh"],
        "diagnostics": diag,
        "n_cases_infout": len(cases),
        "n_history_blocks": len(blocks),
        "files": {
            "infout": str(infout_path),
            "infout_modified": datetime.fromtimestamp(infout_path.stat().st_mtime, tz=timezone.utc).isoformat() if infout_path.exists() else None,
            "fluent_log": str(log_path),
            "fluent_log_modified": datetime.fromtimestamp(log_path.stat().st_mtime, tz=timezone.utc).isoformat() if log_path.exists() else None,
        },
    }
    return summary, rows, history


def process_optional_convergence(case_label, polar_name, polar_path):
    """Keep valid ADF data usable when optional run diagnostics are unavailable."""
    try:
        return process_polar_convergence(case_label, polar_name, polar_path)
    except (OSError, ValueError, KeyError, IndexError) as error:
        message = f"Convergence unavailable: {error}"
        print(f"  WARNING: {message}")
        info = {"meta": {}, "cases": [], "deflections": {}}
        infout_path = polar_path / 'infout'
        try:
            info = parse_infout(infout_path)
        except (OSError, ValueError, KeyError, IndexError):
            pass
        log_path = find_fluent_log(polar_path)
        summary = {
            "case_label": case_label, "polar": polar_name, "polar_path": str(polar_path),
            "meta": info['meta'], "deflections": info.get('deflections', {}),
            "mesh": {"status": "NOT_AVAILABLE", "warnings": []}, "diagnostics": {},
            "n_cases_infout": len(info['cases']), "n_history_blocks": 0,
            "convergence_unavailable": message,
            "files": {"infout": str(infout_path), "fluent_log": str(log_path)},
        }
        return summary, [], {}


def build_provenance(case_configs, summaries, adf_data, drag_rise_data):
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script_version": SCRIPT_VERSION,
        "module_versions": MODULE_VERSIONS,
        "configurations": [
            {
                "label": cfg.label,
                "adf_directory": str(cfg.directory),
                "polars": [polar_name_from_number(p) for p in cfg.polars],
                "drag_rise_dir": cfg.drag_rise_dir,
            }
            for cfg in case_configs
        ],
        "adf_files": [
            {
                "configuration": c.get("case_label"),
                "polar": c.get("polar"),
                "path": c.get("path"),
                "modified": c.get("modified"),
            }
            for c in adf_data.get("curves", [])
        ],
        "run_files": [
            {
                "configuration": s.get("case_label"),
                "polar": s.get("polar"),
                **s.get("files", {}),
            }
            for s in summaries
        ],
        "drag_rise_files": [
            {
                "configuration": c.get("case_label"),
                "file_name": c.get("file_name"),
                "path": c.get("path"),
                "modified": c.get("modified"),
            }
            for c in drag_rise_data.get("curves", [])
        ],
    }


def build_integrity_checks(case_configs, summaries, conv_rows, adf_data, drag_rise_data):
    checks = []
    requested = {
        (cfg.label, polar_name_from_number(p))
        for cfg in case_configs
        for p in cfg.polars
    }
    loaded_adf = {(c.get("case_label"), c.get("polar")) for c in adf_data.get("curves", [])}
    loaded_runs = {(s.get("case_label"), s.get("polar")) for s in summaries}

    for label, polar in sorted(requested):
        if (label, polar) not in loaded_adf:
            checks.append({"severity": "ERROR", "configuration": label, "polar": polar, "check": "ADF file", "details": "Requested ADF file was not loaded."})
        if (label, polar) not in loaded_runs:
            checks.append({"severity": "ERROR", "configuration": label, "polar": polar, "check": "Run files", "details": "infout/FLUENT_LOG data was not processed."})

    summary_map = {(s.get("case_label"), s.get("polar")): s for s in summaries}
    for curve in adf_data.get("curves", []):
        label, polar = curve.get("case_label"), curve.get("polar")
        rows = curve.get("rows", [])
        sweep = curve.get("sweep_var", "ALPHA")
        sweep_values = [r.get(sweep) for r in rows if r.get(sweep) is not None]
        if len(sweep_values) != len(set(sweep_values)):
            checks.append({"severity": "WARNING", "configuration": label, "polar": polar, "check": "Duplicate sweep points", "details": f"Duplicate {sweep} values are present in the ADF table."})
        if len(sweep_values) >= 3:
            diffs = np.diff(np.asarray(sweep_values, dtype=float))
            if not (np.all(diffs >= 0) or np.all(diffs <= 0)):
                checks.append({"severity": "WARNING", "configuration": label, "polar": polar, "check": "Sweep order", "details": f"{sweep} is non-monotonic in file order."})
        mach = [r.get("MACH") for r in rows if r.get("MACH") is not None]
        reynolds = [r.get("REYNOLDS") for r in rows if r.get("REYNOLDS") is not None]
        if mach and max(mach) - min(mach) > 1.0e-6:
            checks.append({"severity": "INFO", "configuration": label, "polar": polar, "check": "Mach variation", "details": f"Mach varies from {min(mach):.5g} to {max(mach):.5g}."})
        if reynolds and max(reynolds) - min(reynolds) > max(1.0, 1.0e-6 * abs(np.mean(reynolds))):
            checks.append({"severity": "INFO", "configuration": label, "polar": polar, "check": "Reynolds variation", "details": f"Reynolds varies from {min(reynolds):.5g} to {max(reynolds):.5g}."})
        summary = summary_map.get((label, polar))
        if summary and summary.get("n_cases_infout") != len(rows):
            checks.append({"severity": "WARNING", "configuration": label, "polar": polar, "check": "Case count mismatch", "details": f"infout has {summary.get('n_cases_infout')} cases while ADF has {len(rows)} rows."})

    for summary in summaries:
        if summary.get('convergence_unavailable'):
            checks.append({"severity": "WARNING", "configuration": summary['case_label'],
                           "polar": summary['polar'], "check": "Convergence unavailable",
                           "details": summary['convergence_unavailable']})
        d = summary.get("diagnostics", {})
        if d.get("history_rows_found", 0) < d.get("total_expected_iters", 0):
            checks.append({"severity": "ERROR", "configuration": summary.get("case_label"), "polar": summary.get("polar"), "check": "Incomplete Fluent history", "details": f"Found {d.get('history_rows_found', 0)} monitor rows; expected {d.get('total_expected_iters', 0)}."})
        if d.get("rows_remaining_after_split", 0) > 0:
            checks.append({"severity": "INFO", "configuration": summary.get("case_label"), "polar": summary.get("polar"), "check": "Extra Fluent rows", "details": f"{d.get('rows_remaining_after_split')} monitor rows remain after infout splitting."})

    for curve in drag_rise_data.get("curves", []):
        target = curve.get("cls_value")
        actual = [r.get("CLS") for r in curve.get("rows", []) if r.get("CLS") is not None]
        if target is not None and actual:
            max_error = max(abs(float(v) - float(target)) for v in actual)
            if max_error > 0.02:
                checks.append({"severity": "WARNING", "configuration": curve.get("case_label"), "polar": curve.get("file_name"), "check": "Drag-rise CLS mismatch", "details": f"Filename target CLS={target:+.3f}; maximum achieved deviation is {max_error:.4f}."})

    if not checks:
        checks.append({"severity": "PASS", "configuration": "All", "polar": "All", "check": "Data integrity", "details": "No integrity problems detected."})
    return checks


# =============================================================================
# OUTPUT WRITERS
# =============================================================================

def write_csv(path, rows):
    fieldnames = [
        "case_label", "polar", "case", "name", "mach", "reynolds", "alpha", "beta",
        "expected_iters", "actual_iters", "status", "overall_status", "quality",
        "numerical_status", "aerodynamic_status", "neighbor_status", "outlier_score",
        "neighbor_dcl", "neighbor_dcd", "neighbor_dcm", "neighbor_sweep", "outlier_source",
        "adf_cls", "adf_cds", "adf_cms25",
        "continuity_final", "x-velocity_final", "y-velocity_final", "z-velocity_final", "energy_final", "nut_final",
        "worst_residual_final", "worst_residual_eq", "worst_residual_p95", "worst_residual_p95_eq", "residual_trend_ratio",
        "clzb_final", "cdxb_final", "cmyb_final", "cnzb_final", "cyyb_final", "crxb_final",
        "clzb_range", "cdxb_range", "cmyb_range", "cpmax_final", "cpmax_window_range", "cpmax_spike_ratio",
        "clzb_mean", "clzb_min", "clzb_max", "clzb_std", "clzb_drift_per_100",
        "cdxb_mean", "cdxb_min", "cdxb_max", "cdxb_std", "cdxb_drift_per_100",
        "cmyb_mean", "cmyb_min", "cmyb_max", "cmyb_std", "cmyb_drift_per_100",
        "cpmax_mean", "cpmax_min", "cpmax_max", "cpmax_std", "cpmax_drift_per_100",
        "case_key", "polar_path", "numerical_reasons", "aerodynamic_reasons", "neighbor_reasons", "reasons",
    ]
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in ("numerical_reasons", "aerodynamic_reasons", "neighbor_reasons", "reasons"):
                out[key] = " | ".join(row.get(key, []))
            writer.writerow(out)


def write_json(path, summaries, rows, history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data=None):
    payload = {
        "summaries": summaries,
        "results": rows,
        "history": history,
        "adf": adf_data,
        "drag_rise": drag_rise_data,
        "script_version": SCRIPT_VERSION,
        "module_versions": MODULE_VERSIONS,
        "provenance": provenance,
        "integrity_checks": integrity_checks,
        "distributions": distribution_data or {"series": [], "issues": []},
    }
    Path(path).write_text(json.dumps(json_safe(payload), indent=2), encoding="utf-8")


def write_static_margin_csv(path, sm_df):
    if sm_df is None or sm_df.empty:
        return
    sm_df.to_csv(path, index=False)


# =============================================================================
# HTML DASHBOARD
# =============================================================================

def make_html(summaries, conv_rows, history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data=None):
    summaries_json = json.dumps(json_safe(summaries))
    conv_json = json.dumps(json_safe(conv_rows))
    history_json = json.dumps(json_safe(history))
    adf_json = json.dumps(json_safe(adf_data))
    drag_json = json.dumps(json_safe(drag_rise_data))
    module_versions_json = json.dumps(json_safe(MODULE_VERSIONS))
    script_version_json = json.dumps(SCRIPT_VERSION)
    provenance_json = json.dumps(json_safe(provenance))
    integrity_json = json.dumps(json_safe(integrity_checks))
    sm_ylim_json = json.dumps(STATIC_MARGIN_YLIM)
    filter_sm_json = json.dumps(bool(FILTER_STATIC_MARGIN_OUTLIERS_FOR_PLOT))

    # JSON is embedded in script context; escape literal script terminators.
    summaries_json = summaries_json.replace('<', '\\u003c')
    conv_json = conv_json.replace('<', '\\u003c')
    history_json = history_json.replace('<', '\\u003c')
    adf_json = adf_json.replace('<', '\\u003c')
    drag_json = drag_json.replace('<', '\\u003c')
    provenance_json = provenance_json.replace('<', '\\u003c')
    integrity_json = integrity_json.replace('<', '\\u003c')
    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>JAMAL · CFD post-processing | Aerodynamic results dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>

:root {{
  --font-words: ui-serif, Georgia, Cambria, "Times New Roman", serif;
  --font-numbers: Arial, Helvetica, sans-serif;
}}
body {{ font-family: Arial, Helvetica, sans-serif; margin: 24px; background: #f6f7f9; color: #222; }}
h1, h2, h3 {{ margin-bottom: 8px; font-family: var(--font-words); font-weight: 600; letter-spacing: -0.015em; }}
.dashboard-hero {{ margin: 2px 0 18px; padding: 4px 2px 6px; }}
.dashboard-brand {{ font-family: var(--font-words); font-size: 15px; font-weight: 600; color: #46515f; letter-spacing: 0.015em; }}
.dashboard-hero h1 {{ margin: 5px 0 4px; font-size: 32px; line-height: 1.1; }}
.dashboard-tagline {{ margin: 0; font-family: var(--font-words); font-size: 15px; color: #5f6873; line-height: 1.45; }}
.card {{ background: white; border-radius: 12px; padding: 18px; margin-bottom: 18px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
.kpi {{ display: inline-block; padding: 12px 18px; border-radius: 10px; background: #f0f2f5; margin-right: 10px; margin-bottom: 10px; }}
.kpi > div {{ font-family: var(--font-words); }}
.kpi strong {{ font-family: var(--font-numbers); font-size: 24px; }}
.status-CONVERGED {{ color: #0a7a28; font-weight: bold; }}
.status-ACCEPTABLE {{ color: #1769aa; font-weight: bold; }}
.status-SUSPICIOUS {{ color: #b26a00; font-weight: bold; }}
.status-DIVERGED {{ color: #b00020; font-weight: bold; }}
.assessment-PASS, .assessment-SUCCESSFUL {{ color: #0a7a28; font-weight: bold; }}
.assessment-ACCEPTABLE {{ color: #1769aa; font-weight: bold; }}
.assessment-WARNING, .assessment-OUTLIER {{ color: #b26a00; font-weight: bold; }}
.assessment-DIVERGED, .assessment-CRITICAL {{ color: #b00020; font-weight: bold; }}
.assessment-NOT_CHECKED, .assessment-NOT_AVAILABLE, .assessment-UNKNOWN {{ color: #666; font-weight: bold; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th, td {{ padding: 7px; border-bottom: 1px solid #ddd; text-align: left; }}
th {{ background: #f0f2f5; font-family: var(--font-words); }}
tr:hover {{ background: #fafafa; }}
.controls {{ margin-bottom: 12px; }}
select, input {{ font-family: var(--font-numbers); padding: 8px 12px; border: 1px solid #ccc; border-radius: 8px; background: white; cursor: pointer; margin-right: 8px; margin-bottom: 8px; }}
button {{ font-family: var(--font-words); padding: 8px 12px; border: 1px solid #ccc; border-radius: 8px; background: white; cursor: pointer; margin-right: 8px; margin-bottom: 8px; }}
label {{ font-family: var(--font-words); }}
.plot {{ height: 440px; }}
.analysis-plot-grid {{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;align-items:start;}}
.analysis-plot-grid > * {{min-width:0;}}
.analysis-plot-grid .plot {{width:100%;min-width:0;}}
#coefficientPlots {{display:block;}}
.coefficient-row {{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px;}}
.coefficient-row > * {{min-width:0;}}
.coefficient-row h2 {{font-size:17px;}}
.numeric-value {{font-family:var(--font-numbers);}}
@media(max-width:760px) {{.coefficient-row {{grid-template-columns:minmax(0,1fr);}}}}
@media(max-width:760px) {{.analysis-plot-grid {{grid-template-columns:minmax(0,1fr);}}}}
.small {{ color: #666; font-size: 12px; }}
.badge {{ display:inline-block; padding: 6px 9px; margin: 3px; border-radius: 7px; color: white; font-size: 12px; }}
.dashboard-nav {{ position: sticky; top: 0; z-index: 20; background: #f6f7f9; padding: 10px 0; margin-bottom: 10px; border-bottom: 1px solid #ddd; }}
.dashboard-nav button {{ font-family: var(--font-words); font-weight: 600; }}
.dashboard-nav button.active {{ background: #2f5f98; color: white; border-color: #2f5f98; }}
.dashboard-section {{ display: none; }}
.dashboard-section.active {{ display: block; }}
.table-scroll {{ overflow-x: auto; }}
.metric-note {{ display: inline-block; margin: 3px 8px 3px 0; padding: 5px 8px; border-radius: 6px; background: #f0f2f5; font-family: var(--font-words); font-size: 12px; }}
.metric-note strong {{ font-family: var(--font-numbers); }}
.analysis-toolbar {{ position: sticky; top: 58px; z-index: 18; border: 1px solid #d7dce2; box-shadow: 0 3px 12px rgba(0,0,0,0.10); }}
.toolbar-row {{ display: flex; flex-wrap: wrap; align-items: center; gap: 6px 10px; }}
.toolbar-separator {{ width: 1px; height: 30px; background: #d6dbe1; margin: 0 4px; }}
.advanced-panel {{ display: none; border-left: 4px solid #2f5f98; }}
.advanced-panel.open {{ display: block; }}
.control-group {{ padding: 10px 0; border-bottom: 1px solid #e6e9ed; }}
.control-group:last-child {{ border-bottom: none; }}
.control-group h3 {{ margin: 0 0 8px; font-size: 16px; }}
.status-panel {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 10px; }}
.status-box {{ background: #f5f7f9; border-radius: 9px; padding: 10px; border: 1px solid #e2e6ea; }}
.status-box .value {{ font-family: var(--font-numbers); font-size: 17px; font-weight: 700; }}
.status-box .label {{ font-family: var(--font-words); color: #58616b; font-size: 12px; }}
.focus-label {{ font-family: var(--font-words); font-size: 12px; color: #46515f; }}
.validation-card {{ border-left: 5px solid #0a7a28; }}
.integrity-ERROR {{ color: #b00020; font-weight: bold; }}
.integrity-WARNING {{ color: #b26a00; font-weight: bold; }}
.integrity-INFO {{ color: #1769aa; font-weight: bold; }}
.integrity-PASS {{ color: #0a7a28; font-weight: bold; }}
.selected-row {{ background: #edf4fb !important; }}
body.density-compact {{ margin: 14px; }}
body.density-compact .card {{ padding: 11px; margin-bottom: 11px; }}
body.density-compact .plot {{ height: 340px; }}
body.density-compact th, body.density-compact td {{ padding: 4px 6px; font-size: 11px; }}
body.density-comfortable .plot {{ height: 440px; }}
body.density-presentation {{ margin: 30px; }}
body.density-presentation .card {{ padding: 24px; margin-bottom: 24px; }}
body.density-presentation .plot {{ height: 560px; }}
body.density-presentation th, body.density-presentation td {{ padding: 9px; font-size: 14px; }}
.comparison-buttons button.active, .metric-buttons button.active {{ background: #2f5f98; color: white; border-color: #2f5f98; }}
.comparison-pair {{ display: flex; flex-wrap: wrap; align-items: center; gap: 6px 10px; margin: 8px 0; padding: 10px; border: 1px solid #dfe4ea; border-radius: 8px; background: #f8f9fb; }}
.comparison-pair .pair-index {{ font-family: var(--font-numbers); min-width: 24px; font-weight: 700; color: #2f5f98; }}
.comparison-pair select {{ min-width: 250px; max-width: 430px; }}
.curve-style-custom {{ display: none; }}
.curve-style-custom.open {{ display: inline-flex; flex-wrap: wrap; align-items: center; gap: 4px; }}
.style-help {{ display: block; margin-top: 7px; color: #666; font-size: 12px; }}
@media (max-width: 900px) {{
  .analysis-toolbar {{ position: static; }}
  .toolbar-separator {{ display: none; }}
  body {{ margin: 12px; }}
}}
.dashboard-hero {{margin-bottom:8px;}}
.dashboard-hero h1 {{font-size:26px;}}
.dashboard-tagline {{display:none;}}
.workspace-details {{padding:12px 16px;}}
.workspace-details summary,.toolbar-options summary {{cursor:pointer;font-family:var(--font-words);font-weight:600;}}
.workspace-details[open] summary {{margin-bottom:14px;}}
#summaryPanel {{margin-bottom:10px;}}
#compactStatus {{margin-left:16px;font-family:var(--font-numbers);font-weight:400;}}
.analysis-toolbar {{padding:10px 14px;margin-bottom:12px;}}
.toolbar-row select,.toolbar-row button {{margin-bottom:0;}}
.toolbar-options {{position:relative;}}
.toolbar-options > .toolbar-row {{position:absolute;right:0;top:30px;width:min(480px,80vw);padding:16px;background:white;box-shadow:0 8px 30px #0002;border:1px solid #d7dce2;border-radius:10px;}}
.coefficient-controls {{padding:12px 16px;margin-bottom:10px;}}
.coefficient-controls h2,.coefficient-controls p {{display:none;}}
.coefficient-controls .controls {{margin:0;}}
.coefficient-controls select {{margin-bottom:0;}}
.workspace-context {{font-family:var(--font-words);color:#46515f;padding:8px 2px 12px;}}
.shared-legend {{display:flex;flex-wrap:wrap;gap:6px;margin:0 0 10px;}}
.shared-legend button {{font-family:var(--font-numbers);font-size:12px;margin:0;}}
.shared-legend button[aria-pressed="false"] {{opacity:.45;text-decoration:line-through;}}
.legend-swatch {{display:inline-block;width:26px;margin-right:8px;vertical-align:middle;border-top:3px solid;}}
.plot-actions {{display:flex;justify-content:flex-end;gap:4px;margin:0 0 2px;}}
.plot-actions button {{padding:4px 8px;font-size:12px;margin:0;}}
#coefficientPlots .card {{padding:12px;margin-bottom:10px;}}
#coefficientPlots h2 {{display:none;}}
#plotDialog {{width:90vw;max-width:1600px;border:1px solid #d7dce2;border-radius:14px;padding:18px;}}
#plotDialog::backdrop {{background:rgba(25,35,50,.55);}}
#plotDialog .plot {{height:72vh;}}
#plotDialog header {{display:flex;justify-content:space-between;align-items:center;}}
#plotDialog h2 {{margin:0;}}
#referencePanel .controls,#momentRefInfo {{overflow-x:auto;}}
@media(max-width:900px) {{.toolbar-options > .toolbar-row {{position:static;width:auto;}}}}
</style>
</head>
<body>
<header class="dashboard-hero">
  <div class="dashboard-brand">JAMAL · CFD post-processing</div>
  <h1>Aerodynamic results dashboard</h1>
  <p class="dashboard-tagline">Interactive aerodynamic, convergence, drag-rise, and reference-moment analysis in a portable standalone report.</p>
</header>
<div class="dashboard-nav" id="dashboardNav">
  <button class="active" data-section="overview" onclick="showSection('overview', this)">Overview</button>
  <button data-section="aero" onclick="showSection('aero', this)">Coefficients</button>
  <button data-section="static-margin" onclick="showSection('static-margin', this)">Static margin</button>
  <button data-section="comparison" onclick="showSection('comparison', this)">Delta</button>
  <button data-section="drag" onclick="showSection('drag', this)">Drag rise</button>
  <button data-section="convergence" onclick="showSection('convergence', this)">Convergence</button>
  <button data-section="quality" onclick="showSection('quality', this)">Quality</button>
  <button data-section="tables" onclick="showSection('tables', this)">Tables</button>
</div>

<details class="card workspace-details" id="summaryPanel">
  <summary>Run summary <span id="compactStatus" class="small"></span></summary>
  <div id="kpis"></div>
  <p class="small">Aerodynamic data files are read from 03-RESULTS/ADF. Convergence histories are read from 02-RUNS/POLAR-XXX/infout and FLUENT_LOG.</p>
</details>

<div class="card analysis-toolbar" id="analysisToolbar">
  <div class="toolbar-row">
    <label>Configuration:</label><select id="labelFilter" onchange="refreshAll()"></select>
    <label>POLAR:</label><select id="polarFilter" onchange="refreshAll()"></select>
    <label>Status:</label><select id="statusFilter" onchange="refreshAll()">
      <option value="ALL">ALL</option><option value="CONVERGED">CONVERGED</option><option value="ACCEPTABLE">ACCEPTABLE</option><option value="SUSPICIOUS">SUSPICIOUS</option><option value="DIVERGED">DIVERGED</option>
    </select>
    <details class="toolbar-options"><summary>View options</summary><div class="toolbar-row">
    <button onclick="setFocusMode('all')">Show all</button>
    <button onclick="setFocusMode('baseline')">Baseline only</button>
    <button onclick="setFocusMode('selected')">Focus selected</button>
    <button onclick="setFocusMode('selected_baseline')">Selected + baseline</button>
    <button onclick="setFocusMode('none')">Hide all</button>
    <button onclick="showSection('convergence', document.querySelector('#dashboardNav button[data-section=&quot;convergence&quot;]'))">Selected history</button>
    <span id="focusModeLabel" class="focus-label"></span>
    <span class="toolbar-separator"></span>
    <label>Density:</label><select id="densitySelect" onchange="setDensity(this.value)">
      <option value="compact">Compact</option><option value="comfortable" selected>Comfortable</option><option value="presentation">Presentation</option>
    </select>
    </div></details>
    <button onclick="toggleAdvancedPanel()">Controls &amp; export</button>
  </div>
</div>
<div class="card advanced-panel" id="advancedPanel">
  <div class="control-group">
    <h3>Classification thresholds</h3>
    <label>Residual successful:</label><input id="thrResSuccessful" type="number" step="0.0001" value="0.0008" onchange="refreshAll()">
    <label>Residual acceptable:</label><input id="thrResAcceptable" type="number" step="0.0001" value="0.0080" onchange="refreshAll()">
    <label>Residual bad:</label><input id="thrResBad" type="number" step="0.0001" value="0.0100" onchange="refreshAll()">
    <label>CL range:</label><input id="thrClRange" type="number" step="0.0001" value="0.0030" onchange="refreshAll()">
    <label>CD range:</label><input id="thrCdRange" type="number" step="0.0001" value="0.0005" onchange="refreshAll()">
    <label>CM range:</label><input id="thrCmRange" type="number" step="0.0001" value="0.0010" onchange="refreshAll()">
    <label>cp-max range:</label><input id="thrCpRange" type="number" step="0.01" value="0.25" onchange="refreshAll()">
    <label>Neighbor ΔCL:</label><input id="thrNeighborCl" type="number" step="0.01" value="0.12" onchange="refreshAll()">
    <label>Neighbor ΔCD:</label><input id="thrNeighborCd" type="number" step="0.001" value="0.03" onchange="refreshAll()">
    <label>Neighbor ΔCM:</label><input id="thrNeighborCm" type="number" step="0.001" value="0.05" onchange="refreshAll()">
    <button onclick="resetThresholds()">Reset JAMAL defaults</button>
  </div>
  <div class="control-group">
    <h3>Curve appearance</h3>
    <label>Style mode:</label>
    <select id="curveStyleMode" onchange="updateCurveStyleControls(); refreshAll()">
      <option value="auto" selected>Auto</option>
      <option value="colors_only">Colors only</option>
      <option value="color_polar">Color by POLAR</option>
      <option value="color_configuration">Color by configuration</option>
      <option value="custom">Custom</option>
    </select>
    <span id="curveStyleCustom" class="curve-style-custom">
      <label>Color:</label><select id="customColorMapping" onchange="refreshAll()">
        <option value="every_curve">Every curve</option>
        <option value="polar">POLAR</option>
        <option value="configuration">Configuration</option>
      </select>
      <label>Line:</label><select id="customLineMapping" onchange="refreshAll()">
        <option value="solid">Solid only</option>
        <option value="configuration">Configuration</option>
        <option value="polar">POLAR</option>
      </select>
      <label>Marker:</label><select id="customMarkerMapping" onchange="refreshAll()">
        <option value="none">None</option>
        <option value="configuration">Configuration</option>
        <option value="polar">POLAR</option>
        <option value="every_curve">Every curve</option>
      </select>
    </span>
    <span class="style-help">Auto uses solid, differently colored POLAR curves for one configuration. With multiple configurations, matching POLARs share a color and configurations use different line styles.</span>
  </div>
  <div class="control-group">
    <h3>Views and presets</h3>
    <button onclick="applyBuiltInPreset('longitudinal')">Longitudinal review</button>
    <button onclick="applyBuiltInPreset('lateral')">Lateral-directional review</button>
    <button onclick="applyBuiltInPreset('convergence')">Convergence audit</button>
    <button onclick="applyBuiltInPreset('drag')">Drag-rise comparison</button><br>
    <label>Preset name:</label><input id="presetNameInput" type="text" placeholder="My review">
    <button onclick="saveNamedPreset()">Save preset</button>
    <select id="presetSelect"></select>
    <button onclick="loadNamedPreset()">Load</button>
    <button onclick="deleteNamedPreset()">Delete</button>
  </div>
  <div class="control-group">
    <h3>Export current view</h3>
    <label>Plot:</label><select id="exportPlotSelect"></select>
    <button onclick="exportSelectedPlot('png')">Export PNG</button>
    <button onclick="exportSelectedPlot('svg')">Export SVG</button>
    <button onclick="exportFilteredCsv()">Export filtered CSV</button>
    <button onclick="copySelectedCaseSummary()">Copy selected-case summary</button>
  </div>
</div>

<section id="section-overview" class="dashboard-section active">
<div class="card"><h2>Status map</h2><div id="statusMap"></div></div>
<div class="card"><h2>Convergence classification 2.0</h2>
  <p class="small">Numerical convergence, aerodynamic final-window stability, and neighbor consistency are assessed independently. Overall status combines the three assessments.</p>
  <div id="assessmentKpis"></div>
</div>
<div class="card"><h2>Module versions</h2>
  <p class="small">Validated modules are versioned independently from the dashboard script.</p>
  <table id="moduleVersionTable"><thead><tr><th>Module</th><th>Version</th></tr></thead><tbody></tbody></table>
</div>
</section>

<section id="section-aero" class="dashboard-section">
<details class="card workspace-details" id="referencePanel"><summary>Moment reference <span id="referenceSummary"></span> · Edit</summary>
  <p class="small">Default mode applies the selected %MAC independently to each POLAR using its own XREF and CREF from infout; Y and Z shifts are added to each POLAR's own YREF/ZREF. Use absolute XYZ only when you want every POLAR shifted to the same aircraft coordinate.</p>
  <div class="controls">
    <label>Reference mode:</label>
    <select id="refModeSelect" onchange="refreshAll()">
      <option value="percent" selected>% MAC per POLAR</option>
      <option value="absolute">Absolute XYZ</option>
    </select>
    <label>X reference [% MAC]:</label><input id="xrefPercentInput" type="number" value="25.0" step="0.1" onchange="refreshAll()" oninput="refreshAll()">
    <label>Y shift [m] (%MAC mode):</label><input id="yOffsetInput" type="number" value="0.0" step="0.001" onchange="refreshAll()" oninput="refreshAll()">
    <label>Z shift [m] (%MAC mode):</label><input id="zOffsetInput" type="number" value="0.0" step="0.001" onchange="refreshAll()" oninput="refreshAll()">
    <label>X absolute [m]:</label><input id="xAbsInput" type="number" value="0.0" step="0.001" onchange="refreshAll()" oninput="refreshAll()">
    <label>Y absolute [m]:</label><input id="yAbsInput" type="number" value="0.0" step="0.001" onchange="refreshAll()" oninput="refreshAll()">
    <label>Z absolute [m]:</label><input id="zAbsInput" type="number" value="0.0" step="0.001" onchange="refreshAll()" oninput="refreshAll()">
    <button onclick="resetMomentReference()">Reset to 25% MAC</button>
  </div>
  <div id="momentRefInfo" class="small"></div>
</details>

<details class="card workspace-details" id="deflectionPanel"><summary>Control-surface and flap deflections</summary><div class="table-scroll"><table id="deflectionTable"><thead><tr><th>Configuration</th><th>POLAR</th><th>RUD1</th><th>RUD2</th><th>RUD3</th><th>RUD4</th><th>ELV1</th><th>ELV2</th><th>ELV3</th><th>ELV4</th><th>AIL1</th><th>AIL2</th><th>AIL3</th><th>AIL4</th><th>FLP1</th><th>FLP2</th><th>FLP3</th><th>FLP4</th></tr></thead><tbody></tbody></table></div></details>

<div class="card coefficient-controls"><h2 id="aeroCoeffHeading">Aerodynamic coefficients vs ALPHA</h2>
  <div class="controls">
    <label>Axis system:</label>
    <select id="stabilityAxisSelect" onchange="drawAdfCoeffPlot()">
      <option value="B">Body axis (B)</option>
      <option value="S" selected>Stability axis (S)</option>
      <option value="W">Wind axis (W)</option>
    </select>
    <label>Moment reference:</label>
    <select id="momentReferenceSelect" onchange="refreshAll()">
      <option value="original" selected>Original ADF reference</option>
      <option value="user">User reference point</option>
    </select>
    <label for="coefficientAbscissa">Plot against:</label>
    <select id="coefficientAbscissa" onchange="drawAdfCoeffPlot();saveLastState();">
      <option value="ALPHA">ALPHA</option><option value="BETA">BETA</option><option value="CL">CL</option><option value="CY">CY</option>
    </select>
    <select id="stabilityCoeffSelect" hidden aria-hidden="true">
      <option value="CD">CD</option>
      <option value="CY">CY</option>
      <option value="CL" selected>CL</option>
      <option value="CR">CR</option>
      <option value="CM">CM</option>
      <option value="CN">CN</option>
    </select>
  </div>
  <p class="small">Choose axes and ALPHA, BETA, CL or CY. Coefficients plotted against themselves are omitted. L/D uses CL/CD in the chosen axes. Moment plots use the selected reference point.</p>
</div>
<div class="workspace-context" id="coefficientConditions"></div>
<div id="coefficientLegend" class="shared-legend" aria-label="Coefficient curve visibility"></div>
<div id="coefficientPlots" class="analysis-plot-grid">
<div class="coefficient-row">
<div class="card" id="coefficientCLCard"><h2 id="coefficientCLHeading">CL</h2><div id="adfCoeffPlot" class="plot"></div></div>
<div class="card"><h2 id="coefficientCDHeading">CD</h2><div id="coefficientCDPlot" class="plot"></div></div>
<div class="card" id="coefficientCYCard"><h2 id="coefficientCYHeading">CY</h2><div id="coefficientCYPlot" class="plot"></div></div>
</div><div class="workspace-context" id="momentRowReference"></div><div class="coefficient-row">
<div class="card"><h2 id="coefficientCMHeading">CM</h2><div id="coefficientCMPlot" class="plot"></div></div>
<div class="card"><h2 id="coefficientCRHeading">CR</h2><div id="coefficientCRPlot" class="plot"></div></div>
<div class="card"><h2 id="coefficientCNHeading">CN</h2><div id="coefficientCNPlot" class="plot"></div></div>
</div><details id="additionalCoefficientPlots" class="workspace-details"><summary>Additional plots · L/D</summary><div class="coefficient-row">
<div class="card"><h2 id="ldHeading">L/D</h2><div id="ldPlot" class="plot"></div></div>
</div></details>
</div>
</section>

<section id="section-static-margin" class="dashboard-section">
<div class="card">
  <h2>Static margin analysis</h2>
  <p class="small">The active moment reference is controlled in the Coefficients tab. Static-margin plots update automatically when that reference changes.</p>
  <button onclick="showSection('aero', document.querySelector('#dashboardNav button[data-section=&quot;aero&quot;]'));document.getElementById('referencePanel').open=true;document.getElementById('referencePanel').scrollIntoView({{block:'center'}})">Open moment-reference controls</button>
</div>
<div class="analysis-plot-grid">
<div class="card"><h2 id="smClHeading">Static margin vs CLS</h2><div id="smClPlot" class="plot"></div></div>
<div class="card"><h2 id="smSweepHeading">Static margin vs ALPHA/BETA</h2><div id="smSweepPlot" class="plot"></div></div>
</div>
</section>

<section id="section-comparison" class="dashboard-section">
<div class="card">
  <h2>Delta comparisons</h2>
  <p class="small">Difference is Comparison − Reference, interpolated onto the reference curve's selected coordinate without extrapolation. Repeated coordinates are averaged. Static margin retains its stability-axis definition; its horizontal coordinate follows your selection.</p>
  <div class="controls">
    <label for="deltaAxis">Axis system:</label><select id="deltaAxis" onchange="drawComparisonPlot();saveLastState();"><option value="W">Wind</option><option value="S" selected>Stability</option><option value="B">Body</option></select>
    <label for="deltaAbscissa">Plot against:</label><select id="deltaAbscissa" onchange="drawComparisonPlot();saveLastState();"><option>ALPHA</option><option>BETA</option><option>CL</option><option>CY</option></select>
  </div>
  <div class="controls comparison-buttons" id="comparisonMetricButtons">
    <button class="active" onclick="setComparisonMetric('CL', this)">ΔCL</button>
    <button onclick="setComparisonMetric('CD', this)">ΔCD</button>
    <button onclick="setComparisonMetric('CY', this)">ΔCY</button>
    <button onclick="setComparisonMetric('CM', this)">ΔCM</button>
    <button onclick="setComparisonMetric('CR', this)">ΔCR</button>
    <button onclick="setComparisonMetric('CN', this)">ΔCN</button>
    <button onclick="setComparisonMetric('LD', this)">ΔL/D</button>
    <button onclick="setComparisonMetric('SM', this)">ΔStatic margin</button>
  </div>
  <div class="controls">
    <button onclick="addComparisonPair()">Add comparison</button>
    <button onclick="removeLastComparisonPair()">Remove last</button>
    <button onclick="clearComparisonPairs()">Clear</button>
  </div>
  <div id="comparisonPairsContainer"></div>
  <div id="comparisonPlot" class="plot"></div>
</div>
</section>

<section id="section-drag" class="dashboard-section">
<div class="card">
  <h2 id="dragRiseHeading">Drag rise: ΔCDS vs Mach</h2>
  <label for="dragRiseAxis">Axis system:</label><select id="dragRiseAxis" onchange="drawDragRisePlots();saveLastState();"><option value="W">Wind</option><option value="S" selected>Stability</option><option value="B">Body</option></select>
  <p id="dragRiseNote" class="small" role="status"></p>
  <p class="small">Optional. To enable, add <code>"drag_rise_dir": "folder_name"</code> to a CASES entry. The script reads <code>03-RESULTS/DRAG-RISE/folder_name/drag_rise_cl*.dat</code>. One plot is created for each CLS file group.</p>
  <div id="dragRisePlots" class="analysis-plot-grid"></div>
</div>
</section>

<section id="section-convergence" class="dashboard-section">
<div class="analysis-plot-grid">
<div class="card"><h2>All final residual equations vs alpha</h2>
  <div class="controls">
    <label>Residual equation:</label>
    <select id="residualEquationSelect" onchange="drawResidualEquationsPlot(getFilteredRows())">
      <option value="ALL" selected>All</option>
      <option value="continuity">Continuity</option>
      <option value="x-velocity">X-velocity</option>
      <option value="y-velocity">Y-velocity</option>
      <option value="z-velocity">Z-velocity</option>
      <option value="energy">Energy</option>
      <option value="nut">Nut</option>
    </select>
    <button onclick="setResidualEquation('ALL')">All</button>
    <button onclick="setResidualEquation('continuity')">Continuity</button>
    <button onclick="setResidualEquation('x-velocity')">X-velocity</button>
    <button onclick="setResidualEquation('y-velocity')">Y-velocity</button>
    <button onclick="setResidualEquation('z-velocity')">Z-velocity</button>
    <button onclick="setResidualEquation('energy')">Energy</button>
    <button onclick="setResidualEquation('nut')">Nut</button>
  </div>
  <div id="residualEquationsPlot" class="plot"></div>
</div>
<div class="card"><h2>cp-max diagnostics</h2><div id="cpmaxPlot" class="plot"></div></div>
</div>

<div class="card">
  <h2>Detailed selected case history</h2>
  <div class="controls"><label>Case:</label><select id="caseHistoryFilter" onchange="drawSelectedHistory(); refreshAll();"></select></div>
  <div id="selectedCaseSummary" class="status-panel"></div>
  <h3>Final-window statistics</h3>
  <div class="table-scroll"><table id="finalWindowStatsTable"><thead><tr><th>Monitor</th><th>Mean</th><th>Minimum</th><th>Maximum</th><th>Standard deviation</th><th>Drift / 100 iters</th></tr></thead><tbody></tbody></table></div>
  <div class="coefficient-row" id="selectedHistoryPlots">
  <div><div id="selectedResidualHistory" class="plot"></div></div>
  <div><div id="selectedAeroHistory" class="plot"></div></div>
  <div><div id="selectedCpHistory" class="plot"></div></div>
  </div>
</div>
</section>

<section id="section-quality" class="dashboard-section">
<div class="card validation-card">
  <h2>Validated analysis modules</h2>
  <p><strong>Moment-reference transfer: VALIDATED</strong></p>
  <p class="small">Body-axis reference transfer followed by Body → Stability/Wind vector transformation. Version 1.0 validated against an independent source.</p>
</div>
<div class="card">
  <h2>Convergence classification 2.0 by case</h2>
  <div class="table-scroll"><table id="classificationTable"><thead><tr>
    <th>Configuration</th><th>POLAR</th><th>Case</th><th>Alpha</th><th>Beta</th><th>Numerical</th><th>Aerodynamic</th><th>Neighbor</th><th>Overall</th><th>Reasons</th>
  </tr></thead><tbody></tbody></table></div>
</div>
<div class="card">
  <h2>Mesh quality</h2>
  <p class="small">Minimum orthogonal quality below 0.01 is flagged as a warning. Maximum aspect ratio is reported but not automatically failed because high boundary-layer aspect ratios can be intentional.</p>
  <div class="table-scroll"><table id="meshQualityTable"><thead><tr>
    <th>Configuration</th><th>POLAR</th><th>Cells</th><th>Nodes</th><th>Min orthogonal quality</th><th>Max aspect ratio</th><th>Negative-volume warnings</th><th>Mesh status</th><th>Warnings</th>
  </tr></thead><tbody></tbody></table></div>
</div>
<div class="card">
  <h2>Automatic neighbor outlier detection</h2>
  <p class="small">Interior sweep points are compared with a linear interpolation between their adjacent points. A normalized score above 1 exceeds at least one configured CL/CD/CM threshold.</p>
  <div id="outlierScorePlot" class="plot"></div>
  <div class="table-scroll"><table id="outlierTable"><thead><tr>
    <th>Configuration</th><th>POLAR</th><th>Case</th><th>Sweep</th><th>Source</th><th>Alpha</th><th>Beta</th><th>Score</th><th>ΔCL</th><th>ΔCD</th><th>ΔCM</th><th>Status</th><th>Reason</th>
  </tr></thead><tbody></tbody></table></div>
</div>
<div class="card">
  <h2>Data integrity checks</h2>
  <div class="table-scroll"><table id="integrityTable"><thead><tr><th>Severity</th><th>Configuration</th><th>POLAR/file</th><th>Check</th><th>Details</th></tr></thead><tbody></tbody></table></div>
</div>
<div class="card">
  <h2>Data provenance</h2>
  <div id="provenanceSummary" class="small"></div>
  <div class="table-scroll"><table id="provenanceTable"><thead><tr><th>Type</th><th>Configuration</th><th>POLAR/file</th><th>Path</th><th>Modified UTC</th></tr></thead><tbody></tbody></table></div>
</div>
<div class="card">
  <h2>Parser diagnostics</h2>
  <table id="diagTable"><thead><tr><th>Configuration</th><th>POLAR</th><th>Headers</th><th>Rows found</th><th>Total expected</th><th>Consumed</th><th>Remaining</th><th>Log path</th></tr></thead><tbody></tbody></table>
</div>
</section>

<section id="section-tables" class="dashboard-section">
<div class="card">
  <h2>Convergence cases</h2>
  <table id="caseTable"><thead><tr>
    <th>Configuration</th><th>POLAR</th><th>Case</th><th>Name</th><th>Alpha</th><th>Beta</th><th>Numerical</th><th>Aerodynamic</th><th>Neighbor</th><th>Overall</th><th>Quality</th>
    <th>Actual iters</th><th>Continuity</th><th>X-vel</th><th>Y-vel</th><th>Z-vel</th><th>Energy</th><th>Nut</th><th>CLZB</th><th>CDXB</th><th>CMYB</th><th>cp-max</th><th>Reasons</th>
  </tr></thead><tbody></tbody></table>
</div>
</section>

<script>
const summaries = {summaries_json};
const convRows = {conv_json};
const historyData = {history_json};
const adfData = {adf_json};
const dragRiseData = {drag_json};
const moduleVersions = {module_versions_json};
const scriptVersion = {script_version_json};
const provenanceData = {provenance_json};
const integrityData = {integrity_json};
const JAMAL_DEFAULTS = {{
  residualSuccessful: {RESIDUAL_SUCCESSFUL}, residualAcceptable: {RESIDUAL_ACCEPTABLE}, residualBad: {RESIDUAL_BAD},
  clRange: {CL_WINDOW_RANGE_MAX}, cdRange: {CD_WINDOW_RANGE_MAX}, cmRange: {CM_WINDOW_RANGE_MAX}, cpRange: {CPMAX_WINDOW_RANGE_MAX},
  neighborCl: {NEIGHBOR_CL_JUMP_MAX}, neighborCd: {NEIGHBOR_CD_JUMP_MAX}, neighborCm: {NEIGHBOR_CM_JUMP_MAX}
}};

// Mixed typography: Claude-inspired serif for words, original sans-serif for numbers.
const WORD_FONT = 'ui-serif, Georgia, Cambria, "Times New Roman", serif';
const NUMBER_FONT = 'Arial, Helvetica, sans-serif';
const hiddenCoefficientCurves = new Set();

function installPlotActions(element) {{
  if(!element?.matches('.dashboard-section .plot')||element.closest('#section-distributions')||element.dataset.actionsReady) return;
  element.dataset.actionsReady='true';
  const actions=document.createElement('div');actions.className='plot-actions';
  ['Expand','PNG','SVG'].forEach(label=>{{
    const button=document.createElement('button');button.textContent=label;
    button.setAttribute('aria-label',`${{label}} ${{element.id}}`);
    button.onclick=()=>{{
      if(!element.data?.length) return;
      const layout=JSON.parse(JSON.stringify(element.layout));layout.showlegend=true;
      if(label!=='Expand') {{
        const previous=element.layout.showlegend;
        Plotly.relayout(element,{{showlegend:true}}).then(()=>Plotly.downloadImage(element,{{format:label.toLowerCase(),filename:`JAMAL_${{element.id}}`,width:1400,height:900}})).finally(()=>Plotly.relayout(element,{{showlegend:previous}}));return;
      }}
      const dialog=document.getElementById('plotDialog');
      document.getElementById('expandedPlotTitle').textContent=element.layout.yaxis?.title?.text||element.id;
      dialog.showModal();delete layout.width;delete layout.height;layout.autosize=true;
      Plotly.newPlot('expandedPlot',JSON.parse(JSON.stringify(element.data)),layout,{{responsive:true}});
    }};
    actions.appendChild(button);
  }});
  element.before(actions);
}}

function setupPlotFirstWorkspace() {{
  const section=document.getElementById('section-aero'), controls=section.querySelector('.coefficient-controls');
  section.prepend(controls);
  section.appendChild(document.getElementById('deflectionPanel'));
  const dialog=document.createElement('dialog');dialog.id='plotDialog';
  dialog.innerHTML='<header><h2 id="expandedPlotTitle">Expanded plot</h2><button type="button">Close</button></header><div id="expandedPlot" class="plot"></div>';
  document.body.appendChild(dialog);
  dialog.querySelector('button').onclick=()=>dialog.close();
  dialog.addEventListener('close',()=>Plotly.purge('expandedPlot'));
  document.getElementById('additionalCoefficientPlots').addEventListener('toggle',()=>{{if(document.getElementById('additionalCoefficientPlots').open) Plotly.Plots.resize('ldPlot');}});
}}

function mixedTypographyLayout(layout) {{
  const out = Object.assign({{}}, layout || {{}});
  out.font = Object.assign({{}}, out.font || {{}}, {{family: NUMBER_FONT}});

  if (typeof out.title === 'string') {{
    out.title = {{text: out.title, font: {{family: WORD_FONT}}}};
  }} else if (out.title && typeof out.title === 'object') {{
    out.title = Object.assign({{}}, out.title, {{font: Object.assign({{}}, out.title.font || {{}}, {{family: WORD_FONT}})}});
  }}

  ['xaxis', 'yaxis', 'yaxis2', 'xaxis2'].forEach(axisName => {{
    if (!out[axisName]) return;
    const axis = Object.assign({{}}, out[axisName]);
    axis.tickfont = Object.assign({{}}, axis.tickfont || {{}}, {{family: NUMBER_FONT}});
    if (typeof axis.title === 'string') {{
      axis.title = {{text: axis.title, font: {{family: WORD_FONT}}}};
    }} else if (axis.title && typeof axis.title === 'object') {{
      axis.title = Object.assign({{}}, axis.title, {{font: Object.assign({{}}, axis.title.font || {{}}, {{family: WORD_FONT}})}});
    }}
    out[axisName] = axis;
  }});

  // Legends and hover labels often contain POLAR numbers, so keep the numeric font.
  out.legend = Object.assign({{}}, out.legend || {{}}, {{font: Object.assign({{}}, (out.legend || {{}}).font || {{}}, {{family: NUMBER_FONT}})}});
  out.hoverlabel = Object.assign({{}}, out.hoverlabel || {{}}, {{font: Object.assign({{}}, (out.hoverlabel || {{}}).font || {{}}, {{family: NUMBER_FONT}})}});
  return out;
}}

const originalPlotlyNewPlot = Plotly.newPlot.bind(Plotly);
const originalPlotlyReact = Plotly.react.bind(Plotly);
Plotly.newPlot = function(gd, data, layout, config) {{
  const element = typeof gd === "string" ? document.getElementById(gd) : gd;
  const styled = mixedTypographyLayout(layout);
  installPlotActions(element);
  styled.uirevision = JSON.stringify([styled.xaxis?.title, styled.yaxis?.title]);
  return element?._fullLayout
    ? originalPlotlyReact(gd, data, styled, config)
    : originalPlotlyNewPlot(gd, data, styled, config);
}};
const staticMarginYlim = {sm_ylim_json};
const filterStaticMargin = {filter_sm_json};
const residualCols = ["continuity", "x-velocity", "y-velocity", "z-velocity", "energy", "nut"];
const stabilityCoeffs = ["CDS", "CYS", "CLS", "CRS25", "CMS25", "CNS25"];

function unique(arr) {{ return [...new Set(arr)].sort(); }}
function statusColor(status) {{
  if (status === "CONVERGED") return "#1a8f3a";
  if (status === "ACCEPTABLE") return "#2878b5";
  if (status === "SUSPICIOUS") return "#d98200";
  return "#d62728";
}}
function fmt(x, digits=5) {{ if (x === null || x === undefined || isNaN(x)) return ""; if (Math.abs(x) < 1e-3 || Math.abs(x) > 1e4) return Number(x).toExponential(3); return Number(x).toFixed(digits); }}
function sciAxis(title) {{ return {{title: title, type: "log", exponentformat: "e", showexponent: "all", tickformat: ".0e"}}; }}
function assessmentClass(value) {{ return "assessment-" + String(value || "UNKNOWN").replace(/\\s+/g, "_"); }}
function showSection(sectionName, button) {{
  document.querySelectorAll(".dashboard-section").forEach(el => el.classList.remove("active"));
  const target = document.getElementById("section-" + sectionName);
  if (target) target.classList.add("active");
  document.querySelectorAll("#dashboardNav button").forEach(el => el.classList.remove("active"));
  const navButton = button || document.querySelector(`#dashboardNav button[data-section="${{sectionName}}"]`);
  if (navButton) navButton.classList.add("active");
  setTimeout(() => {{
    if (!target) return;
    target.querySelectorAll(".plot").forEach(div => {{
      try {{ Plotly.Plots.resize(div); }} catch (err) {{}}
    }});
  }}, 30);
}}

function summaryFor(label, polar) {{
  return summaries.find(s => s.case_label === label && s.polar === polar) || null;
}}

function firstSummary() {{
  return summaries.length ? summaries[0] : null;
}}

function initializeMomentReferenceInputs() {{
  const s = firstSummary();
  if (!s || !s.meta) return;
  document.getElementById("xrefPercentInput").value = 25.0;
  document.getElementById("yOffsetInput").value = 0.0;
  document.getElementById("zOffsetInput").value = 0.0;
  document.getElementById("xAbsInput").value = Number(s.meta.xref ?? 0).toFixed(6);
  document.getElementById("yAbsInput").value = Number(s.meta.yref ?? 0).toFixed(6);
  document.getElementById("zAbsInput").value = Number(s.meta.zref ?? 0).toFixed(6);
}}

function resetMomentReference() {{
  const s = firstSummary();
  document.getElementById("refModeSelect").value = "percent";
  document.getElementById("xrefPercentInput").value = 25.0;
  document.getElementById("yOffsetInput").value = 0.0;
  document.getElementById("zOffsetInput").value = 0.0;
  if (s && s.meta) {{
    document.getElementById("xAbsInput").value = Number(s.meta.xref ?? 0).toFixed(6);
    document.getElementById("yAbsInput").value = Number(s.meta.yref ?? 0).toFixed(6);
    document.getElementById("zAbsInput").value = Number(s.meta.zref ?? 0).toFixed(6);
  }}
  refreshAll();
}}

function currentMomentInputs() {{
  return {{
    mode: document.getElementById("refModeSelect").value,
    xPercent: Number(document.getElementById("xrefPercentInput").value),
    yOffset: Number(document.getElementById("yOffsetInput").value),
    zOffset: Number(document.getElementById("zOffsetInput").value),
    xAbs: Number(document.getElementById("xAbsInput").value),
    yAbs: Number(document.getElementById("yAbsInput").value),
    zAbs: Number(document.getElementById("zAbsInput").value)
  }};
}}

function bodyToStabilityVector(dx, dy, dz, alphaDeg) {{
  // Body -> Stability vector transform.
  // For beta = 0 this is a rotation about the common Y axis.
  const a = (alphaDeg || 0) * Math.PI / 180.0;
  const ca = Math.cos(a), sa = Math.sin(a);
  return {{
    x: dx * ca + dz * sa,
    y: dy,
    z: -dx * sa + dz * ca
  }};
}}

function bodyToWindVector(dx, dy, dz, alphaDeg, betaDeg) {{
  // Body -> Wind vector transform. This matches the earlier combined alpha/beta
  // transform used by the dashboard for a vector quantity.
  const a = (alphaDeg || 0) * Math.PI / 180.0;
  const b = (betaDeg || 0) * Math.PI / 180.0;
  const ca = Math.cos(a), sa = Math.sin(a);
  const cb = Math.cos(b), sb = Math.sin(b);
  return {{
    x: dx * ca * cb + dy * sb + dz * sa * cb,
    y: -dx * ca * sb + dy * cb - dz * sa * sb,
    z: -dx * sa + dz * ca
  }};
}}

function newReferencePoint(meta, inp) {{
  const xref = Number(meta.xref ?? 0);
  const yref = Number(meta.yref ?? 0);
  const zref = Number(meta.zref ?? 0);
  const cref = Number(meta.cref ?? 1);

  if (inp.mode === "absolute") {{
    return {{xnew: inp.xAbs, ynew: inp.yAbs, znew: inp.zAbs}};
  }}

  // Per-POLAR %MAC mode. XREF is assumed to be 25% MAC for each POLAR.
  // 0% MAC = XREF - 0.25*CREF, 25% MAC = XREF, 50% MAC = XREF + 0.25*CREF.
  const xnew = xref + ((inp.xPercent - 25.0) / 100.0) * cref;
  return {{xnew: xnew, ynew: yref + inp.yOffset, znew: zref + inp.zOffset}};
}}

function shiftedRow(row, meta) {{
  const out = Object.assign({{}}, row);
  if (!meta) return out;

  const inp = currentMomentInputs();
  const xref = Number(meta.xref ?? 0);
  const yref = Number(meta.yref ?? 0);
  const zref = Number(meta.zref ?? 0);
  const cref = Number(meta.cref ?? 1);
  const bref = Number(meta.bref ?? 1);
  const p = newReferencePoint(meta, inp);

  // Reference displacement in the infout geometry frame:
  //   X_info positive backward, Z_info positive upward.
  const dxInfo = p.xnew - xref;
  const dyInfo = p.ynew - yref;
  const dzInfo = p.znew - zref;

  // Convert to the body/aero calculation frame:
  //   X_body positive forward, Y same, Z_body positive downward.
  const dxBody = -dxInfo;
  const dyBody =  dyInfo;
  const dzBody = -dzInfo;

  out.XREF_USER_PERCENT_MAC = inp.xPercent;
  out.XREF_USER_ABSOLUTE = p.xnew;
  out.YREF_USER_ABSOLUTE = p.ynew;
  out.ZREF_USER_ABSOLUTE = p.znew;
  out.DX_USER_REF_INFO = dxInfo;
  out.DY_USER_REF_INFO = dyInfo;
  out.DZ_USER_REF_INFO = dzInfo;
  out.DX_USER_REF = dxBody;
  out.DY_USER_REF = dyBody;
  out.DZ_USER_REF = dzBody;

  // If Body-axis force/moment columns are not available, keep original row values.
  const hasBody = row.CDB !== undefined && row.CYB !== undefined && row.CLB !== undefined &&
                  row.CRB25 !== undefined && row.CMB25 !== undefined && row.CNB25 !== undefined;
  if (!hasBody) return out;

  // Store original moments for debugging/hover if needed.
  out.CRB25_ORIGINAL = Number(row.CRB25);
  out.CMB25_ORIGINAL = Number(row.CMB25);
  out.CNB25_ORIGINAL = Number(row.CNB25);
  if (row.CRS25 !== undefined) out.CRS25_ORIGINAL = Number(row.CRS25);
  if (row.CMS25 !== undefined) out.CMS25_ORIGINAL = Number(row.CMS25);
  if (row.CNS25 !== undefined) out.CNS25_ORIGINAL = Number(row.CNS25);
  if (row.CRW25 !== undefined) out.CRW25_ORIGINAL = Number(row.CRW25);
  if (row.CMW25 !== undefined) out.CMW25_ORIGINAL = Number(row.CMW25);
  if (row.CNW25 !== undefined) out.CNW25_ORIGINAL = Number(row.CNW25);

  // Convert intuitive Body-axis force coefficients back to physical force components.
  // Body axes: +X forward, +Z downward. Drag and lift are reported positive.
  const fxB = -Number(row.CDB);
  const fyB =  Number(row.CYB);
  const fzB = -Number(row.CLB);

  // Body-axis moment transfer: M_new = M_old - r x F.
  const crb = Number(row.CRB25) - (dyBody * fzB - dzBody * fyB) / bref;
  const cmb = Number(row.CMB25) - (dzBody * fxB - dxBody * fzB) / cref;
  const cnb = Number(row.CNB25) - (dxBody * fyB - dyBody * fxB) / bref;

  out.CRB25 = crb;
  out.CMB25 = cmb;
  out.CNB25 = cnb;

  // Forces are unchanged by a reference-point shift. Keep ADF force coefficients as-is.
  // Rotate the corrected Body moment vector to Stability and Wind axes.
  const ms = bodyToStabilityVector(crb, cmb, cnb, row.ALPHA || 0);
  const mw = bodyToWindVector(crb, cmb, cnb, row.ALPHA || 0, row.BETA || 0);

  out.CRS25 = ms.x;
  out.CMS25 = ms.y;
  out.CNS25 = ms.z;
  out.CRW25 = mw.x;
  out.CMW25 = mw.y;
  out.CNW25 = mw.z;

  out.RX_BODY = dxBody;
  out.RY_BODY = dyBody;
  out.RZ_BODY = dzBody;
  const rs = bodyToStabilityVector(dxBody, dyBody, dzBody, row.ALPHA || 0);
  const rw = bodyToWindVector(dxBody, dyBody, dzBody, row.ALPHA || 0, row.BETA || 0);
  out.RX_STABILITY = rs.x;
  out.RY_STABILITY = rs.y;
  out.RZ_STABILITY = rs.z;
  out.RX_WIND = rw.x;
  out.RY_WIND = rw.y;
  out.RZ_WIND = rw.z;
  return out;
}}

function currentMomentReferenceMode() {{
  const el = document.getElementById("momentReferenceSelect");
  return el ? el.value : "original";
}}

function currentAdfCurves() {{
  const refMode = currentMomentReferenceMode();
  return filteredAdfCurves().map(c => {{
    if (refMode === "original") {{
      return Object.assign({{}}, c, {{rows: c.rows.map(r => Object.assign({{}}, r))}});
    }}
    const s = summaryFor(c.case_label, c.polar);
    const meta = s ? s.meta : null;
    return Object.assign({{}}, c, {{rows: c.rows.map(r => shiftedRow(r, meta))}});
  }});
}}

function derivativeNonuniform(xs, ys) {{
  const n = xs.length;
  if (n < 3) return Array(n).fill(null);
  const d = Array(n).fill(null);
  function derivAt(i0, i1, i2, x) {{
    const x0=xs[i0], x1=xs[i1], x2=xs[i2];
    const y0=ys[i0], y1=ys[i1], y2=ys[i2];
    const w0 = (2*x - x1 - x2) / ((x0-x1)*(x0-x2));
    const w1 = (2*x - x0 - x2) / ((x1-x0)*(x1-x2));
    const w2 = (2*x - x0 - x1) / ((x2-x0)*(x2-x1));
    return w0*y0 + w1*y1 + w2*y2;
  }}
  d[0] = derivAt(0,1,2,xs[0]);
  for (let i=1; i<n-1; i++) d[i] = derivAt(i-1,i,i+1,xs[i]);
  d[n-1] = derivAt(n-3,n-2,n-1,xs[n-1]);
  return d;
}}

function currentStaticMarginRows() {{
  const out = [];
  currentAdfCurves().forEach(c => {{
    let rows = c.rows.filter(r => r.CLS !== undefined && r.CMS25 !== undefined).slice().sort((a,b) => a.CLS - b.CLS);
    // Average duplicate CLS values to avoid zero spacing.
    const reduced = [];
    rows.forEach(r => {{
      const last = reduced.length ? reduced[reduced.length-1] : null;
      if (last && Math.abs(last.CLS - r.CLS) < 1e-12) {{
        last.CMS25 = 0.5*(last.CMS25 + r.CMS25);
        last.ALPHA = 0.5*((last.ALPHA ?? 0) + (r.ALPHA ?? 0));
        last.BETA = 0.5*((last.BETA ?? 0) + (r.BETA ?? 0));
      }} else {{
        reduced.push(Object.assign({{}}, r));
      }}
    }});
    if (reduced.length < 3) return;
    const xs = reduced.map(r => Number(r.CLS));
    const ys = reduced.map(r => Number(r.CMS25));
    const d = derivativeNonuniform(xs, ys);
    reduced.forEach((r,i) => {{
      const sm = -d[i];
      if (!Number.isFinite(sm)) return;
      out.push(Object.assign({{}}, r, {{
        case_label: c.case_label, polar: c.polar, curve_label: c.curve_label, sweep_var: c.sweep_var,
        STATIC_MARGIN: sm, STATIC_MARGIN_PERCENT: 100.0*sm, dCMS25_dCLS: d[i]
      }}));
    }});
  }});
  let filtered = out;
  if (filterStaticMargin && staticMarginYlim !== null) {{
    filtered = filtered.filter(r => r.STATIC_MARGIN_PERCENT >= staticMarginYlim[0] && r.STATIC_MARGIN_PERCENT <= staticMarginYlim[1]);
  }}
  return filtered;
}}

function updateMomentRefInfo() {{
  const inp = currentMomentInputs();
  const rowsInfo = summaries.map(s => {{
    const meta = s.meta || {{}};
    const xref = Number(meta.xref ?? 0);
    const yref = Number(meta.yref ?? 0);
    const zref = Number(meta.zref ?? 0);
    const cref = Number(meta.cref ?? 1);
    const xleMac = xref - 0.25 * cref;
    const p = newReferencePoint(meta, inp);
    const dxInfo = p.xnew - xref;
    const dyInfo = p.ynew - yref;
    const dzInfo = p.znew - zref;
    const dxCalc = -dxInfo;
    const dyCalc =  dyInfo;
    const dzCalc = -dzInfo;
    return `<tr><td>${{s.case_label}}</td><td>${{s.polar}}</td><td>${{xref.toFixed(4)}}</td><td>${{yref.toFixed(4)}}</td><td>${{zref.toFixed(4)}}</td><td>${{cref.toFixed(4)}}</td><td>${{xleMac.toFixed(4)}}</td><td>${{p.xnew.toFixed(4)}}</td><td>${{p.ynew.toFixed(4)}}</td><td>${{p.znew.toFixed(4)}}</td><td>${{dxInfo.toFixed(4)}}</td><td>${{dyInfo.toFixed(4)}}</td><td>${{dzInfo.toFixed(4)}}</td><td>${{dxCalc.toFixed(4)}}</td><td>${{dyCalc.toFixed(4)}}</td><td>${{dzCalc.toFixed(4)}}</td></tr>`;
  }}).join("");
  const modeText = inp.mode === "absolute"
    ? `Absolute XYZ mode: all POLARs use X=${{inp.xAbs.toFixed(4)}} m, Y=${{inp.yAbs.toFixed(4)}} m, Z=${{inp.zAbs.toFixed(4)}} m in the infout coordinate frame.`
    : `%MAC per POLAR mode: each POLAR uses XNEW = XREF + ((${{inp.xPercent.toFixed(2)}} - 25)/100)*CREF, with Y/Z offsets ${{inp.yOffset.toFixed(4)}} / ${{inp.zOffset.toFixed(4)}} m in the infout coordinate frame.`;
  const plotRefMode = currentMomentReferenceMode && currentMomentReferenceMode() === "user" ? "User reference point" : "Original ADF reference";
  document.getElementById("momentRefInfo").innerHTML =
    `${{modeText}} Current stability-coefficient plot mode: <b>${{plotRefMode}}</b>. Infout frame: X positive backward, Z positive upward. Moment calculation frame: X positive forward, Z positive downward. When User reference is selected, moment correction is applied once in Body axes, then the corrected moment vector is rotated to Stability and Wind axes.` +
    `<table style="margin-top:8px;"><thead><tr><th>Configuration</th><th>POLAR</th><th>Original XREF</th><th>Original YREF</th><th>Original ZREF</th><th>CREF</th><th>MAC LE</th><th>New XREF</th><th>New YREF</th><th>New ZREF</th><th>dx info</th><th>dy info</th><th>dz info</th><th>dx calc</th><th>dy calc</th><th>dz calc</th></tr></thead><tbody>${{rowsInfo}}</tbody></table>`;
}}

function setupFilters() {{
  const labels = ["ALL"].concat(unique([...convRows,...adfData.curves].map(r => r.case_label)));
  const polars = ["ALL"].concat(unique([...convRows,...adfData.curves].map(r => r.polar)));
  document.getElementById("labelFilter").innerHTML = labels.map(v => `<option value="${{v}}">${{v}}</option>`).join("");
  document.getElementById("polarFilter").innerHTML = polars.map(v => `<option value="${{v}}">${{v}}</option>`).join("");
  setupCaseHistoryFilter();
}}

function setupCaseHistoryFilter() {{
  const select = document.getElementById("caseHistoryFilter");
  const sorted = convRows.slice().sort((a,b) => (a.case_label + a.polar + a.case).localeCompare(b.case_label + b.polar + b.case));
  select.innerHTML = sorted.map(r => {{
    const label = `${{r.case_label}} | ${{r.polar}} | CASE ${{r.case}} | α=${{r.alpha}} | ${{r.status}}`;
    return `<option value="${{r.case_key}}">${{label}}</option>`;
  }}).join("");
}}

function getFilteredRows() {{
  const label = document.getElementById("labelFilter").value;
  const polar = document.getElementById("polarFilter").value;
  const status = document.getElementById("statusFilter").value;
  return convRows.filter(r => {{
    if (label !== "ALL" && r.case_label !== label) return false;
    if (polar !== "ALL" && r.polar !== polar) return false;
    if (status !== "ALL" && r.status !== status) return false;
    return true;
  }});
}}

function filteredSummaries() {{
  const label = document.getElementById("labelFilter").value;
  const polar = document.getElementById("polarFilter").value;
  return summaries.filter(s => {{
    if (label !== "ALL" && s.case_label !== label) return false;
    if (polar !== "ALL" && s.polar !== polar) return false;
    return true;
  }});
}}


function filteredAdfCurves() {{
  const label = document.getElementById("labelFilter").value;
  const polar = document.getElementById("polarFilter").value;
  return adfData.curves.filter(c => {{
    if (label !== "ALL" && c.case_label !== label) return false;
    if (polar !== "ALL" && c.polar !== polar) return false;
    return true;
  }});
}}

function filteredSmRows() {{
  return currentStaticMarginRows().filter(r => {{
    const label = document.getElementById("labelFilter").value;
    const polar = document.getElementById("polarFilter").value;
    if (label !== "ALL" && r.case_label !== label) return false;
    if (polar !== "ALL" && r.polar !== polar) return false;
    return true;
  }});
}}
function groupedByLabelPolar(data) {{
  const groups = {{}};
  data.forEach(r => {{ const key = r.case_label + " | " + r.polar; if (!groups[key]) groups[key] = []; groups[key].push(r); }});
  Object.keys(groups).forEach(k => groups[k].sort((a,b) => a.alpha - b.alpha));
  return groups;
}}

function drawKpis(data) {{
  const total = data.length;
  const conv = data.filter(r => r.status === "CONVERGED").length;
  const acc = data.filter(r => r.status === "ACCEPTABLE").length;
  const susp = data.filter(r => r.status === "SUSPICIOUS").length;
  const divg = data.filter(r => r.status === "DIVERGED").length;
  const meshWarn = filteredSummaries().filter(s => ["WARNING", "CRITICAL"].includes((s.mesh || {{}}).status)).length;
  document.getElementById('compactStatus').textContent=convRows.length?`${{total}} cases · ${{conv}} converged · ${{acc}} acceptable · ${{susp}} suspicious · ${{divg}} diverged`:`${{filteredAdfCurves().length}} ADF polars · Convergence history unavailable`;
  document.getElementById("kpis").innerHTML = `
    <div class="kpi"><div>Total cases</div><strong>${{total}}</strong></div>
    <div class="kpi"><div>Converged</div><strong>${{conv}}</strong></div>
    <div class="kpi"><div>Acceptable</div><strong>${{acc}}</strong></div>
    <div class="kpi"><div>Suspicious</div><strong>${{susp}}</strong></div>
    <div class="kpi"><div>Diverged</div><strong>${{divg}}</strong></div>
    <div class="kpi"><div>Mesh warnings</div><strong>${{meshWarn}}</strong></div>`;
}}

function drawAssessmentKpis(data) {{
  const nSuccessful = data.filter(r => r.numerical_status === "SUCCESSFUL").length;
  const nAcceptable = data.filter(r => r.numerical_status === "ACCEPTABLE").length;
  const nWarning = data.filter(r => r.numerical_status === "WARNING").length;
  const nDiverged = data.filter(r => r.numerical_status === "DIVERGED").length;
  const aeroWarnings = data.filter(r => r.aerodynamic_status === "WARNING").length;
  const outliers = data.filter(r => r.neighbor_status === "OUTLIER").length;
  document.getElementById("assessmentKpis").innerHTML = `
    <span class="metric-note">Numerical successful: <strong>${{nSuccessful}}</strong></span>
    <span class="metric-note">Numerical acceptable: <strong>${{nAcceptable}}</strong></span>
    <span class="metric-note">Numerical warnings: <strong>${{nWarning}}</strong></span>
    <span class="metric-note">Numerical diverged: <strong>${{nDiverged}}</strong></span>
    <span class="metric-note">Aerodynamic warnings: <strong>${{aeroWarnings}}</strong></span>
    <span class="metric-note">Neighbor outliers: <strong>${{outliers}}</strong></span>`;
}}

function drawStatusMap(data) {{
  document.getElementById("statusMap").innerHTML = data.map(r => `<span class="badge" title="${{r.case_label}} ${{r.polar}} ${{r.case}} ${{r.name}} | Numerical: ${{r.numerical_status}} | Aero: ${{r.aerodynamic_status}} | Neighbor: ${{r.neighbor_status}} | Overall: ${{r.status}}" style="background:${{statusColor(r.status)}};">${{r.case_label}} | ${{r.polar}} | α=${{r.alpha}} | ${{r.status}}</span>`).join("");
}}

function drawModuleVersions() {{
  const rows = Object.entries(moduleVersions).map(([name, version]) => `<tr><td>${{name}}</td><td>${{version}}</td></tr>`).join("");
  document.querySelector("#moduleVersionTable tbody").innerHTML = `<tr><td><strong>Dashboard script</strong></td><td><strong>${{scriptVersion}}</strong></td></tr>` + rows;
}}


function coefficientKey(axis, baseCoeff) {{
  if (["CD", "CY", "CL"].includes(baseCoeff)) return baseCoeff + axis;
  return baseCoeff + axis + "25";
}}

function coefficientDisplayName(axis, baseCoeff) {{
  return coefficientKey(axis, baseCoeff);
}}

function setStabilityCoeff(baseCoeff) {{
  document.getElementById("stabilityCoeffSelect").value = baseCoeff;
  drawAdfCoeffPlot();
}}

function drawAdfCoeffPlot() {{
  const curves = currentAdfCurves();
  const selectedAxis = document.getElementById("stabilityAxisSelect").value || "S";
  const selectedBaseCoeff = document.getElementById("stabilityCoeffSelect").value || "CL";
  const selectedCoeff = coefficientKey(selectedAxis, selectedBaseCoeff);
  const coeffLabel = coefficientDisplayName(selectedAxis, selectedBaseCoeff);
  const traces = [];
  let sweepName = "ALPHA";

  curves.forEach(c => {{
    sweepName = c.sweep_var || sweepName;
    if (!c.rows.some(r => r[selectedCoeff] !== undefined)) return;

    const x = c.rows.map(r => r[c.sweep_var]);
    const y = c.rows.map(r => r[selectedCoeff]);
    const hover = c.rows.map(r =>
      `Configuration: ${{c.case_label}}<br>` +
      `POLAR: ${{c.polar}}<br>` +
      `Axis: ${{selectedAxis}}<br>` +
      `Moment reference: ${{currentMomentReferenceMode() === "user" ? "User" : "Original ADF"}}<br>` +
      `${{c.sweep_var}}: ${{fmt(r[c.sweep_var],3)}}<br>` +
      `${{coeffLabel}}: ${{fmt(r[selectedCoeff],6)}}<br>` +
      `X ref: ${{fmt(r.XREF_USER_PERCENT_MAC ?? 25,2)}}% MAC<br>` +
      `dx/dy/dz body: ${{fmt(r.DX_USER_REF ?? 0,4)}} / ${{fmt(r.DY_USER_REF ?? 0,4)}} / ${{fmt(r.DZ_USER_REF ?? 0,4)}} m<extra></extra>`
    );

    traces.push({{
      x:x,
      y:y,
      mode:"lines+markers",
      name: c.case_label + " " + c.polar.replace("POLAR-", "P"),
      text:hover,
      hovertemplate:"%{{text}}"
    }});
  }});

  const refTitle = currentMomentReferenceMode() === "user" ? "User reference" : "Original ADF reference";
  const title = `${{coeffLabel}} vs ${{sweepName}} (${{refTitle}})`;
  Plotly.newPlot("adfCoeffPlot", traces, {{
    title:title,
    xaxis:{{title:`${{sweepName}} [deg]`}},
    yaxis:{{title:coeffLabel}},
    margin:{{l:60,r:30,t:50,b:50}}
  }}, {{responsive:true}});
}}
function drawAdfXY(divId, xKey, yKey, title, yFunc=null, yLabel=null) {{
  const traces = [];
  currentAdfCurves().forEach(c => {{
    const rows = c.rows.filter(r => r[xKey] !== undefined && r[yKey] !== undefined);
    const yvals = rows.map(r => yFunc ? yFunc(r) : r[yKey]);
    traces.push({{
      x: rows.map(r => r[xKey]),
      y: yvals,
      mode:"lines+markers",
      name: c.case_label + " " + c.polar.replace("POLAR-", "P"),
      text: rows.map((r,i) => `Configuration: ${{c.case_label}}<br>POLAR: ${{c.polar}}<br>${{xKey}}: ${{fmt(r[xKey],6)}}<br>${{yLabel || yKey}}: ${{fmt(yvals[i],6)}}<br>X ref: ${{fmt(r.XREF_USER_PERCENT_MAC ?? 25,2)}}% MAC<br>dx/dy/dz body: ${{fmt(r.DX_USER_REF ?? 0,4)}} / ${{fmt(r.DY_USER_REF ?? 0,4)}} / ${{fmt(r.DZ_USER_REF ?? 0,4)}} m<extra></extra>`),
      hovertemplate:"%{{text}}"
    }});
  }});
  Plotly.newPlot(divId, traces, {{title:title, xaxis:{{title:xKey}}, yaxis:{{title:yLabel || yKey}}, margin:{{l:60,r:30,t:50,b:50}}}}, {{responsive:true}});
}}
function drawSmPlots() {{
  const sm = filteredSmRows();
  const groups = {{}};
  sm.forEach(r => {{ const key = r.curve_label; if (!groups[key]) groups[key] = []; groups[key].push(r); }});
  const tracesCl = [];
  const tracesSweep = [];
  let sweepName = "ALPHA";
  Object.keys(groups).forEach(key => {{
    const g1 = groups[key].slice().sort((a,b) => a.CLS - b.CLS);
    const shortName = key.replace("POLAR-", "P");
    tracesCl.push({{
      x:g1.map(r => r.CLS), y:g1.map(r => r.STATIC_MARGIN_PERCENT), mode:"lines+markers", name:shortName,
      text:g1.map(r => `Configuration: ${{r.case_label}}<br>POLAR: ${{r.polar}}<br>CLS: ${{fmt(r.CLS,6)}}<br>SM: ${{fmt(r.STATIC_MARGIN_PERCENT,4)}}% CREF<br>X ref: ${{fmt(r.XREF_USER_PERCENT_MAC ?? 25,2)}}% MAC<extra></extra>`),
      hovertemplate:"%{{text}}"
    }});
    const g2 = groups[key].slice().sort((a,b) => (a.ALPHA ?? a.BETA) - (b.ALPHA ?? b.BETA));
    if (g2.length && g2[0].sweep_var) sweepName = g2[0].sweep_var;
    tracesSweep.push({{
      x:g2.map(r => r.ALPHA !== undefined ? r.ALPHA : r.BETA), y:g2.map(r => r.STATIC_MARGIN_PERCENT), mode:"lines+markers", name:shortName,
      text:g2.map(r => `Configuration: ${{r.case_label}}<br>POLAR: ${{r.polar}}<br>${{r.sweep_var}}: ${{fmt(r.ALPHA !== undefined ? r.ALPHA : r.BETA,3)}}<br>SM: ${{fmt(r.STATIC_MARGIN_PERCENT,4)}}% CREF<br>X ref: ${{fmt(r.XREF_USER_PERCENT_MAC ?? 25,2)}}% MAC<extra></extra>`),
      hovertemplate:"%{{text}}"
    }});
  }});
  const smLayout = {{title:"Static margin [% CREF]", yaxis:{{title:"Static margin [% CREF]"}}, margin:{{l:70,r:30,t:50,b:50}}, shapes:[{{type:"line",xref:"paper",x0:0,x1:1,y0:0,y1:0,line:{{dash:"dash"}}}}]}};
  if (staticMarginYlim !== null && !filterStaticMargin) smLayout.yaxis.range = staticMarginYlim;
  Plotly.newPlot("smClPlot", tracesCl, Object.assign({{}}, smLayout, {{xaxis:{{title:"CLS"}}}}), {{responsive:true}});
  Plotly.newPlot("smSweepPlot", tracesSweep, Object.assign({{}}, smLayout, {{xaxis:{{title:`${{sweepName}} [deg]`}}}}), {{responsive:true}});
}}
function filteredDragRiseCurves() {{
  const label = document.getElementById("labelFilter").value;
  const curves = (dragRiseData && dragRiseData.curves) ? dragRiseData.curves : [];
  return curves.filter(c => {{
    if (label !== "ALL" && c.case_label !== label) return false;
    return true;
  }});
}}

function drawDragRisePlots() {{
  const container = document.getElementById("dragRisePlots");
  if (!container) return;
  const axis=document.getElementById("dragRiseAxis").value, cdKey="CD"+axis;
  document.getElementById("dragRiseHeading").textContent=`Drag rise: Δ${{cdKey}} vs Mach`;
  const missing=[];
  document.getElementById("dragRiseNote").textContent="Each increment uses the selected-axis drag at the lowest Mach point. File groups retain their original CLS target.";

  const curves = filteredDragRiseCurves();

  if (!curves.length) {{
    container.innerHTML = "<p class='small'>No drag-rise directory configured or no drag-rise files found.</p>";
    return;
  }}

  const groups = {{}};
  curves.forEach(c => {{
    const key = c.cls_label || c.file_name;
    if (!groups[key]) groups[key] = [];
    groups[key].push(c);
  }});

  const groupKeys = Object.keys(groups).sort();
  container.innerHTML = groupKeys.map((key, idx) =>
    `<div style="margin-top: 10px;"><h3>${{dragRiseLiftLabel(groups[key],axis)}}</h3><div id="dragRisePlot_${{idx}}" class="plot"></div></div>`
  ).join("");

  groupKeys.forEach((key, idx) => {{
    const traces = groups[key].map(c => {{
      const rows = (c.rows || []).slice().sort((a,b) => Number(a.MACH) - Number(b.MACH));
      const baseline=rows.length?rows[0][cdKey]:null;
      if(!Number.isFinite(baseline)||rows.some(r=>!Number.isFinite(r[cdKey]))) missing.push(c.file_name);
      const st = traceStyle(c.case_label, key);
      return {{
        x: rows.map(r => r.MACH),
        y: dragRiseValues(rows,axis),connectgaps:false,
        mode: st.mode,
        name: c.case_label,
        line: {{color:st.color,dash:st.dash}},
        marker: {{color:st.color,symbol:st.symbol,size:st.markerSize}},
        text: rows.map(r =>
          `Configuration: ${{c.case_label}}<br>` +
          `Drag-rise dir: ${{c.drag_rise_dir}}<br>` +
          `File: ${{c.file_name}}<br>` +
          `MACH: ${{fmt(r.MACH,4)}}<br>` +
          `${{cdKey}}: ${{fmt(r[cdKey],6)}}<br>` +
          `Δ${{cdKey}}: ${{fmt(Number.isFinite(baseline)&&Number.isFinite(r[cdKey])?r[cdKey]-baseline:null,6)}}<br>` +
          `Actual CL${{axis}}: ${{fmt(r["CL"+axis],6)}}<br>` +
          `POLAR: ${{r.POLAR ?? ""}}<extra></extra>`
        ),
        hovertemplate: "%{{text}}"
      }};
    }});

    Plotly.newPlot(`dragRisePlot_${{idx}}`, traces, {{
      title: `Drag rise ${{dragRiseLiftLabel(groups[key],axis)}}`,
      xaxis: {{title: "Mach"}},
      yaxis: {{title: `Δ${{cdKey}}`}},
      margin: {{l:70,r:30,t:50,b:50}}
    }}, {{responsive:true}});
  }});
  if(missing.length) document.getElementById("dragRiseNote").textContent+=` Missing ${{cdKey}} values in ${{[...new Set(missing)].join(", ")}}; unavailable points are omitted.`;
}}

function dragRiseLiftLabel(curves,axis) {{
  const key="CL"+axis;
  const values=curves.flatMap(c=>(c.rows||[]).map(r=>r[key])).filter(Number.isFinite);
  if(!values.length) return `${{key}} unavailable`;
  const low=Math.min(...values), high=Math.max(...values);
  const value=fmt(low,3)===fmt(high,3)?fmt(low,3):`${{fmt(low,3)}} to ${{fmt(high,3)}}`;
  return `${{key}} = ${{value}}`;
}}

function dragRiseValues(rows,axis) {{
  const key="CD"+axis, baseline=rows.length?rows[0][key]:null;
  return rows.map(r=>Number.isFinite(baseline)&&Number.isFinite(r[key])?r[key]-baseline:null);
}}

function drawCoeffPlot(data, divId, yKey, title, yTitle) {{
  const groups = groupedByLabelPolar(data);
  const traces = [];
  Object.keys(groups).forEach(key => {{
    const g = groups[key];
    traces.push({{
      x:g.map(r=>r.alpha), y:g.map(r=>r[yKey]), mode:"lines+markers",
      name:key.replace("POLAR-", "P"),
      text:g.map(r=>`Configuration: ${{r.case_label}}<br>POLAR: ${{r.polar}}<br>Case: ${{r.case}}<br>Alpha: ${{fmt(r.alpha,3)}}<br>${{yTitle}}: ${{fmt(r[yKey],6)}}<br>Status: ${{r.status}}<extra></extra>`),
      hovertemplate:"%{{text}}"
    }});
  }});
  Plotly.newPlot(divId, traces, {{title:title, xaxis:{{title:"Alpha [deg]"}}, yaxis:{{title:yTitle}}, margin:{{l:60,r:30,t:50,b:50}}}}, {{responsive:true}});
}}
function drawResidualSummaryPlot(data) {{
  const groups = groupedByLabelPolar(data); const traces = [];
  Object.keys(groups).forEach(key => {{
    const g = groups[key];
    traces.push({{
      x:g.map(r=>r.alpha), y:g.map(r=>r.worst_residual_final), mode:"lines+markers",
      name:key.replace("POLAR-", "P"),
      text:g.map(r=>`Configuration: ${{r.case_label}}<br>POLAR: ${{r.polar}}<br>Case: ${{r.case}}<br>Alpha: ${{fmt(r.alpha,3)}}<br>Worst residual: ${{Number(r.worst_residual_final).toExponential(3)}}<br>Equation: ${{r.worst_residual_eq}}<br>Status: ${{r.status}}<extra></extra>`),
      hovertemplate:"%{{text}}"
    }});
  }});
  Plotly.newPlot("residualSummaryPlot", traces, {{title:"Worst final residual", xaxis:{{title:"Alpha [deg]"}}, yaxis:sciAxis("Residual"), shapes:[{{type:"line",xref:"paper",x0:0,x1:1,y0:8e-4,y1:8e-4,line:{{dash:"dot"}}}},{{type:"line",xref:"paper",x0:0,x1:1,y0:8e-3,y1:8e-3,line:{{dash:"dot"}}}},{{type:"line",xref:"paper",x0:0,x1:1,y0:1e-2,y1:1e-2,line:{{dash:"dot"}}}}], margin:{{l:70,r:40,t:50,b:50}}}}, {{responsive:true}});
}}
function setResidualEquation(eq) {{
  document.getElementById("residualEquationSelect").value = eq;
  drawResidualEquationsPlot(getFilteredRows());
}}

function drawResidualEquationsPlot(data) {{
  const selectedResidual = document.getElementById("residualEquationSelect") ? document.getElementById("residualEquationSelect").value : "ALL";
  const colsToPlot = selectedResidual === "ALL" ? residualCols : [selectedResidual];
  const groups = groupedByLabelPolar(data); const traces = [];
  Object.keys(groups).forEach(key => {{
    const g = groups[key];
    colsToPlot.forEach(col => traces.push({{
      x:g.map(r=>r.alpha), y:g.map(r=>r[col + "_final"]), mode:"lines+markers",
      name:key.replace("POLAR-", "P") + (selectedResidual === "ALL" ? " " + col : ""),
      text:g.map(r=>`Configuration: ${{r.case_label}}<br>POLAR: ${{r.polar}}<br>Case: ${{r.case}}<br>Alpha: ${{fmt(r.alpha,3)}}<br>${{col}}: ${{Number(r[col + "_final"]).toExponential(3)}}<extra></extra>`),
      hovertemplate:"%{{text}}"
    }}));
  }});
  const plotTitle = selectedResidual === "ALL" ? "Final residuals by equation" : `Final residual: ${{selectedResidual}}`;
  Plotly.newPlot("residualEquationsPlot", traces, {{title:plotTitle, xaxis:{{title:"Alpha [deg]"}}, yaxis:sciAxis("Residual"), margin:{{l:70,r:40,t:50,b:50}}}}, {{responsive:true}});
}}
function drawCpmaxPlot(data) {{
  const groups = groupedByLabelPolar(data); const traces = [];
  Object.keys(groups).forEach(key => {{ const g = groups[key]; traces.push({{x:g.map(r=>r.alpha), y:g.map(r=>r.cpmax_final), mode:"lines+markers", name:key}}); }});
  Plotly.newPlot("cpmaxPlot", traces, {{title:"Final cp-max", xaxis:{{title:"Alpha [deg]"}}, yaxis:{{title:"cp-max"}}, margin:{{l:60,r:30,t:50,b:50}}}}, {{responsive:true}});
}}

function drawSelectedHistory() {{
  const caseKey = document.getElementById("caseHistoryFilter").value;
  const hist = historyData[caseKey] || [];
  if (hist.length === 0) {{ Plotly.newPlot("selectedResidualHistory", [], {{title:"No history found"}}); Plotly.newPlot("selectedAeroHistory", [], {{title:"No history found"}}); Plotly.newPlot("selectedCpHistory", [], {{title:"No history found"}}); return; }}
  const x = hist.map(r => r.iter);
  Plotly.newPlot("selectedResidualHistory", residualCols.map(col => ({{x:x, y:hist.map(r => Math.abs(r[col])), mode:"lines", name:col}})), {{title:"Residual history", xaxis:{{title:"Iteration"}}, yaxis:sciAxis("Residual"), margin:{{l:70,r:40,t:50,b:50}}}}, {{responsive:true}});
  Plotly.newPlot("selectedAeroHistory", [{{x:x,y:hist.map(r=>r.clzb),mode:"lines",name:"CLZB"}},{{x:x,y:hist.map(r=>r.cdxb),mode:"lines",name:"CDXB"}},{{x:x,y:hist.map(r=>r.cmyb),mode:"lines",name:"CMYB"}}], {{title:"Aero history", xaxis:{{title:"Iteration"}}, yaxis:{{title:"Coefficient"}}, margin:{{l:70,r:40,t:50,b:50}}}}, {{responsive:true}});
  Plotly.newPlot("selectedCpHistory", [{{x:x,y:hist.map(r=>r["cp-max"]),mode:"lines",name:"cp-max"}},{{x:x,y:hist.map(r=>r["tstep-ave"]),mode:"lines",name:"tstep-ave",yaxis:"y2"}}], {{title:"cp-max and time-step history", xaxis:{{title:"Iteration"}}, yaxis:{{title:"cp-max"}}, yaxis2:{{title:"tstep-ave",overlaying:"y",side:"right"}}, margin:{{l:70,r:70,t:50,b:50}}}}, {{responsive:true}});
}}

function drawDeflections() {{
  const keys = ["RUD1","RUD2","RUD3","RUD4","ELV1","ELV2","ELV3","ELV4","AIL1","AIL2","AIL3","AIL4","FLP1","FLP2","FLP3","FLP4"];
  document.querySelector("#deflectionTable tbody").innerHTML = filteredSummaries().map(s => {{
    const d = s.deflections || {{}};
    const cells = keys.map(k => `<td>${{d[k] ?? "-"}}</td>`).join("");
    return `<tr><td>${{s.case_label}}</td><td>${{s.polar}}</td>${{cells}}</tr>`;
  }}).join("");
}}

function drawClassificationTable(data) {{
  document.querySelector("#classificationTable tbody").innerHTML = data.map(r => `
    <tr>
      <td>${{r.case_label ?? ""}}</td><td>${{r.polar ?? ""}}</td><td>${{r.case ?? ""}}</td>
      <td>${{fmt(r.alpha,2)}}</td><td>${{fmt(r.beta,2)}}</td>
      <td class="${{assessmentClass(r.numerical_status)}}">${{r.numerical_status ?? ""}}</td>
      <td class="${{assessmentClass(r.aerodynamic_status)}}">${{r.aerodynamic_status ?? ""}}</td>
      <td class="${{assessmentClass(r.neighbor_status)}}">${{r.neighbor_status ?? ""}}</td>
      <td class="status-${{r.status}}">${{r.status ?? ""}}</td>
      <td>${{(r.reasons || []).join(" | ")}}</td>
    </tr>`).join("");
}}

function drawMeshQuality() {{
  document.querySelector("#meshQualityTable tbody").innerHTML = filteredSummaries().map(s => {{
    const m = s.mesh || {{}};
    const warningText = (m.warnings || []).join(" | ");
    return `<tr>
      <td>${{s.case_label}}</td><td>${{s.polar}}</td>
      <td>${{m.cell_count ?? ""}}</td><td>${{m.node_count ?? ""}}</td>
      <td>${{fmt(m.min_orthogonal_quality,6)}}</td><td>${{fmt(m.max_aspect_ratio,4)}}</td>
      <td>${{m.negative_volume_warnings ?? 0}}</td>
      <td class="${{assessmentClass(m.status)}}">${{m.status ?? "UNKNOWN"}}</td>
      <td>${{warningText}}</td>
    </tr>`;
  }}).join("");
}}

function drawOutlierDiagnostics(data) {{
  const checked = data.filter(r => r.outlier_score !== null && r.outlier_score !== undefined);
  const groups = groupedByLabelPolar(checked);
  const traces = [];
  Object.keys(groups).forEach(key => {{
    const g = groups[key];
    traces.push({{
      x: g.map(r => r.neighbor_sweep === "BETA" ? r.beta : r.alpha),
      y: g.map(r => r.outlier_score),
      mode: "lines+markers",
      name: key.replace("POLAR-", "P"),
      marker: {{color: g.map(r => r.neighbor_status === "OUTLIER" ? "#d98200" : "#1a8f3a")}},
      text: g.map(r => `Configuration: ${{r.case_label}}<br>POLAR: ${{r.polar}}<br>Case: ${{r.case}}<br>Sweep: ${{r.neighbor_sweep}}<br>Score: ${{fmt(r.outlier_score,3)}}<br>ΔCL: ${{fmt(r.neighbor_dcl,5)}}<br>ΔCD: ${{fmt(r.neighbor_dcd,5)}}<br>ΔCM: ${{fmt(r.neighbor_dcm,5)}}<br>Status: ${{r.neighbor_status}}<extra></extra>`),
      hovertemplate: "%{{text}}"
    }});
  }});
  Plotly.newPlot("outlierScorePlot", traces, {{
    title: "Neighbor interpolation outlier score",
    xaxis: {{title: "Sweep value [deg]"}},
    yaxis: {{title: "Normalized outlier score"}},
    shapes: [{{type:"line", xref:"paper", x0:0, x1:1, y0:1, y1:1, line:{{dash:"dot", color:"#d98200"}}}}],
    margin: {{l:70,r:30,t:50,b:50}}
  }}, {{responsive:true}});

  const outliers = data.filter(r => r.neighbor_status === "OUTLIER");
  document.querySelector("#outlierTable tbody").innerHTML = outliers.map(r => `
    <tr><td>${{r.case_label}}</td><td>${{r.polar}}</td><td>${{r.case}}</td><td>${{r.neighbor_sweep}}</td><td>${{r.outlier_source ?? ""}}</td>
    <td>${{fmt(r.alpha,2)}}</td><td>${{fmt(r.beta,2)}}</td><td>${{fmt(r.outlier_score,3)}}</td>
    <td>${{fmt(r.neighbor_dcl,5)}}</td><td>${{fmt(r.neighbor_dcd,5)}}</td><td>${{fmt(r.neighbor_dcm,5)}}</td>
    <td class="${{assessmentClass(r.neighbor_status)}}">${{r.neighbor_status}}</td><td>${{(r.neighbor_reasons || []).join(" | ")}}</td></tr>`).join("");
  if (!outliers.length) document.querySelector("#outlierTable tbody").innerHTML = `<tr><td colspan="13">No neighbor outliers detected for the current filters.</td></tr>`;
}}

function drawDiagnostics() {{
  document.querySelector("#diagTable tbody").innerHTML = filteredSummaries().map(s => {{ const d = s.diagnostics || {{}}; return `<tr><td>${{s.case_label}}</td><td>${{s.polar}}</td><td>${{d.history_headers_found ?? ""}}</td><td>${{d.history_rows_found ?? ""}}</td><td>${{d.total_expected_iters ?? ""}}</td><td>${{d.rows_consumed_by_infout ?? ""}}</td><td>${{d.rows_remaining_after_split ?? ""}}</td><td>${{d.log_path ?? ""}}</td></tr>`; }}).join("");
}}

function drawTable(data) {{
  document.querySelector("#caseTable tbody").innerHTML = data.map(r => `<tr>
    <td>${{r.case_label ?? ""}}</td><td>${{r.polar ?? ""}}</td><td>${{r.case ?? ""}}</td><td>${{r.name ?? ""}}</td>
    <td>${{fmt(r.alpha,2)}}</td><td>${{fmt(r.beta,2)}}</td>
    <td class="${{assessmentClass(r.numerical_status)}}">${{r.numerical_status ?? ""}}</td>
    <td class="${{assessmentClass(r.aerodynamic_status)}}">${{r.aerodynamic_status ?? ""}}</td>
    <td class="${{assessmentClass(r.neighbor_status)}}">${{r.neighbor_status ?? ""}}</td>
    <td class="status-${{r.status}}">${{r.status}}</td><td>${{r.quality ?? ""}}</td>
    <td>${{r.actual_iters ?? ""}}</td><td>${{fmt(r["continuity_final"],3)}}</td><td>${{fmt(r["x-velocity_final"],3)}}</td>
    <td>${{fmt(r["y-velocity_final"],3)}}</td><td>${{fmt(r["z-velocity_final"],3)}}</td><td>${{fmt(r["energy_final"],3)}}</td>
    <td>${{fmt(r["nut_final"],3)}}</td><td>${{fmt(r.clzb_final,5)}}</td><td>${{fmt(r.cdxb_final,5)}}</td>
    <td>${{fmt(r.cmyb_final,5)}}</td><td>${{fmt(r.cpmax_final,4)}}</td><td>${{(r.reasons || []).join(" | ")}}</td>
  </tr>`).join("");
}}


function refreshAll() {{
  const data = getFilteredRows();
  updateMomentRefInfo();
  drawKpis(data); drawAssessmentKpis(data); drawStatusMap(data);
  drawAdfCoeffPlot(); drawAdfXY("dragPolarPlot", "CLS", "CDS", "CDS vs CLS"); drawAdfXY("cmClPlot", "CLS", "CMS25", "CMS25 vs CLS"); drawAdfXY("ldPlot", "CLS", "CDS", "L/D vs CLS", r => r.CLS / r.CDS, "L/D"); drawSmPlots(); drawDragRisePlots();
   drawResidualEquationsPlot(data); drawCpmaxPlot(data);
  drawClassificationTable(data); drawOutlierDiagnostics(data); drawMeshQuality(); drawDeflections(); drawDiagnostics(); drawTable(data);
}}


// =============================================================================
// V19 ENGINEERING UI AND ANALYSIS ENHANCEMENTS
// =============================================================================
let currentSectionName = "overview";
let focusMode = "all";
let comparisonMetric = "CL";
let selectedPlotId = "adfCoeffPlot";
const CURVE_COLORS = [
  "#2f5f98", "#c44e52", "#55a868", "#8172b2", "#cc8963", "#4c4c4c",
  "#64b5cd", "#dd8452", "#937860", "#da8bc3", "#8c8c8c", "#ccb974",
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"
];
const CONFIG_COLORS = CURVE_COLORS;
const POLAR_DASHES = ["solid", "dash", "dot", "dashdot", "longdash", "longdashdot"];
const MARKER_SYMBOLS = ["circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "star", "hexagon", "pentagon"];

const baseShowSection = showSection;
showSection = function(sectionName, button) {{
  currentSectionName = sectionName;
  baseShowSection(sectionName, button);
  saveLastState();
}};

function toggleAdvancedPanel() {{
  document.getElementById("advancedPanel").classList.toggle("open");
}}

function numericInput(id, fallback) {{
  const el = document.getElementById(id);
  const v = el ? Number(el.value) : fallback;
  return Number.isFinite(v) ? v : fallback;
}}

function currentThresholds() {{
  return {{
    residualSuccessful: numericInput("thrResSuccessful", JAMAL_DEFAULTS.residualSuccessful),
    residualAcceptable: numericInput("thrResAcceptable", JAMAL_DEFAULTS.residualAcceptable),
    residualBad: numericInput("thrResBad", JAMAL_DEFAULTS.residualBad),
    clRange: numericInput("thrClRange", JAMAL_DEFAULTS.clRange),
    cdRange: numericInput("thrCdRange", JAMAL_DEFAULTS.cdRange),
    cmRange: numericInput("thrCmRange", JAMAL_DEFAULTS.cmRange),
    cpRange: numericInput("thrCpRange", JAMAL_DEFAULTS.cpRange),
    neighborCl: numericInput("thrNeighborCl", JAMAL_DEFAULTS.neighborCl),
    neighborCd: numericInput("thrNeighborCd", JAMAL_DEFAULTS.neighborCd),
    neighborCm: numericInput("thrNeighborCm", JAMAL_DEFAULTS.neighborCm)
  }};
}}

function resetThresholds() {{
  const mapping = {{thrResSuccessful:"residualSuccessful",thrResAcceptable:"residualAcceptable",thrResBad:"residualBad",thrClRange:"clRange",thrCdRange:"cdRange",thrCmRange:"cmRange",thrCpRange:"cpRange",thrNeighborCl:"neighborCl",thrNeighborCd:"neighborCd",thrNeighborCm:"neighborCm"}};
  Object.entries(mapping).forEach(([id,key]) => document.getElementById(id).value = JAMAL_DEFAULTS[key]);
  refreshAll();
}}

function evaluateCase(raw) {{
  const r = Object.assign({{}}, raw);
  const t = currentThresholds();
  let numerical = "SUCCESSFUL", numericalReasons = [];
  if (!r.has_history) {{ numerical = "DIVERGED"; numericalReasons = ["No iteration history found"]; }}
  else if ((r.worst_residual_final ?? 0) >= t.residualBad || (r.worst_residual_p95 ?? 0) >= t.residualBad || (r.residual_trend_ratio ?? 0) > 10) {{
    numerical = "DIVERGED";
    if ((r.worst_residual_final ?? 0) >= t.residualBad) numericalReasons.push(`final residual ${{fmt(r.worst_residual_final,3)}} ≥ ${{t.residualBad.toExponential(1)}}`);
    if ((r.worst_residual_p95 ?? 0) >= t.residualBad) numericalReasons.push(`p95 residual ${{fmt(r.worst_residual_p95,3)}} ≥ ${{t.residualBad.toExponential(1)}}`);
    if ((r.residual_trend_ratio ?? 0) > 10) numericalReasons.push(`residual grew ${{fmt(r.residual_trend_ratio,2)}}× near the end`);
  }} else if ((r.worst_residual_final ?? 0) >= t.residualAcceptable || ((r.actual_iters ?? 0) < 0.8 * (r.expected_iters ?? 0))) {{
    numerical = "WARNING";
    if ((r.worst_residual_final ?? 0) >= t.residualAcceptable) numericalReasons.push(`final residual ${{fmt(r.worst_residual_final,3)}} is above acceptable threshold`);
    if ((r.actual_iters ?? 0) < 0.8 * (r.expected_iters ?? 0)) numericalReasons.push(`actual iterations below 80% of expected`);
  }} else if ((r.worst_residual_final ?? 0) >= t.residualSuccessful) {{
    numerical = "ACCEPTABLE"; numericalReasons = ["residual below acceptable threshold"];
  }} else {{ numericalReasons = ["residual below successful threshold"]; }}

  const aeroReasons = [];
  if ((r.clzb_range ?? 0) > t.clRange) aeroReasons.push(`CLZB range ${{fmt(r.clzb_range,5)}} > ${{fmt(t.clRange,5)}}`);
  if ((r.cdxb_range ?? 0) > t.cdRange) aeroReasons.push(`CDXB range ${{fmt(r.cdxb_range,6)}} > ${{fmt(t.cdRange,6)}}`);
  if ((r.cmyb_range ?? 0) > t.cmRange) aeroReasons.push(`CMYB range ${{fmt(r.cmyb_range,6)}} > ${{fmt(t.cmRange,6)}}`);
  if ((r.cpmax_window_range ?? 0) > t.cpRange) aeroReasons.push(`cp-max range ${{fmt(r.cpmax_window_range,4)}} > ${{fmt(t.cpRange,4)}}`);
  const aerodynamic = aeroReasons.length ? "WARNING" : "PASS";

  const neighborReasons = [];
  if (r.neighbor_dcl !== null && r.neighbor_dcl !== undefined && r.neighbor_dcl > t.neighborCl) neighborReasons.push(`neighbor ΔCL ${{fmt(r.neighbor_dcl,5)}}`);
  if (r.neighbor_dcd !== null && r.neighbor_dcd !== undefined && r.neighbor_dcd > t.neighborCd) neighborReasons.push(`neighbor ΔCD ${{fmt(r.neighbor_dcd,5)}}`);
  if (r.neighbor_dcm !== null && r.neighbor_dcm !== undefined && r.neighbor_dcm > t.neighborCm) neighborReasons.push(`neighbor ΔCM ${{fmt(r.neighbor_dcm,5)}}`);
  const neighbor = (r.outlier_score === null || r.outlier_score === undefined) ? "NOT_CHECKED" : (neighborReasons.length ? "OUTLIER" : "PASS");
  const scores = [];
  if (r.neighbor_dcl !== null && r.neighbor_dcl !== undefined) scores.push(r.neighbor_dcl / t.neighborCl);
  if (r.neighbor_dcd !== null && r.neighbor_dcd !== undefined) scores.push(r.neighbor_dcd / t.neighborCd);
  if (r.neighbor_dcm !== null && r.neighbor_dcm !== undefined) scores.push(r.neighbor_dcm / t.neighborCm);

  let overall = "CONVERGED", quality = "successful";
  if (numerical === "DIVERGED") {{ overall = "DIVERGED"; quality = "bad"; }}
  else if (numerical === "WARNING" || aerodynamic === "WARNING" || neighbor === "OUTLIER") {{ overall = "SUSPICIOUS"; quality = "warning"; }}
  else if (numerical === "ACCEPTABLE") {{ overall = "ACCEPTABLE"; quality = "acceptable"; }}

  r.numerical_status = numerical; r.numerical_reasons = numericalReasons;
  r.aerodynamic_status = aerodynamic; r.aerodynamic_reasons = aeroReasons.length ? aeroReasons : ["final-window monitors stable"];
  r.neighbor_status = neighbor; r.neighbor_reasons = neighborReasons.length ? neighborReasons : (neighbor === "PASS" ? ["consistent with adjacent sweep points"] : []);
  r.outlier_score = scores.length ? Math.max(...scores) : null;
  r.status = overall; r.overall_status = overall; r.quality = quality;
  r.reasons = [...r.numerical_reasons, ...r.aerodynamic_reasons, ...r.neighbor_reasons];
  return r;
}}

function evaluatedRows() {{ return convRows.map(evaluateCase); }}

getFilteredRows = function() {{
  const label = document.getElementById("labelFilter").value;
  const polar = document.getElementById("polarFilter").value;
  const status = document.getElementById("statusFilter").value;
  return evaluatedRows().filter(r => {{
    if (label !== "ALL" && r.case_label !== label) return false;
    if (polar !== "ALL" && r.polar !== polar) return false;
    if (status !== "ALL" && r.status !== status) return false;
    if (!curveAllowed(r.case_label, r.polar)) return false;
    return true;
  }});
}};

function selectedCaseRow() {{
  const key = document.getElementById("caseHistoryFilter")?.value;
  return evaluatedRows().find(r => r.case_key === key) || evaluatedRows()[0] || null;
}}

function selectedCurveKey() {{
  const r = selectedCaseRow();
  return r ? `${{r.case_label}}|${{r.polar}}` : null;
}}

function baselineLabel() {{
  const first = adfData.curves && adfData.curves.length ? adfData.curves[0].case_label : null;
  return first;
}}

function curveAllowed(label, polar) {{
  if (focusMode === "all") return true;
  if (focusMode === "none") return false;
  if (focusMode === "baseline") return label === baselineLabel();
  const selected = selectedCurveKey();
  if (!selected) return true;
  const key = `${{label}}|${{polar}}`;
  if (focusMode === "selected") return key === selected;
  if (focusMode === "selected_baseline") return key === selected || label === baselineLabel();
  return true;
}}

function setFocusMode(mode) {{
  focusMode = mode;
  const labels = {{all:"All curves", baseline:"Baseline only", selected:"Selected curve only", selected_baseline:"Selected + baseline", none:"All curves hidden"}};
  document.getElementById("focusModeLabel").textContent = labels[mode] || mode;
  refreshAll();
}}

const originalFilteredAdfCurvesV19 = filteredAdfCurves;
filteredAdfCurves = function() {{ return originalFilteredAdfCurvesV19().filter(c => curveAllowed(c.case_label, c.polar)); }};
const originalFilteredDragRiseCurvesV19 = filteredDragRiseCurves;
filteredDragRiseCurves = function() {{
  const curves = originalFilteredDragRiseCurvesV19();
  if (focusMode === "all") return curves;
  if (focusMode === "none") return [];
  if (focusMode === "baseline") return curves.filter(c => c.case_label === baselineLabel());
  const r = selectedCaseRow();
  if (!r) return curves;
  return curves.filter(c => c.case_label === r.case_label || (focusMode === "selected_baseline" && c.case_label === baselineLabel()));
}};

function curveStyleCatalog() {{
  const items = [];
  (adfData.curves || []).forEach(c => items.push({{label:c.case_label, polar:c.polar}}));
  (convRows || []).forEach(r => items.push({{label:r.case_label, polar:r.polar}}));
  (dragRiseData.curves || []).forEach(c => items.push({{label:c.case_label, polar:c.cls_label || c.file_name || "Drag rise"}}));
  const labels = unique(items.map(x => x.label));
  const polars = unique(items.map(x => x.polar));
  const curves = unique(items.map(x => `${{x.label}}|${{x.polar}}`));
  return {{labels, polars, curves}};
}}

function currentCurveStyleSettings() {{
  const catalog = curveStyleCatalog();
  const mode = document.getElementById("curveStyleMode")?.value || "auto";
  if (mode === "colors_only") return {{color:"every_curve", line:"solid", marker:"every_curve"}};
  if (mode === "color_polar") return {{color:"polar", line:"configuration", marker:"configuration"}};
  if (mode === "color_configuration") return {{color:"configuration", line:"polar", marker:"polar"}};
  if (mode === "custom") return {{
    color:document.getElementById("customColorMapping")?.value || "every_curve",
    line:document.getElementById("customLineMapping")?.value || "solid",
    marker:document.getElementById("customMarkerMapping")?.value || "none"
  }};
  // Auto
  if (catalog.labels.length <= 1) return {{color:"polar", line:"solid", marker:"polar"}};
  return {{color:"polar", line:"configuration", marker:"configuration"}};
}}

function updateCurveStyleControls() {{
  const custom = document.getElementById("curveStyleMode")?.value === "custom";
  document.getElementById("curveStyleCustom")?.classList.toggle("open", custom);
}}

function mappedIndex(mapping, label, polar, catalog) {{
  if (mapping === "configuration") return Math.max(0, catalog.labels.indexOf(label));
  if (mapping === "polar") return Math.max(0, catalog.polars.indexOf(polar));
  if (mapping === "every_curve") return Math.max(0, catalog.curves.indexOf(`${{label}}|${{polar}}`));
  return 0;
}}

function traceStyle(label, polar) {{
  const catalog = curveStyleCatalog();
  const settings = currentCurveStyleSettings();
  const colorIndex = mappedIndex(settings.color, label, polar, catalog);
  const dashIndex = mappedIndex(settings.line, label, polar, catalog);
  const markerIndex = mappedIndex(settings.marker, label, polar, catalog);
  return {{
    color: CURVE_COLORS[colorIndex % CURVE_COLORS.length],
    dash: settings.line === "solid" ? "solid" : POLAR_DASHES[dashIndex % POLAR_DASHES.length],
    symbol: settings.marker === "none" ? "circle" : MARKER_SYMBOLS[markerIndex % MARKER_SYMBOLS.length],
    markerSize: settings.marker === "none" ? 0 : 6,
    mode: settings.marker === "none" ? "lines" : "lines+markers"
  }};
}}

function formatScientificInteger(value) {{
  const v = Number(value);
  if (!Number.isFinite(v)) return "n/a";
  const parts = v.toExponential(0).split("e");
  const exponent = Number(parts[1]);
  const sign = exponent >= 0 ? "+" : "-";
  return `${{parts[0]}}E${{sign}}${{String(Math.abs(exponent)).padStart(2,"0")}}`;
}}

function conditionText(curves) {{
  const mach = [], re = [];
  (curves || []).forEach(c => (c.rows || []).forEach(r => {{ if (Number.isFinite(Number(r.MACH))) mach.push(Number(r.MACH)); if (Number.isFinite(Number(r.REYNOLDS))) re.push(Number(r.REYNOLDS)); }}));
  function machText(values) {{
    if (!values.length) return "M=n/a";
    const mn=Math.min(...values), mx=Math.max(...values);
    return Math.abs(mx-mn)<1e-6 ? `M=${{mn.toFixed(2)}}` : `M=${{mn.toFixed(2)}}–${{mx.toFixed(2)}}`;
  }}
  function reynoldsText(values) {{
    if (!values.length) return "Re=n/a";
    const mn=Math.min(...values), mx=Math.max(...values);
    return Math.abs(mx-mn)<=Math.max(1,1e-6*Math.abs((mn+mx)/2)) ? `Re=${{formatScientificInteger(mn)}}` : `Re=${{formatScientificInteger(mn)}}–${{formatScientificInteger(mx)}}`;
  }}
  return `${{machText(mach)}} · ${{reynoldsText(re)}}`;
}}

function fixedAngleText(curves) {{
  const labels=(curves||[]).map(c=>{{
    const fixed=c.sweep_var==='BETA'?'ALPHA':'BETA';
    const values=(c.rows||[]).map(r=>r[fixed]).filter(Number.isFinite);
    if(!values.length) return `${{fixed}}=n/a`;
    const low=Math.min(...values),high=Math.max(...values);
    const value=Math.abs(high-low)<1e-8?low.toFixed(2):`${{low.toFixed(2)}}–${{high.toFixed(2)}} (varies)`;
    return `${{fixed}}=${{value}}°`;
  }});
  const unique=[...new Set(labels)];
  if(unique.length<=1) return unique[0]||'';
  return labels.map((label,i)=>`${{curves[i].case_label}} ${{curves[i].polar}}: ${{label}}`).join(' · ');
}}

function referenceTitle() {{
  if (currentMomentReferenceMode() === "original") return "Xref = 25% MAC";
  const inp = currentMomentInputs();
  return inp.mode === "percent" ? `Xref = ${{Number(inp.xPercent).toFixed(0)}}% MAC` : `Xref = ${{Number(inp.xAbs).toFixed(0)}} m`;
}}

function axisName(axis) {{ return axis === "B" ? "Body-axis" : (axis === "W" ? "Wind-axis" : "Stability-axis"); }}
coefficientDisplayName = function(axis, baseCoeff) {{ return baseCoeff + axis; }};

function bindPlotInteractions(divId) {{
  const div = document.getElementById(divId);
  if (!div || typeof div.on !== "function") return;
  try {{ div.removeAllListeners("plotly_click"); }} catch (err) {{}}
  div.on("plotly_click", ev => {{
    const cd = ev?.points?.[0]?.customdata;
    if (!cd || !cd.case_label || !cd.polar) return;
    const candidates = evaluatedRows().filter(r => r.case_label === cd.case_label && r.polar === cd.polar);
    if (!candidates.length) return;
    let best = candidates[0], bestDist = Infinity;
    candidates.forEach(r => {{
      const da = Math.abs(Number(r.alpha ?? 0) - Number(cd.alpha ?? 0));
      const db = Math.abs(Number(r.beta ?? 0) - Number(cd.beta ?? 0));
      if (da + db < bestDist) {{ best = r; bestDist = da + db; }}
    }});
    document.getElementById("caseHistoryFilter").value = best.case_key;
    drawSelectedHistory();
    drawSelectedCaseSummary();
    saveLastState();
  }});
  div.onclick = () => {{ selectedPlotId = divId; const sel = document.getElementById("exportPlotSelect"); if (sel && [...sel.options].some(o => o.value === divId)) sel.value = divId; }};
}}

const plotlyWithInteraction = Plotly.newPlot.bind(Plotly);
Plotly.newPlot = function(gd, data, layout, config) {{
  const id = typeof gd === "string" ? gd : gd.id;
  return plotlyWithInteraction(gd, data, layout, config).then(result => {{ bindPlotInteractions(id); return result; }});
}};

function customDataForCurve(c, row) {{ return {{case_label:c.case_label, polar:c.polar, alpha:row.ALPHA, beta:row.BETA}}; }}

function currentAxisKeys() {{
  const axis = document.getElementById("stabilityAxisSelect").value || "S";
  return {{axis, cl:`CL${{axis}}`, cd:`CD${{axis}}`, cm:`CM${{axis}}25`, cmLabel:`CM${{axis}}`, clLabel:`CL${{axis}}`, cdLabel:`CD${{axis}}`}};
}}

drawAdfCoeffPlot = function() {{
  drawStandardAeroPlots();
  saveLastState();
}};

function coefficientPlotSpecs(axis, abscissa) {{
  const xKey=["CL","CY"].includes(abscissa)?coefficientKey(axis,abscissa):abscissa;
  const specs=["CL","CD","CY","CM","CR","CN"]
    .filter(base=>base!==abscissa)
    .map(base=>({{base,xKey,yKey:coefficientKey(axis,base),label:coefficientDisplayName(axis,base),
      id:base==="CL"?"adfCoeffPlot":`coefficient${{base}}Plot`,heading:`coefficient${{base}}Heading`}}));
  specs.push({{base:"LD",xKey,yKey:coefficientKey(axis,"CD"),clKey:coefficientKey(axis,"CL"),
    label:"L/D",id:"ldPlot",heading:"ldHeading"}});
  return specs;
}}

function coefficientPlotPoints(rows, spec) {{
  return rows.filter(r=>Number.isFinite(r[spec.xKey])).map(r=>{{
    const value=spec.base==="LD"
      ? (Number.isFinite(r[spec.clKey])&&Number.isFinite(r[spec.yKey])&&r[spec.yKey]!==0?r[spec.clKey]/r[spec.yKey]:null)
      : r[spec.yKey];
    return {{row:r,x:r[spec.xKey],y:Number.isFinite(value)?value:null}};
  }}).sort((a,b)=>a.x-b.x);
}}

function drawStandardAeroPlots() {{
  const missingReference=filteredAdfCurves().some(c=>{{const meta=summaryFor(c.case_label,c.polar)?.meta;return !meta||!Number.isFinite(meta.cref)||!Number.isFinite(meta.bref);}});
  const referenceSelect=document.getElementById('momentReferenceSelect');
  referenceSelect.querySelector('option[value="user"]').disabled=missingReference;
  if(missingReference) referenceSelect.value='original';
  const curves = currentAdfCurves();
  const axis=document.getElementById("stabilityAxisSelect").value||"S";
  const abscissa=document.getElementById("coefficientAbscissa").value;
  const specs=coefficientPlotSpecs(axis,abscissa);
  document.getElementById('coefficientConditions').textContent=[conditionText(curves).split(' · ').reverse().join(' · '),fixedAngleText(curves)].filter(Boolean).join(' · ');
  const reference=currentMomentReferenceMode()==='original'?'Original ADF reference':referenceTitle();
  document.getElementById('referenceSummary').textContent=reference+(missingReference?' · infout reference unavailable':'');
  document.getElementById('momentRowReference').textContent=`Moment coefficients · ${{reference}}`;
  const legend=document.getElementById('coefficientLegend');legend.replaceChildren();
  curves.forEach(c=>{{
    const key=JSON.stringify([c.case_label,c.polar]),style=traceStyle(c.case_label,c.polar);
    const button=document.createElement('button');button.setAttribute('aria-pressed',String(!hiddenCoefficientCurves.has(key)));
    const swatch=document.createElement('span');swatch.className='legend-swatch';swatch.style.borderColor=style.color;
    if(style.dash!=='solid') swatch.style.borderTopStyle='dashed';
    button.append(swatch,document.createTextNode(`${{c.case_label}} ${{c.polar.replace('POLAR-','P')}}`));
    button.onclick=()=>{{hiddenCoefficientCurves.has(key)?hiddenCoefficientCurves.delete(key):hiddenCoefficientCurves.add(key);drawAdfCoeffPlot();}};
    legend.appendChild(button);
  }});
  document.getElementById("coefficientCLCard").hidden=abscissa==="CL";
  document.getElementById("coefficientCYCard").hidden=abscissa==="CY";
  if(abscissa==="CL") Plotly.purge("adfCoeffPlot");
  if(abscissa==="CY") Plotly.purge("coefficientCYPlot");
  document.getElementById("aeroCoeffHeading").textContent=`${{axisName(axis)}} coefficients vs ${{specs[0].xKey}}`;
  specs.forEach(spec=>{{
    const title=conditionText(curves).split(" · ").reverse().join(" · ")+
      (["CM","CR","CN"].includes(spec.base)?` · ${{referenceTitle()}}`:"");
    const heading=document.getElementById(spec.heading);
    heading.replaceChildren(...title.split(/([+−–-]?\\d[\\d.Ee+−–%-]*)/g).filter(Boolean).map(part=>{{
      if(!/\\d/.test(part)) return document.createTextNode(part);
      const span=document.createElement("span");span.className="numeric-value";span.textContent=part;return span;
    }}));
    const traces = [];
    curves.forEach(c => {{
      const points=coefficientPlotPoints(c.rows,spec);
      if (!points.some(p=>p.y!==null)) return;
      const st = traceStyle(c.case_label,c.polar);
      traces.push({{x:points.map(p=>p.x),y:points.map(p=>p.y),connectgaps:false,mode:st.mode,visible:hiddenCoefficientCurves.has(JSON.stringify([c.case_label,c.polar]))?'legendonly':true,
        name:`${{c.case_label}} ${{c.polar.replace("POLAR-","P")}}`,line:{{color:st.color,dash:st.dash}},
        marker:{{color:st.color,symbol:st.symbol,size:st.markerSize}},customdata:points.map(p=>customDataForCurve(c,p.row)),
        hovertemplate:`${{spec.xKey}}: %{{x:.5f}}<br>${{spec.label}}: %{{y:.6f}}<extra>%{{fullData.name}}</extra>`}});
    }});
    Plotly.newPlot(spec.id,traces,{{xaxis:{{title:["ALPHA","BETA"].includes(spec.xKey)?`${{spec.xKey}} [deg]`:spec.xKey}},
      yaxis:{{title:spec.label}},showlegend:false,legend:{{orientation:"h",y:-.22}},margin:{{l:65,r:20,t:15,b:55}}}},{{responsive:true}});
  }});
  populateExportPlots();
}}

const baseDrawSmPlotsV19 = drawSmPlots;
drawSmPlots = function() {{
  const curves = currentAdfCurves();
  const sm = filteredSmRows();
  const groups = {{}}; sm.forEach(r => {{if(!groups[r.curve_label]) groups[r.curve_label]=[]; groups[r.curve_label].push(r);}});
  const tracesCl=[], tracesSweep=[]; let sweepName="ALPHA";
  Object.keys(groups).forEach(key=>{{
    const g=groups[key]; const first=g[0]; const st=traceStyle(first.case_label,first.polar);
    const g1=g.slice().sort((a,b)=>a.CLS-b.CLS); const g2=g.slice().sort((a,b)=>(a.ALPHA??a.BETA)-(b.ALPHA??b.BETA)); if(g2.length) sweepName=g2[0].sweep_var||sweepName;
    tracesCl.push({{x:g1.map(r=>r.CLS),y:g1.map(r=>r.STATIC_MARGIN_PERCENT),mode:st.mode,name:key.replace("POLAR-","P"),line:{{color:st.color,dash:st.dash}},marker:{{color:st.color,symbol:st.symbol,size:st.markerSize}},customdata:g1.map(r=>({{case_label:r.case_label,polar:r.polar,alpha:r.ALPHA,beta:r.BETA}})),hovertemplate:"CLS: %{{x:.5f}}<br>SM: %{{y:.3f}}%<extra>%{{fullData.name}}</extra>"}});
    tracesSweep.push({{x:g2.map(r=>r.ALPHA!==undefined?r.ALPHA:r.BETA),y:g2.map(r=>r.STATIC_MARGIN_PERCENT),mode:st.mode,name:key.replace("POLAR-","P"),line:{{color:st.color,dash:st.dash}},marker:{{color:st.color,symbol:st.symbol,size:st.markerSize}},customdata:g2.map(r=>({{case_label:r.case_label,polar:r.polar,alpha:r.ALPHA,beta:r.BETA}})),hovertemplate:`${{sweepName}}: %{{x:.3f}}<br>SM: %{{y:.3f}}%<extra>%{{fullData.name}}</extra>`}});
  }});
  const cond=conditionText(curves), ref=referenceTitle();
  const base={{yaxis:{{title:"Static margin [% CREF]"}},margin:{{l:70,r:30,t:55,b:50}},shapes:[{{type:"line",xref:"paper",x0:0,x1:1,y0:0,y1:0,line:{{dash:"dash"}}}}]}};
  if(staticMarginYlim!==null&&!filterStaticMargin) base.yaxis.range=staticMarginYlim;
  Plotly.newPlot("smClPlot",tracesCl,Object.assign({{}},base,{{title:`Static margin vs CLS · ${{ref}} · ${{cond}}`,xaxis:{{title:"CLS"}}}}),{{responsive:true}});
  Plotly.newPlot("smSweepPlot",tracesSweep,Object.assign({{}},base,{{title:`Static margin vs ${{sweepName}} · ${{ref}} · ${{cond}}`,xaxis:{{title:`${{sweepName}} [deg]`}}}}),{{responsive:true}});
  document.getElementById("smClHeading").textContent="Static margin vs CLS";
  document.getElementById("smSweepHeading").textContent=`Static margin vs ${{sweepName}}`;
}};

function interpolationPoints(rows, xKey, yFunc) {{
  const buckets = new Map();
  (rows || []).forEach(r => {{
    const x = Number(r[xKey]);
    const y = Number(yFunc(r));
    if (!Number.isFinite(x) || !Number.isFinite(y)) return;
    const key = x.toPrecision(15);
    if (!buckets.has(key)) buckets.set(key, {{x:x, values:[]}});
    buckets.get(key).values.push(y);
  }});
  return [...buckets.values()]
    .map(p => ({{x:p.x, y:p.values.reduce((a,b)=>a+b,0)/p.values.length}}))
    .sort((a,b)=>a.x-b.x);
}}

function interpolateAt(rows, xKey, yFunc, x) {{
  const pts = interpolationPoints(rows, xKey, yFunc);
  if (!pts.length || !Number.isFinite(Number(x)) || x < pts[0].x || x > pts[pts.length-1].x) return null;
  for (let i=0; i<pts.length; i++) if (Math.abs(pts[i].x-x) < 1e-10) return pts[i].y;
  for (let i=0; i<pts.length-1; i++) {{
    if (pts[i].x <= x && x <= pts[i+1].x) {{
      const dx = pts[i+1].x - pts[i].x;
      if (Math.abs(dx) < 1e-14) return 0.5 * (pts[i].y + pts[i+1].y);
      return pts[i].y + (x-pts[i].x) * (pts[i+1].y-pts[i].y) / dx;
    }}
  }}
  return null;
}}

function setComparisonMetric(metric, button) {{
  comparisonMetric=metric;
  document.querySelectorAll("#comparisonMetricButtons button").forEach(b=>b.classList.remove("active"));
  if(button)button.classList.add("active");
  drawComparisonPlot();
}}

const MAX_COMPARISON_PAIRS = 8;
let comparisonPairCounter = 0;

function comparisonCurveOptions(selectedKey="") {{
  return adfData.curves.map(c => {{
    const key = `${{c.case_label}}|${{c.polar}}`;
    return `<option value="${{key}}" ${{key===selectedKey?"selected":""}}>${{c.curve_label}}</option>`;
  }}).join("");
}}

function defaultComparisonKeys(index=0) {{
  const curves = adfData.curves || [];
  if (!curves.length) return ["",""];
  const ref = curves[Math.min(index, curves.length-1)];
  const cmp = curves[Math.min(index+1, curves.length-1)];
  return [`${{ref.case_label}}|${{ref.polar}}`, `${{cmp.case_label}}|${{cmp.polar}}`];
}}

function addComparisonPair(referenceKey=null, comparisonKey=null) {{
  const container = document.getElementById("comparisonPairsContainer");
  if (!container || container.querySelectorAll(".comparison-pair").length >= MAX_COMPARISON_PAIRS) return;
  const idx = comparisonPairCounter++;
  const defaults = defaultComparisonKeys(container.querySelectorAll(".comparison-pair").length);
  const refKey = referenceKey || defaults[0];
  const cmpKey = comparisonKey || defaults[1];
  const row = document.createElement("div");
  row.className = "comparison-pair";
  row.dataset.pairId = String(idx);
  row.innerHTML = `
    <span class="pair-index">${{container.children.length+1}}</span>
    <label>Reference:</label><select class="comparison-reference" onchange="drawComparisonPlot();saveLastState();">${{comparisonCurveOptions(refKey)}}</select>
    <label>Comparison:</label><select class="comparison-curve" onchange="drawComparisonPlot();saveLastState();">${{comparisonCurveOptions(cmpKey)}}</select>
    <button onclick="removeComparisonPair(this)">Remove</button>`;
  container.appendChild(row);
  renumberComparisonPairs();
  drawComparisonPlot();
}}

function renumberComparisonPairs() {{
  document.querySelectorAll("#comparisonPairsContainer .comparison-pair").forEach((row,i) => {{
    const label = row.querySelector(".pair-index"); if (label) label.textContent = String(i+1);
  }});
}}

function removeComparisonPair(button) {{
  const row = button.closest(".comparison-pair"); if (row) row.remove();
  renumberComparisonPairs(); drawComparisonPlot(); saveLastState();
}}

function removeLastComparisonPair() {{
  const rows = document.querySelectorAll("#comparisonPairsContainer .comparison-pair");
  if (rows.length) rows[rows.length-1].remove();
  renumberComparisonPairs(); drawComparisonPlot(); saveLastState();
}}

function clearComparisonPairs() {{
  const container = document.getElementById("comparisonPairsContainer");
  if (container) container.innerHTML = "";
  drawComparisonPlot(); saveLastState();
}}

function getComparisonPairs() {{
  return [...document.querySelectorAll("#comparisonPairsContainer .comparison-pair")].map(row => ({{
    reference: row.querySelector(".comparison-reference")?.value || "",
    comparison: row.querySelector(".comparison-curve")?.value || ""
  }})).filter(p => p.reference && p.comparison);
}}

function setupComparisonControls(pairs=null) {{
  const container = document.getElementById("comparisonPairsContainer");
  if (!container) return;
  container.innerHTML = "";
  comparisonPairCounter = 0;
  const initial = Array.isArray(pairs) && pairs.length ? pairs : [{{reference:null,comparison:null}}];
  initial.slice(0,MAX_COMPARISON_PAIRS).forEach(p => addComparisonPair(p.reference,p.comparison));
}}

function curveFromKey(key) {{ return (adfData.curves || []).find(c=>`${{c.case_label}}|${{c.polar}}`===key); }}

function selectedComparisonSeries(ref,cmp,metric,smMap,axis,abscissa) {{
  const xKey=["CL","CY"].includes(abscissa)?coefficientKey(axis,abscissa):abscissa;
  const yKey=metric==="SM"?"STATIC_MARGIN_PERCENT":coefficientKey(axis,metric);
  const yFunc=metric==="LD"?r=>Number.isFinite(r["CL"+axis])&&Number.isFinite(r["CD"+axis])&&r["CD"+axis]!==0?r["CL"+axis]/r["CD"+axis]:null:r=>r[yKey];
  const source=c=>metric==="SM"?(smMap[`${{c.case_label}}|${{c.polar}}`]||[]):c.rows;
  const valid=rows=>rows.filter(r=>Number.isFinite(r[xKey])&&Number.isFinite(yFunc(r)));
  const reference=interpolationPoints(valid(source(ref)),xKey,yFunc);
  const comparison=valid(source(cmp)),xs=[],ys=[];
  reference.forEach(p=>{{const value=interpolateAt(comparison,xKey,yFunc,p.x);if(value!==null){{xs.push(p.x);ys.push(value-p.y);}}}});
  return {{xs,ys,yLabel:metric==="SM"?"ΔStatic margin [% CREF]":metric==="LD"?"ΔL/D":`Δ${{metric}}${{axis}}`,
    xAxisTitle:["ALPHA","BETA"].includes(xKey)?`${{xKey}} [deg]`:xKey,
    comparisonNote:`Comparison interpolated on reference ${{xKey}}; duplicate coordinates averaged`}};
}}

function comparisonSeries(ref, cmp, metric, smMap) {{
  let xs=[],ys=[],yLabel="",xAxisTitle="",comparisonNote="";
  if(metric==="CL"){{
    const xKey="ALPHA";
    const refRows=ref.rows.slice().sort((a,b)=>Number(a[xKey])-Number(b[xKey]));
    refRows.forEach(r=>{{
      const x=Number(r[xKey]), refY=Number(r.CLS), cmpY=interpolateAt(cmp.rows,xKey,q=>q.CLS,x);
      if(cmpY!==null&&Number.isFinite(refY)){{xs.push(x);ys.push(cmpY-refY);}}
    }});
    yLabel="ΔCLS"; xAxisTitle="ALPHA [deg]"; comparisonNote="Comparison interpolated on the reference ALPHA grid";
  }} else if(metric==="SM"){{
    const refSm=(smMap[`${{ref.case_label}}|${{ref.polar}}`]||[]).slice().sort((a,b)=>Number(a.CLS)-Number(b.CLS));
    const cmpSm=smMap[`${{cmp.case_label}}|${{cmp.polar}}`]||[];
    refSm.forEach(r=>{{
      const cl=Number(r.CLS), refY=Number(r.STATIC_MARGIN_PERCENT), cmpY=interpolateAt(cmpSm,"CLS",q=>q.STATIC_MARGIN_PERCENT,cl);
      if(cmpY!==null&&Number.isFinite(refY)){{xs.push(cl);ys.push(cmpY-refY);}}
    }});
    yLabel="ΔStatic margin [% CREF]"; xAxisTitle="CLS"; comparisonNote="Comparison interpolated on the reference CLS grid";
  }} else {{
    const refRows=ref.rows.slice().filter(r=>Number.isFinite(Number(r.CLS))).sort((a,b)=>Number(a.CLS)-Number(b.CLS));
    let yFuncRef,yFuncCmp;
    if(metric==="CD"){{yFuncRef=r=>r.CDS;yFuncCmp=r=>r.CDS;yLabel="ΔCDS";}}
    if(metric==="CM"){{yFuncRef=r=>r.CMS25;yFuncCmp=r=>r.CMS25;yLabel="ΔCMS";}}
    if(metric==="LD"){{yFuncRef=r=>Number(r.CLS)/Number(r.CDS);yFuncCmp=yFuncRef;yLabel="ΔL/D";}}
    refRows.forEach(r=>{{
      const cl=Number(r.CLS), refY=Number(yFuncRef(r)), cmpY=interpolateAt(cmp.rows,"CLS",yFuncCmp,cl);
      if(cmpY!==null&&Number.isFinite(refY)){{xs.push(cl);ys.push(cmpY-refY);}}
    }});
    xAxisTitle="CLS"; comparisonNote="Comparison interpolated on the reference CLS grid";
  }}
  return {{xs,ys,yLabel,xAxisTitle,comparisonNote}};
}}

function drawComparisonPlot() {{
  const pairs=getComparisonPairs();
  if(!pairs.length){{Plotly.newPlot("comparisonPlot",[],{{title:"Add one or more comparison pairs"}});return;}}
  const smMap={{}};
  currentStaticMarginRows().forEach(r=>{{const k=`${{r.case_label}}|${{r.polar}}`;if(!smMap[k])smMap[k]=[];smMap[k].push(r);}});
  const traces=[], allCurves=[];
  let yLabel="",xAxisTitle="";
  pairs.forEach((pair,index)=>{{
    const transform=c=>{{
      if(!c||currentMomentReferenceMode()==="original")return c;
      const meta=summaryFor(c.case_label,c.polar)?.meta;
      return {{...c,rows:c.rows.map(r=>shiftedRow(r,meta))}};
    }};
    const ref=transform(curveFromKey(pair.reference)), cmp=transform(curveFromKey(pair.comparison));
    if(!ref||!cmp) return;
    const series=selectedComparisonSeries(ref,cmp,comparisonMetric,smMap,document.getElementById("deltaAxis").value,document.getElementById("deltaAbscissa").value);
    if(!series.xs.length) return;
    yLabel=series.yLabel; xAxisTitle=series.xAxisTitle;
    allCurves.push(ref,cmp);
    const color=CONFIG_COLORS[index%CONFIG_COLORS.length];
    traces.push({{
      x:series.xs,y:series.ys,mode:"lines+markers",name:`${{cmp.curve_label}} − ${{ref.curve_label}}`,
      line:{{color,dash:POLAR_DASHES[index%POLAR_DASHES.length]}},marker:{{color}},
      text:series.xs.map((x,i)=>`${{xAxisTitle}}: ${{fmt(x,6)}}<br>${{yLabel}}: ${{fmt(series.ys[i],6)}}<br>${{series.comparisonNote}}<extra></extra>`),hovertemplate:"%{{text}}"
    }});
  }});
  if(!traces.length){{Plotly.newPlot("comparisonPlot",[],{{title:"No overlapping interpolation range for the selected comparisons"}});return;}}
  Plotly.newPlot("comparisonPlot",traces,{{
    title:`${{yLabel}} vs ${{xAxisTitle.replace(" [deg]","")}} · ${{referenceTitle()}} · ${{conditionText(allCurves)}}`,
    xaxis:{{title:xAxisTitle}},yaxis:{{title:yLabel,zeroline:true}},margin:{{l:70,r:30,t:55,b:50}}
  }},{{responsive:true}});
}}

function drawSelectedCaseSummary() {{
  const r=selectedCaseRow(), box=document.getElementById("selectedCaseSummary"); if(!r){{box.innerHTML="";return;}}
  const summary=summaryFor(r.case_label,r.polar), mesh=summary?.mesh||{{}}, defl=summary?.deflections||{{}};
  box.innerHTML=`
    <div class="status-box"><div class="label">Selected case</div><div class="value">${{r.case_label}} · ${{r.polar}} · CASE ${{r.case}}</div><div>α=${{fmt(r.alpha,2)}}°, β=${{fmt(r.beta,2)}}°</div></div>
    <div class="status-box"><div class="label">Numerical</div><div class="value ${{assessmentClass(r.numerical_status)}}">${{r.numerical_status}}</div><div>${{(r.numerical_reasons||[]).join(" · ")}}</div></div>
    <div class="status-box"><div class="label">Aerodynamic</div><div class="value ${{assessmentClass(r.aerodynamic_status)}}">${{r.aerodynamic_status}}</div><div>${{(r.aerodynamic_reasons||[]).join(" · ")}}</div></div>
    <div class="status-box"><div class="label">Neighbor</div><div class="value ${{assessmentClass(r.neighbor_status)}}">${{r.neighbor_status}}</div><div>${{(r.neighbor_reasons||[]).join(" · ")}}</div></div>
    <div class="status-box"><div class="label">Mesh</div><div class="value ${{assessmentClass(mesh.status)}}">${{mesh.status||"UNKNOWN"}}</div><div>Min OQ=${{fmt(mesh.min_orthogonal_quality,5)}} · ARmax=${{fmt(mesh.max_aspect_ratio,3)}}</div></div>
    <div class="status-box"><div class="label">Control surfaces</div><div>FLP1=${{defl.FLP1??"-"}} · AIL1=${{defl.AIL1??"-"}} · ELV1=${{defl.ELV1??"-"}} · RUD1=${{defl.RUD1??"-"}}</div></div>`;
  const stats=[
    ["CLZB",r.clzb_mean,r.clzb_min,r.clzb_max,r.clzb_std,r.clzb_drift_per_100],
    ["CDXB",r.cdxb_mean,r.cdxb_min,r.cdxb_max,r.cdxb_std,r.cdxb_drift_per_100],
    ["CMYB",r.cmyb_mean,r.cmyb_min,r.cmyb_max,r.cmyb_std,r.cmyb_drift_per_100],
    ["cp-max",r.cpmax_mean,r.cpmax_min,r.cpmax_max,r.cpmax_std,r.cpmax_drift_per_100]
  ];
  document.querySelector("#finalWindowStatsTable tbody").innerHTML=stats.map(s=>`<tr><td>${{s[0]}}</td><td>${{fmt(s[1],6)}}</td><td>${{fmt(s[2],6)}}</td><td>${{fmt(s[3],6)}}</td><td>${{fmt(s[4],6)}}</td><td>${{fmt(s[5],6)}}</td></tr>`).join("");
}}

const baseDrawSelectedHistoryV19=drawSelectedHistory;
drawSelectedHistory=function(){{baseDrawSelectedHistoryV19();drawSelectedCaseSummary();}};

function setDensity(value) {{
  document.body.classList.remove("density-compact","density-comfortable","density-presentation"); document.body.classList.add(`density-${{value}}`);
  setTimeout(()=>document.querySelectorAll(".plot").forEach(d=>{{try{{Plotly.Plots.resize(d)}}catch(e){{}}}}),20); saveLastState();
}}

function populateExportPlots() {{
  const ids=[...(document.getElementById("coefficientAbscissa").value==="CL"?[]:["adfCoeffPlot"]),"coefficientCDPlot",...(document.getElementById("coefficientAbscissa").value==="CY"?[]:["coefficientCYPlot"]),"coefficientCMPlot","coefficientCRPlot","coefficientCNPlot","ldPlot","smClPlot","smSweepPlot","comparisonPlot","residualEquationsPlot","cpmaxPlot","selectedResidualHistory","selectedAeroHistory","selectedCpHistory","outlierScorePlot"];
  const previous=document.getElementById("exportPlotSelect").value;
  document.getElementById("exportPlotSelect").innerHTML=ids.map(id=>`<option value="${{id}}">${{id}}</option>`).join("");
  if(ids.includes(previous)) document.getElementById("exportPlotSelect").value=previous;
}}
function exportSelectedPlot(format) {{
  const id=document.getElementById("exportPlotSelect").value||selectedPlotId; const div=document.getElementById(id); if(!div||!div.data)return alert("Plot is not available in the current view.");
  Plotly.downloadImage(div,{{format,filename:`JAMAL_${{id}}`,height:900,width:1400,scale:1.2}});
}}
function csvEscape(v){{const s=Array.isArray(v)?v.join(" | "):String(v??"");return /[",\\n]/.test(s)?`"${{s.replace(/"/g,'""')}}"`:s;}}
function exportFilteredCsv(){{const rows=getFilteredRows();if(!rows.length)return;const cols=Object.keys(rows[0]);const csv=[cols.join(","),...rows.map(r=>cols.map(c=>csvEscape(r[c])).join(","))].join("\\n");const a=document.createElement("a");a.href=URL.createObjectURL(new Blob([csv],{{type:"text/csv"}}));a.download="JAMAL_filtered_cases.csv";a.click();URL.revokeObjectURL(a.href);}}
async function copySelectedCaseSummary(){{const r=selectedCaseRow();if(!r)return;const txt=`${{r.case_label}} | ${{r.polar}} | CASE ${{r.case}} | α=${{r.alpha}} β=${{r.beta}}\nNumerical: ${{r.numerical_status}} — ${{(r.numerical_reasons||[]).join("; ")}}\nAerodynamic: ${{r.aerodynamic_status}} — ${{(r.aerodynamic_reasons||[]).join("; ")}}\nNeighbor: ${{r.neighbor_status}} — ${{(r.neighbor_reasons||[]).join("; ")}}\nOverall: ${{r.status}}`;try{{await navigator.clipboard.writeText(txt);alert("Selected-case summary copied.");}}catch(e){{prompt("Copy summary:",txt);}}}}

function collectState(){{return{{section:currentSectionName,label:document.getElementById("labelFilter").value,polar:document.getElementById("polarFilter").value,status:document.getElementById("statusFilter").value,axis:document.getElementById("stabilityAxisSelect").value,coeff:document.getElementById("stabilityCoeffSelect").value,momentRef:document.getElementById("momentReferenceSelect").value,refMode:document.getElementById("refModeSelect").value,xPct:document.getElementById("xrefPercentInput").value,yOff:document.getElementById("yOffsetInput").value,zOff:document.getElementById("zOffsetInput").value,xAbs:document.getElementById("xAbsInput").value,yAbs:document.getElementById("yAbsInput").value,zAbs:document.getElementById("zAbsInput").value,residual:document.getElementById("residualEquationSelect").value,density:document.getElementById("densitySelect").value,curveStyleMode:document.getElementById("curveStyleMode").value,customColorMapping:document.getElementById("customColorMapping").value,customLineMapping:document.getElementById("customLineMapping").value,customMarkerMapping:document.getElementById("customMarkerMapping").value,focusMode,thresholds:currentThresholds(),caseKey:document.getElementById("caseHistoryFilter").value,comparisonPairs:getComparisonPairs(),comparisonMetric}};}}
function applyState(s){{if(!s)return;const set=(id,v)=>{{const e=document.getElementById(id);if(e&&v!==undefined&&v!==null)e.value=v;}};set("labelFilter",s.label);set("polarFilter",s.polar);set("statusFilter",s.status);set("stabilityAxisSelect",s.axis);set("stabilityCoeffSelect",s.coeff);set("momentReferenceSelect",s.momentRef);set("refModeSelect",s.refMode);set("xrefPercentInput",s.xPct);set("yOffsetInput",s.yOff);set("zOffsetInput",s.zOff);set("xAbsInput",s.xAbs);set("yAbsInput",s.yAbs);set("zAbsInput",s.zAbs);set("residualEquationSelect",s.residual);set("densitySelect",s.density);set("curveStyleMode",s.curveStyleMode||"auto");set("customColorMapping",s.customColorMapping||"every_curve");set("customLineMapping",s.customLineMapping||"solid");set("customMarkerMapping",s.customMarkerMapping||"none");set("caseHistoryFilter",s.caseKey);updateCurveStyleControls();if(s.thresholds){{const t=s.thresholds;set("thrResSuccessful",t.residualSuccessful);set("thrResAcceptable",t.residualAcceptable);set("thrResBad",t.residualBad);set("thrClRange",t.clRange);set("thrCdRange",t.cdRange);set("thrCmRange",t.cmRange);set("thrCpRange",t.cpRange);set("thrNeighborCl",t.neighborCl);set("thrNeighborCd",t.neighborCd);set("thrNeighborCm",t.neighborCm);}}focusMode=s.focusMode||"all";comparisonMetric=s.comparisonMetric||"CL";setupComparisonControls(s.comparisonPairs || (s.comparisonReference ? [{{reference:s.comparisonReference,comparison:s.comparisonCurve}}] : null));setDensity(s.density||"comfortable");showSection(s.section||"overview",null);refreshAll();}}
function storageGet(key,fallback){{try{{return JSON.parse(localStorage.getItem(key))??fallback;}}catch(e){{return fallback;}}}}
function saveLastState(){{try{{localStorage.setItem("JAMAL_v24_last_state",JSON.stringify(collectState()));}}catch(e){{}}}}
function refreshPresetSelect(){{const p=storageGet("JAMAL_v24_presets",{{}});document.getElementById("presetSelect").innerHTML=Object.keys(p).sort().map(n=>`<option value="${{n}}">${{n}}</option>`).join("");}}
function saveNamedPreset(){{const n=document.getElementById("presetNameInput").value.trim();if(!n)return alert("Enter a preset name.");const p=storageGet("JAMAL_v24_presets",{{}});p[n]=collectState();localStorage.setItem("JAMAL_v24_presets",JSON.stringify(p));refreshPresetSelect();document.getElementById("presetSelect").value=n;}}
function loadNamedPreset(){{const n=document.getElementById("presetSelect").value,p=storageGet("JAMAL_v24_presets",{{}});applyState(p[n]);}}
function deleteNamedPreset(){{const n=document.getElementById("presetSelect").value,p=storageGet("JAMAL_v24_presets",{{}});delete p[n];localStorage.setItem("JAMAL_v24_presets",JSON.stringify(p));refreshPresetSelect();}}
function applyBuiltInPreset(name){{
  const nav = section => document.querySelector(`#dashboardNav button[data-section="${{section}}"]`);
  if(name==="longitudinal"){{document.getElementById("stabilityAxisSelect").value="S";document.getElementById("stabilityCoeffSelect").value="CM";showSection("aero",nav("aero"));}}
  if(name==="lateral"){{document.getElementById("stabilityAxisSelect").value="S";document.getElementById("stabilityCoeffSelect").value="CN";showSection("aero",nav("aero"));}}
  if(name==="convergence")showSection("convergence",nav("convergence"));
  if(name==="drag")showSection("drag",nav("drag"));
  refreshAll();
}}

function drawIntegrity(){{document.querySelector("#integrityTable tbody").innerHTML=(integrityData||[]).map(i=>`<tr><td class="integrity-${{i.severity}}">${{i.severity}}</td><td>${{i.configuration}}</td><td>${{i.polar}}</td><td>${{i.check}}</td><td>${{i.details}}</td></tr>`).join("");}}
function drawProvenance(){{document.getElementById("provenanceSummary").innerHTML=`Generated UTC: <strong>${{provenanceData.generated_utc}}</strong> · Script: <strong>${{provenanceData.script_version}}</strong>`;const rows=[];(provenanceData.adf_files||[]).forEach(f=>rows.push(["ADF",f.configuration,f.polar,f.path,f.modified]));(provenanceData.run_files||[]).forEach(f=>{{rows.push(["infout",f.configuration,f.polar,f.infout,f.infout_modified]);rows.push(["FLUENT_LOG",f.configuration,f.polar,f.fluent_log,f.fluent_log_modified]);}});(provenanceData.drag_rise_files||[]).forEach(f=>rows.push(["Drag rise",f.configuration,f.file_name,f.path,f.modified]));document.querySelector("#provenanceTable tbody").innerHTML=rows.map(r=>`<tr><td>${{r[0]}}</td><td>${{r[1]}}</td><td>${{r[2]}}</td><td>${{r[3]??""}}</td><td>${{r[4]??""}}</td></tr>`).join("");}}


drawCoeffPlot = function(data, divId, yKey, title, yTitle) {{
  const groups=groupedByLabelPolar(data), traces=[];
  Object.keys(groups).forEach(key=>{{const g=groups[key];if(!g.length)return;const st=traceStyle(g[0].case_label,g[0].polar);const xKey=(Math.max(...g.map(r=>Number(r.alpha??0)))-Math.min(...g.map(r=>Number(r.alpha??0)))) >= (Math.max(...g.map(r=>Number(r.beta??0)))-Math.min(...g.map(r=>Number(r.beta??0)))) ? "alpha" : "beta";traces.push({{x:g.map(r=>r[xKey]),y:g.map(r=>r[yKey]),mode:st.mode,name:key.replace("POLAR-","P"),line:{{color:st.color,dash:st.dash}},marker:{{color:st.color,symbol:st.symbol,size:st.markerSize}},customdata:g.map(r=>({{case_label:r.case_label,polar:r.polar,alpha:r.alpha,beta:r.beta}})),hovertemplate:`${{xKey.toUpperCase()}}: %{{x:.3f}}<br>${{yTitle}}: %{{y:.6f}}<extra>%{{fullData.name}}</extra>`}});}});
  const sweep=data.length&&data.some(r=>Math.abs(Number(r.beta??0))>1e-10)?"ALPHA/BETA":"ALPHA";
  Plotly.newPlot(divId,traces,{{title,xaxis:{{title:`${{sweep}} [deg]`}},yaxis:{{title:yTitle}},margin:{{l:65,r:30,t:50,b:50}}}},{{responsive:true}});
}};

drawResidualSummaryPlot = function(data) {{
  const groups=groupedByLabelPolar(data),traces=[];Object.keys(groups).forEach(key=>{{const g=groups[key];if(!g.length)return;const st=traceStyle(g[0].case_label,g[0].polar);traces.push({{x:g.map(r=>r.alpha),y:g.map(r=>r.worst_residual_final),mode:st.mode,name:key.replace("POLAR-","P"),line:{{color:st.color,dash:st.dash}},marker:{{color:st.color,symbol:st.symbol,size:st.markerSize}},customdata:g.map(r=>({{case_label:r.case_label,polar:r.polar,alpha:r.alpha,beta:r.beta}})),text:g.map(r=>`${{r.worst_residual_eq}} · ${{r.status}}`),hovertemplate:"α: %{{x:.3f}}<br>Residual: %{{y:.3e}}<br>%{{text}}<extra>%{{fullData.name}}</extra>"}});}});
  const t=currentThresholds();Plotly.newPlot("residualSummaryPlot",traces,{{title:"Worst final residual",xaxis:{{title:"Alpha [deg]"}},yaxis:sciAxis("Residual"),shapes:[{{type:"line",xref:"paper",x0:0,x1:1,y0:t.residualSuccessful,y1:t.residualSuccessful,line:{{dash:"dot"}}}},{{type:"line",xref:"paper",x0:0,x1:1,y0:t.residualAcceptable,y1:t.residualAcceptable,line:{{dash:"dot"}}}},{{type:"line",xref:"paper",x0:0,x1:1,y0:t.residualBad,y1:t.residualBad,line:{{dash:"dot"}}}}],margin:{{l:70,r:40,t:50,b:50}}}},{{responsive:true}});
}};

drawResidualEquationsPlot = function(data) {{
  const selected=document.getElementById("residualEquationSelect").value||"ALL",groups=groupedByLabelPolar(data),traces=[];const equations=selected==="ALL"?residualCols:[selected];
  Object.keys(groups).forEach(key=>{{const g=groups[key];if(!g.length)return;const base=traceStyle(g[0].case_label,g[0].polar);equations.forEach((eq,i)=>{{const all=selected==="ALL";traces.push({{x:g.map(r=>r.alpha),y:g.map(r=>r[eq+"_final"]),mode:all?"lines+markers":base.mode,name:`${{key.replace("POLAR-","P")}} · ${{eq}}`,line:{{color:base.color,dash:all?POLAR_DASHES[i%POLAR_DASHES.length]:base.dash}},marker:{{color:base.color,symbol:all?MARKER_SYMBOLS[i%MARKER_SYMBOLS.length]:base.symbol,size:all?6:base.markerSize}},customdata:g.map(r=>({{case_label:r.case_label,polar:r.polar,alpha:r.alpha,beta:r.beta}})),hovertemplate:`α: %{{x:.3f}}<br>${{eq}}: %{{y:.3e}}<extra>%{{fullData.name}}</extra>`}});}});}});
  Plotly.newPlot("residualEquationsPlot",traces,{{title:selected==="ALL"?"All final residual equations":`${{selected}} final residual`,xaxis:{{title:"Alpha [deg]"}},yaxis:sciAxis("Residual"),margin:{{l:70,r:40,t:50,b:50}}}},{{responsive:true}});
}};

drawCpmaxPlot = function(data) {{const groups=groupedByLabelPolar(data),traces=[];Object.keys(groups).forEach(key=>{{const g=groups[key];if(!g.length)return;const st=traceStyle(g[0].case_label,g[0].polar);traces.push({{x:g.map(r=>r.alpha),y:g.map(r=>r.cpmax_final),mode:st.mode,name:key.replace("POLAR-","P"),line:{{color:st.color,dash:st.dash}},marker:{{color:st.color,symbol:st.symbol,size:st.markerSize}},customdata:g.map(r=>({{case_label:r.case_label,polar:r.polar,alpha:r.alpha,beta:r.beta}})),hovertemplate:"α: %{{x:.3f}}<br>cp-max: %{{y:.5f}}<extra>%{{fullData.name}}</extra>"}});}});Plotly.newPlot("cpmaxPlot",traces,{{title:"Final cp-max",xaxis:{{title:"Alpha [deg]"}},yaxis:{{title:"cp-max"}},margin:{{l:60,r:30,t:50,b:50}}}},{{responsive:true}});}};

const baseDrawDragRisePlotsV19=drawDragRisePlots;
drawDragRisePlots=function(){{baseDrawDragRisePlotsV19();setTimeout(()=>document.querySelectorAll('[id^="dragRisePlot_"]').forEach(d=>bindPlotInteractions(d.id)),10);}};

const baseRefreshAllV19=refreshAll;
refreshAll=function(){{
  const data=getFilteredRows();updateMomentRefInfo();updateCaseHistoryOptions();drawKpis(data);drawAssessmentKpis(data);drawStatusMap(data);
  drawAdfCoeffPlot();drawSmPlots();drawComparisonPlot();drawDragRisePlots();
  drawResidualEquationsPlot(data);drawCpmaxPlot(data);drawClassificationTable(data);drawOutlierDiagnostics(data);drawMeshQuality();drawDeflections();drawDiagnostics();drawTable(data);drawSelectedCaseSummary();drawIntegrity();drawProvenance();const fm={{all:"All curves",baseline:"Baseline only",selected:"Selected curve only",selected_baseline:"Selected + baseline",none:"All curves hidden"}};document.getElementById("focusModeLabel").textContent=fm[focusMode]||focusMode;saveLastState();
}};

function updateCaseHistoryOptions(){{const select=document.getElementById("caseHistoryFilter"),old=select.value;const rows=evaluatedRows().slice().sort((a,b)=>(a.case_label+a.polar+a.case).localeCompare(b.case_label+b.polar+b.case));select.innerHTML=rows.map(r=>`<option value="${{r.case_key}}">${{r.case_label}} | ${{r.polar}} | CASE ${{r.case}} | α=${{r.alpha}} | ${{r.status}}</option>`).join("");if(rows.some(r=>r.case_key===old))select.value=old;}}

const collectStateBeforeAbscissa=collectState;
collectState=function(){{return {{...collectStateBeforeAbscissa(),coefficientAbscissa:document.getElementById("coefficientAbscissa").value,
  deltaAxis:document.getElementById("deltaAxis").value,deltaAbscissa:document.getElementById("deltaAbscissa").value,dragRiseAxis:document.getElementById("dragRiseAxis").value}};}};
const applyStateBeforeAbscissa=applyState;
applyState=function(s){{if(s){{
  ["coefficientAbscissa","deltaAbscissa"].forEach(id=>document.getElementById(id).value=["ALPHA","BETA","CL","CY"].includes(s[id])?s[id]:"ALPHA");
  ["deltaAxis","dragRiseAxis"].forEach(id=>document.getElementById(id).value=["W","S","B"].includes(s[id])?s[id]:"S");
}}return applyStateBeforeAbscissa(s);}};

setupPlotFirstWorkspace();
setupComparisonControls();populateExportPlots();refreshPresetSelect();drawIntegrity();drawProvenance();
const lastState=storageGet("JAMAL_v24_last_state",null);
if(lastState){{setTimeout(()=>applyState(convRows.length?lastState:{{...lastState,section:'aero'}}),10);}}else{{setDensity("comfortable");if(!convRows.length)showSection('aero',null);}}

initializeMomentReferenceInputs(); setupFilters(); drawDeflections(); drawDiagnostics(); drawMeshQuality(); drawModuleVersions(); refreshAll(); drawSelectedHistory();
</script>
</body>
</html>'''
    return jamal_distributions.add_to_html(html, distribution_data)


def write_html(path, summaries, conv_rows, history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data=None):
    Path(path).write_text(make_html(summaries, conv_rows, history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data), encoding="utf-8")


# =============================================================================
# MAIN
# =============================================================================

def main():
    case_configs = normalize_cases(CASES)
    if not case_configs:
        raise RuntimeError("CASES is empty. Add at least one configuration.")
    if len(case_configs) > MAX_CONFIGURATIONS:
        raise RuntimeError(f"Too many configurations: {len(case_configs)}. Maximum allowed is {MAX_CONFIGURATIONS}.")

    all_conv_rows = []
    all_summaries = []
    all_history = {}
    output_dir = None

    print("\nReading ADF polar files...")
    adf_polars = load_all_adf_polars(case_configs)
    sm_df = compute_all_static_margin(adf_polars)
    adf_data = make_adf_plot_rows(adf_polars, sm_df)

    print("\nReading optional drag-rise files...")
    drag_rise_data = load_drag_rise_data(case_configs)

    print("\nReading Fluent convergence logs...")
    for case_cfg in case_configs:
        _, runs_dir, convergence_dir = get_paths_from_adf_dir(case_cfg.directory)
        if output_dir is None:
            output_dir = convergence_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        for polar_number in case_cfg.polars:
            polar_name = polar_name_from_number(polar_number)
            polar_path = runs_dir / polar_name
            print(f"\nProcessing {case_cfg.label} | {polar_name}")
            print(f"  Run directory: {polar_path}")
            try:
                summary, rows, history = process_optional_convergence(case_cfg.label, polar_name, polar_path)
            except Exception as err:
                print(f"  ERROR: {err}")
                continue
            all_summaries.append(summary)
            all_conv_rows.extend(rows)
            all_history.update(history)

    all_conv_rows = attach_adf_coefficients_to_convergence(all_conv_rows, adf_data)
    all_conv_rows = add_neighbor_consistency(all_conv_rows)
    provenance = build_provenance(case_configs, all_summaries, adf_data, drag_rise_data)
    integrity_checks = build_integrity_checks(case_configs, all_summaries, all_conv_rows, adf_data, drag_rise_data)

    if output_dir is None:
        raise RuntimeError("No output directory could be determined from CASES.")

    distribution_configs = []
    for cfg in case_configs:
        results_dir, runs_dir, _ = get_paths_from_adf_dir(cfg.directory)
        distribution_configs.append({'label': cfg.label, 'base_directory': str(results_dir.parent),
                                     'runs_directory': str(runs_dir), 'polars': cfg.polars})
    distribution_data = jamal_distributions.read_distributions(
        distribution_configs, all_summaries, parse_infout, output_dir / '.jamal_cache' / 'distributions')
    if distribution_data['issues']:
        integrity_checks = [row for row in integrity_checks if row['severity'] != 'PASS'] + distribution_data['issues']
    provenance['distribution_sources'] = distribution_data['sources']

    html_path = output_dir / "dashboard.html"
    json_path = output_dir / "dashboard.json"

    write_html(html_path, all_summaries, all_conv_rows, all_history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data)
    write_json(json_path, all_summaries, all_conv_rows, all_history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data)

    print("\nDone.")
    print(f"HTML: {html_path}")
    print(f"JSON: {json_path}")


if __name__ == "__main__":
    main()
