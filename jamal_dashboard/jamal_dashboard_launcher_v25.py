#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JAMAL Local Dashboard Launcher v25
==================================

This launcher provides a two-phase local browser interface for selecting JAMAL data sources.
Phase 1 performs fast, shallow discovery only in 03-RESULTS/ADF and
03-RESULTS/DRAG-RISE. Phase 2 processes only the selected POLARs from
02-RUNS and generates the validated v25 standalone HTML dashboard.

No Flask or Dash dependency is required. The server uses only Python's standard
library. The validated dashboard engine remains in the sibling file:

    jamal_polar_convergence_dashboard_v25.py

Run:

    python jamal_dashboard_launcher_v25.py

Then use the browser interface to add up to five configurations, scan folders,
select POLARs, choose optional drag-rise folders, and generate the standalone
report.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import pickle
import re
import socket
import shutil
import subprocess
import sys
import threading
import tempfile
import time
import traceback
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urlparse


APP_VERSION = "v25.6"
ENGINE_FILENAME = "jamal_polar_convergence_dashboard_v25.py"
MAX_CONFIGURATIONS = 5
DEFAULT_OUTPUT_NAME = "dashboard"

STATE: Dict[str, Any] = {
    "last_report": None,
    "last_log": "",
    "last_error": None,
    "jobs": {},
}
JOB_LOCK = threading.Lock()
GENERATION_LOCK = threading.Lock()
REPORT_LOCK = threading.Lock()
CACHE_SCHEMA_VERSION = "1.0"


def load_engine():
    engine_path = Path(__file__).resolve().with_name(ENGINE_FILENAME)
    if not engine_path.exists():
        raise FileNotFoundError(
            f"Dashboard engine not found: {engine_path}\n"
            f"Keep {ENGINE_FILENAME} beside this launcher."
        )

    spec = importlib.util.spec_from_file_location("jamal_dashboard_engine_v25", engine_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load dashboard engine: {engine_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE = load_engine()


def sanitize_output_name(value: str) -> str:
    value = (value or DEFAULT_OUTPUT_NAME).strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("._") or DEFAULT_OUTPUT_NAME


def natural_key(text: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def resolve_base_directory(raw_path: str) -> Path:
    """Resolve the JAMAL base folder containing 02-RUNS and 03-RESULTS.

    The preferred input is the base folder itself. For convenience, paths inside
    02-RUNS or 03-RESULTS are also accepted and resolved back to the base folder.
    No recursive filesystem search is performed.
    """
    if not raw_path or not str(raw_path).strip():
        raise ValueError("JAMAL base directory is empty.")

    path = Path(str(raw_path).strip().strip('"')).expanduser()
    if path.is_file():
        path = path.parent

    candidates: List[Path] = [path]

    # Direct children of the base folder.
    if path.name.upper() in {"02-RUNS", "03-RESULTS"}:
        candidates.insert(0, path.parent)

    # Common locations inside 03-RESULTS, including a selected drag-rise case.
    current = path
    for _ in range(4):
        if current.name.upper() == "03-RESULTS":
            candidates.insert(0, current.parent)
            break
        if current.parent == current:
            break
        current = current.parent

    seen = set()
    for candidate in candidates:
        try:
            key = str(candidate.resolve())
        except Exception:
            key = str(candidate)
        if key in seen:
            continue
        seen.add(key)

        runs_dir = candidate / "02-RUNS"
        results_dir = candidate / "03-RESULTS"
        if candidate.is_dir() and (results_dir / 'ADF').is_dir():
            return candidate.resolve()

    raise FileNotFoundError(
        "Could not locate a JAMAL base folder containing 03-RESULTS/ADF. "
        "Select the base folder or its ADF directory."
    )


def resolve_adf_directory(raw_path: str) -> Path:
    """Return 03-RESULTS/ADF from a JAMAL base-folder path."""
    base_dir = resolve_base_directory(raw_path)
    adf_dir = base_dir / "03-RESULTS" / "ADF"
    if not adf_dir.is_dir():
        raise FileNotFoundError(f"ADF directory not found: {adf_dir}")
    return adf_dir.resolve()


def shallow_polar_listing(adf_dir: Path) -> List[Dict[str, Any]]:
    """List only POLAR-XXX.adf files directly inside ADF.

    This deliberately ignores ADF_COMP and every other subdirectory. os.scandir
    performs a single shallow directory enumeration, which is substantially
    faster on mapped/network drives than per-POLAR validation.
    """
    polars: List[Dict[str, Any]] = []
    pattern = re.compile(r"POLAR-(\d+)\.adf$", flags=re.IGNORECASE)
    with os.scandir(adf_dir) as entries:
        for entry in entries:
            match = pattern.fullmatch(entry.name)
            if not match:
                continue
            try:
                is_file = entry.is_file(follow_symlinks=False)
            except OSError:
                is_file = True
            if not is_file:
                continue
            number = int(match.group(1))
            polars.append({
                "number": number,
                "name": f"POLAR-{number:03d}",
                "adf_path": str(adf_dir / entry.name),
            })
    polars.sort(key=lambda item: item["number"])
    return polars


def shallow_drag_rise_listing(drag_rise_root: Path) -> List[str]:
    """List only first-level directories in 03-RESULTS/DRAG-RISE."""
    if not drag_rise_root.is_dir():
        return []
    names: List[str] = []
    with os.scandir(drag_rise_root) as entries:
        for entry in entries:
            try:
                if entry.is_dir(follow_symlinks=False):
                    names.append(entry.name)
            except OSError:
                continue
    return sorted(names, key=natural_key)


def scan_base_directory(raw_path: str) -> Dict[str, Any]:
    """Phase 1: fast, non-recursive discovery from the JAMAL base folder."""
    base_dir = resolve_base_directory(raw_path)
    runs_dir = base_dir / "02-RUNS"
    results_dir = base_dir / "03-RESULTS"
    adf_dir = results_dir / "ADF"
    drag_rise_root = results_dir / "DRAG-RISE"
    dashboard_dir = results_dir / "DASHBOARD"

    if not adf_dir.is_dir():
        raise FileNotFoundError(f"ADF directory not found: {adf_dir}")

    polars = shallow_polar_listing(adf_dir)
    drag_rise_dirs = shallow_drag_rise_listing(drag_rise_root)

    return {
        "base_directory": str(base_dir),
        "adf_directory": str(adf_dir),
        "results_directory": str(results_dir),
        "runs_directory": str(runs_dir),
        "drag_rise_root": str(drag_rise_root),
        "dashboard_directory": str(dashboard_dir),
        "polars": polars,
        "drag_rise_directories": drag_rise_dirs,
        "discovery_mode": "shallow",
    }


def choose_directory(initial_directory: Optional[str] = None) -> str:
    """Open a native folder dialog when Tk is available."""
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        initial = None
        if initial_directory:
            candidate = Path(initial_directory).expanduser()
            if candidate.exists():
                initial = str(candidate)
        selected = filedialog.askdirectory(
            title="Select JAMAL base directory",
            initialdir=initial,
            mustexist=True,
        )
        root.destroy()
        return selected or ""
    except Exception as exc:
        raise RuntimeError(
            "Native folder browsing is unavailable. Enter or paste the JAMAL base-folder path manually. "
            f"Details: {exc}"
        ) from exc


def output_directory(configurations: List[Dict[str, Any]]) -> Path:
    first_base = resolve_base_directory(configurations[0]["base_directory"])
    return first_base / "03-RESULTS" / "DASHBOARD"


def output_report_path(configurations: List[Dict[str, Any]]) -> Path:
    return output_directory(configurations) / "dashboard.html"


def _set_job(job_id: str, **updates: Any) -> None:
    with JOB_LOCK:
        job = STATE["jobs"].setdefault(job_id, {})
        job.update(updates)
        job["updated_at"] = time.time()


def _job_log(job_id: str, message: str) -> None:
    with JOB_LOCK:
        job = STATE["jobs"].setdefault(job_id, {})
        current = job.get("log", "")
        job["log"] = current + ("\n" if current else "") + str(message)
        job["updated_at"] = time.time()


def _file_fingerprint(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        return {"path": str(path.resolve()), "missing": True}
    stat = path.stat()
    try:
        resolved = str(path.resolve())
    except Exception:
        resolved = str(path)
    return {"path": resolved, "size": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


def _cache_slug(*parts: Any) -> str:
    raw = "|".join(str(part) for part in parts)
    return hashlib.sha1(raw.encode("utf-8", errors="replace")).hexdigest()


def _load_pickle(path: Path) -> Optional[Dict[str, Any]]:
    try:
        with path.open("rb") as handle:
            value = pickle.load(handle)
        return value if isinstance(value, dict) else None
    except Exception:
        return None


def _save_pickle(path: Path, value: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temp.replace(path)


def _normalize_configurations(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    configurations = payload.get("configurations") or []
    if not configurations:
        raise ValueError("Add at least one configuration.")
    if len(configurations) > MAX_CONFIGURATIONS:
        raise ValueError(f"A maximum of {MAX_CONFIGURATIONS} configurations is supported.")

    normalized: List[Dict[str, Any]] = []
    labels = set()
    for index, cfg in enumerate(configurations, start=1):
        label = str(cfg.get("label") or f"Config_{index}").strip()
        if not label or '|' in label:
            raise ValueError("Configuration labels must be non-empty and cannot contain '|'.")
        if label.casefold() in labels:
            raise ValueError(f"Duplicate configuration label: {label}. Choose a unique label.")
        labels.add(label.casefold())
        base_dir = resolve_base_directory(str(cfg.get("base_directory") or cfg.get("directory") or ""))
        adf_dir = base_dir / "03-RESULTS" / "ADF"
        runs_dir = base_dir / "02-RUNS"
        polars = sorted({int(value) for value in (cfg.get("polars") or [])})
        folder = cfg.get("drag_rise_dir")
        if folder and (Path(str(folder)).name != str(folder) or str(folder) in {'.', '..'}):
            raise ValueError(f"{label}: select a first-level drag-rise folder name.")
        if not polars:
            raise ValueError(f"{label}: select at least one POLAR.")
        normalized.append({
            "label": label,
            "base_directory": str(base_dir),
            "adf_directory": str(adf_dir),
            "runs_directory": str(runs_dir),
            "polars": polars,
            "drag_rise_dir": cfg.get("drag_rise_dir") or None,
        })
    return normalized


def _polar_source_info(cfg: Dict[str, Any], polar_number: int) -> Dict[str, Any]:
    polar_name = f"POLAR-{polar_number:03d}"
    adf_path = Path(cfg["adf_directory"]) / f"{polar_name}.adf"
    run_dir = Path(cfg["runs_directory"]) / polar_name
    infout_path = run_dir / "infout"
    log_path = ENGINE.find_fluent_log(run_dir)
    missing = [str(path) for path in (adf_path,) if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"{cfg['label']} | {polar_name}: missing source files: {missing}")
    return {
        "polar_name": polar_name,
        "adf_path": adf_path,
        "run_dir": run_dir,
        "infout_path": infout_path,
        "log_path": log_path,
        "fingerprints": {
            "adf": _file_fingerprint(adf_path),
            "infout": _file_fingerprint(infout_path),
            "fluent_log": _file_fingerprint(log_path),
        },
    }


def _cache_valid(payload: Optional[Dict[str, Any]], fingerprints: Dict[str, Any]) -> bool:
    if not payload:
        return False
    return (
        payload.get("cache_schema") == CACHE_SCHEMA_VERSION
        and payload.get("engine_version") == getattr(ENGINE, "SCRIPT_VERSION", None)
        and payload.get("fingerprints") == fingerprints
    )


def _build_update_plan(normalized: List[Dict[str, Any]], cache_dir: Path, force: bool) -> List[Dict[str, Any]]:
    plan: List[Dict[str, Any]] = []
    polar_cache_dir = cache_dir / "polars"
    for cfg in normalized:
        for polar_number in cfg["polars"]:
            source = _polar_source_info(cfg, polar_number)
            cache_file = polar_cache_dir / f"{_cache_slug(cfg['base_directory'], cfg['label'], source['polar_name'])}.pkl"
            cached = None if force else _load_pickle(cache_file)
            if force:
                action = "rebuild"
            elif _cache_valid(cached, source["fingerprints"]):
                action = "cached"
            elif cached:
                action = "modified"
            else:
                action = "new"
            plan.append({"cfg": cfg, "polar_number": polar_number, "source": source, "cache_file": cache_file, "cached_payload": cached, "action": action})
    return plan


def _drag_rise_incremental(normalized: List[Dict[str, Any]], cache_dir: Path, force: bool, job_id: str) -> tuple[Dict[str, Any], Dict[str, int]]:
    curves: List[Dict[str, Any]] = []
    counts = {"new": 0, "modified": 0, "cached": 0}
    cache_root = cache_dir / "drag_rise"
    for cfg in normalized:
        folder_name = cfg.get("drag_rise_dir")
        if not folder_name:
            continue
        drag_dir = Path(cfg["base_directory"]) / "03-RESULTS" / "DRAG-RISE" / str(folder_name)
        if not drag_dir.is_dir():
            _job_log(job_id, f"WARNING: drag-rise directory not found: {drag_dir}")
            continue
        files = sorted([p for p in drag_dir.iterdir() if p.is_file() and re.match(r"drag[-_]rise[-_]cl.*\.dat$", p.name, re.I)], key=lambda p: natural_key(p.name))
        for path in files:
            fingerprint = _file_fingerprint(path)
            cache_file = cache_root / f"{_cache_slug(cfg['base_directory'], cfg['label'], folder_name, path.name)}.pkl"
            cached = None if force else _load_pickle(cache_file)
            if _cache_valid(cached, {"drag_rise": fingerprint}):
                curve = cached["curve"]
                counts["cached"] += 1
            else:
                action = "modified" if cached else "new"
                cls_label, cls_value = ENGINE.parse_drag_rise_cls_from_filename(path)
                df = ENGINE.read_drag_rise_file(path)
                cols = [c for c in ["POLAR", "MACH", "REYNOLDS", "ALPHA", "BETA", "CDB", "CDW", "CDS", "CLB", "CLW", "CLS", "DELTA_CDS"] if c in df.columns]
                curve = {
                    "case_label": cfg["label"], "drag_rise_dir": str(folder_name), "path": str(path),
                    "modified": ENGINE.datetime.fromtimestamp(path.stat().st_mtime, tz=ENGINE.timezone.utc).isoformat(),
                    "file_name": path.name, "cls_label": cls_label, "cls_value": cls_value,
                    "rows": df[cols].to_dict(orient="records"),
                }
                _save_pickle(cache_file, {
                    "cache_schema": CACHE_SCHEMA_VERSION, "engine_version": getattr(ENGINE, "SCRIPT_VERSION", None),
                    "fingerprints": {"drag_rise": fingerprint}, "curve": curve,
                })
                counts[action] += 1
            curves.append(curve)
    return {"curves": curves}, counts


def _run_generation_job(job_id: str, payload: Dict[str, Any]) -> None:
    started = time.time()
    try:
        _set_job(job_id, status="running", percent=2, phase="Validating configuration", message="Resolving JAMAL base folders and selected POLARs.")
        normalized = _normalize_configurations(payload)
        force = bool(payload.get("force_full_rebuild"))
        output_dir = output_directory(normalized)
        output_dir.mkdir(parents=True, exist_ok=True)
        cache_dir = output_dir / ".jamal_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)

        _set_job(job_id, percent=7, phase="Building update plan", message="Checking selected source timestamps against the persistent cache.")
        plan = _build_update_plan(normalized, cache_dir, force)
        plan_counts = {name: sum(1 for item in plan if item["action"] == name) for name in ("new", "modified", "cached", "rebuild")}
        _set_job(job_id, plan=plan_counts)
        _job_log(job_id, "Update plan: " + " · ".join(f"{key}={value}" for key, value in plan_counts.items()))

        adf_polars = []
        summaries = []
        conv_rows = []
        history: Dict[str, Any] = {}
        total = max(1, len(plan))
        parsed_count = 0
        reused_count = 0

        for index, item in enumerate(plan, start=1):
            cfg, source = item["cfg"], item["source"]
            progress = 10 + int(62 * (index - 1) / total)
            _set_job(job_id, percent=progress, phase="Processing selected POLARs", message=f"{cfg['label']} · {source['polar_name']} · {item['action']}")
            if item["action"] == "cached":
                cached = item["cached_payload"]
                df = cached["adf_dataframe"]
                sweep_var = cached["sweep_var"]
                summary, rows, polar_history = cached["summary"], cached["rows"], cached["history"]
                reused_count += 1
            else:
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
                    df = ENGINE.read_polar_file(source["adf_path"])
                    sweep_var = ENGINE.detect_sweep_variable(df)
                    summary, rows, polar_history = ENGINE.process_optional_convergence(cfg["label"], source["polar_name"], source["run_dir"])
                text = buffer.getvalue().strip()
                if text:
                    _job_log(job_id, text)
                if _polar_source_info(cfg, item['polar_number'])['fingerprints'] != source['fingerprints']:
                    raise RuntimeError(f"{cfg['label']} · {source['polar_name']}: source changed during parsing; retry generation.")
                cache_payload = {
                    "cache_schema": CACHE_SCHEMA_VERSION,
                    "engine_version": getattr(ENGINE, "SCRIPT_VERSION", None),
                    "fingerprints": source["fingerprints"],
                    "adf_dataframe": df,
                    "sweep_var": sweep_var,
                    "summary": summary,
                    "rows": rows,
                    "history": polar_history,
                }
                _save_pickle(item["cache_file"], cache_payload)
                parsed_count += 1

            adf_polars.append(ENGINE.PolarData(case_label=cfg["label"], polar_number=item["polar_number"], path=source["adf_path"], data=df, sweep_var=sweep_var))
            summaries.append(summary)
            conv_rows.extend(rows)
            history.update(polar_history)

        _set_job(job_id, percent=73, phase="Processing drag rise", message="Loading only new or modified drag-rise files.")
        drag_rise_data, drag_counts = _drag_rise_incremental(normalized, cache_dir, force, job_id)

        _set_job(job_id, percent=78, phase="Processing distributions", message="Reading selected components and changed distribution files.")
        distribution_data = ENGINE.jamal_distributions.read_distributions(
            normalized, summaries, ENGINE.parse_infout, cache_dir / 'distributions', force,
            lambda message: _set_job(job_id, message=message))

        _set_job(job_id, percent=82, phase="Building derived results", message="Computing static margin, classifications, outliers, provenance, and integrity checks.")
        sm_df = ENGINE.compute_all_static_margin(adf_polars)
        adf_data = ENGINE.make_adf_plot_rows(adf_polars, sm_df)
        conv_rows = ENGINE.attach_adf_coefficients_to_convergence(conv_rows, adf_data)
        conv_rows = ENGINE.add_neighbor_consistency(conv_rows)
        case_configs = [ENGINE.CaseConfig(label=cfg["label"], directory=Path(cfg["adf_directory"]), polars=cfg["polars"], drag_rise_dir=cfg.get("drag_rise_dir")) for cfg in normalized]
        provenance = ENGINE.build_provenance(case_configs, summaries, adf_data, drag_rise_data)
        integrity_checks = ENGINE.build_integrity_checks(case_configs, summaries, conv_rows, adf_data, drag_rise_data)
        if distribution_data['issues']:
            integrity_checks = [row for row in integrity_checks if row['severity'] != 'PASS'] + distribution_data['issues']
        provenance['distribution_sources'] = distribution_data['sources']
        provenance['generation_id'] = job_id

        json_path = output_dir / "dashboard.json"
        html_path = output_dir / "dashboard.html"
        _set_job(job_id, percent=91, phase="Writing dashboard data", message=str(json_path))
        with tempfile.TemporaryDirectory(prefix='.jamal_stage_', dir=output_dir) as stage_name:
            stage = Path(stage_name)
            ENGINE.write_json(stage / json_path.name, summaries, conv_rows, history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data)
            _set_job(job_id, percent=96, phase="Writing standalone dashboard", message=str(html_path))
            ENGINE.write_html(stage / html_path.name, summaries, conv_rows, history, adf_data, drag_rise_data, provenance, integrity_checks, distribution_data)
            json.loads((stage / json_path.name).read_text(encoding='utf-8'))
            if not (stage / html_path.name).read_text(encoding='utf-8').rstrip().endswith('</html>'):
                raise ValueError('Generated HTML is incomplete.')
            _publish_report_pair(stage, output_dir)

        manifest = {
            "cache_schema": CACHE_SCHEMA_VERSION,
            "engine_version": getattr(ENGINE, "SCRIPT_VERSION", None),
            "generated_at": time.time(),
            "configurations": normalized,
            "plan": plan_counts,
            "parsed_polars": parsed_count,
            "reused_polars": reused_count,
            "drag_rise": drag_counts,
            "distributions": distribution_data['counts'],
            "dashboard_html": str(html_path),
            "dashboard_json": str(json_path),
        }
        manifest_temp = cache_dir / f'manifest.{job_id}.tmp'
        manifest_temp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        manifest_temp.replace(cache_dir / 'manifest.json')

        elapsed = time.time() - started
        result = {
            "report_path": str(html_path), "json_path": str(json_path), "report_url": "/report",
            "configurations": normalized, "parsed_polars": parsed_count, "reused_polars": reused_count,
            "drag_rise": drag_counts, "elapsed_seconds": elapsed, "plan": plan_counts,
            "distributions": distribution_data['counts'],
        }
        STATE["last_report"] = html_path
        STATE["last_log"] = STATE["jobs"][job_id].get("log", "")
        STATE["last_error"] = None
        _job_log(job_id, f"Dashboard updated · parsed={parsed_count} · reused={reused_count} · elapsed={elapsed:.1f}s")
        _set_job(job_id, status="complete", percent=100, phase="Complete", message="dashboard.html and dashboard.json are ready.", result=result)
    except Exception as exc:
        traceback.print_exc()
        _job_log(job_id, traceback.format_exc())
        STATE["last_error"] = str(exc)
        _set_job(job_id, status="error", phase="Failed", message=str(exc), error=str(exc))


def _publish_report_pair(stage: Path, output_dir: Path) -> None:
    """Publish completed files with rollback on ordinary write failures.

    The lock also protects HTTP readers. OS crashes and readers opening files
    directly cannot be made transactional across two paths by os.replace.
    """
    names = ('dashboard.json', 'dashboard.html')
    with REPORT_LOCK:
        existing = {name for name in names if (output_dir / name).exists()}
        for name in existing:
            shutil.copy2(output_dir / name, stage / (name + '.previous'))
        replaced = []
        try:
            for name in names:
                (stage / name).replace(output_dir / name)
                replaced.append(name)
        except Exception:
            for name in reversed(replaced):
                if name in existing:
                    (stage / (name + '.previous')).replace(output_dir / name)
                else:
                    (output_dir / name).unlink(missing_ok=True)
            raise


def start_generation_job(payload: Dict[str, Any]) -> str:
    if not GENERATION_LOCK.acquire(blocking=False):
        raise RuntimeError('Generation or cache clearing is already active. Wait for it to finish.')
    job_id = uuid.uuid4().hex
    with JOB_LOCK:
        STATE["jobs"][job_id] = {
            "job_id": job_id, "status": "queued", "percent": 0, "phase": "Queued",
            "message": "Waiting to start.", "log": "", "created_at": time.time(), "updated_at": time.time(),
        }
    def worker():
        try:
            _run_generation_job(job_id, payload)
        finally:
            GENERATION_LOCK.release()
    try:
        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
    except Exception:
        GENERATION_LOCK.release()
        _set_job(job_id, status='error', error='Could not start generation worker')
        raise
    return job_id


def job_snapshot(job_id: str) -> Dict[str, Any]:
    with JOB_LOCK:
        job = STATE["jobs"].get(job_id)
        if not job:
            raise KeyError(f"Unknown generation job: {job_id}")
        return dict(job)


def clear_cache(payload: Dict[str, Any]) -> Dict[str, Any]:
    if not GENERATION_LOCK.acquire(blocking=False):
        raise RuntimeError('Cannot clear cache during generation. Wait for it to finish.')
    try:
        normalized = _normalize_configurations(payload)
        output_dir = output_directory(normalized).resolve()
        cache_dir = output_dir / '.jamal_cache'
        if cache_dir.resolve().parent != output_dir or cache_dir.is_symlink():
            raise ValueError('Cache path must stay within the dashboard directory.')
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        return {"ok": True, "cache_directory": str(cache_dir)}
    finally:
        GENERATION_LOCK.release()

def open_in_file_manager(path: Path) -> None:
    target = path if path.is_dir() else path.parent
    if sys.platform.startswith("win"):
        os.startfile(str(target))  # type: ignore[attr-defined]
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(target)])
    else:
        subprocess.Popen(["xdg-open", str(target)])


INDEX_HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>JAMAL · Data source launcher</title>
<style>
:root {
  --bg:#f5f7fa; --card:#ffffff; --ink:#1f2937; --muted:#667085;
  --line:#d9dee8; --blue:#2563eb; --blue2:#1d4ed8; --green:#15803d;
  --amber:#b45309; --red:#b91c1c; --soft:#eef3fb;
}
* { box-sizing:border-box; }
body { margin:0; font-family:Arial,Helvetica,sans-serif; background:var(--bg); color:var(--ink); }
header { background:#fff; border-bottom:1px solid var(--line); padding:22px 28px; position:sticky; top:0; z-index:10; }
header h1 { margin:0 0 5px; font-size:24px; }
header p { margin:0; color:var(--muted); }
main { max-width:1380px; margin:0 auto; padding:22px; }
.toolbar { display:flex; gap:10px; flex-wrap:wrap; align-items:center; margin-bottom:18px; }
button,.button { border:1px solid #c9d2e2; background:white; color:#253047; border-radius:8px; padding:9px 13px; font-weight:600; cursor:pointer; }
button:hover { background:#f2f6fc; }
button.primary { background:var(--blue); border-color:var(--blue); color:white; }
button.primary:hover { background:var(--blue2); }
button.danger { color:var(--red); }
button:disabled { cursor:not-allowed; opacity:.55; }
input,select { border:1px solid #c9d2e2; border-radius:7px; padding:8px 10px; background:white; min-height:36px; }
input[type="text"] { width:100%; }
.config-card { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:17px; margin-bottom:16px; box-shadow:0 1px 3px rgba(16,24,40,.05); }
.config-head { display:flex; align-items:center; justify-content:space-between; gap:12px; margin-bottom:13px; }
.config-head h2 { margin:0; font-size:17px; }
.grid { display:grid; grid-template-columns:minmax(160px,.55fr) minmax(420px,2fr); gap:12px; align-items:end; }
.field label { display:block; color:#475467; font-size:12px; font-weight:700; margin-bottom:5px; }
.path-row { display:grid; grid-template-columns:1fr auto auto; gap:7px; }
.scan-meta { margin-top:12px; display:none; }
.meta-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:7px; font-size:12px; color:#475467; margin-bottom:12px; }
.meta-item { background:#f8fafc; padding:8px; border-radius:7px; overflow-wrap:anywhere; }
.polar-toolbar { display:flex; align-items:center; gap:8px; flex-wrap:wrap; margin-bottom:8px; }
.polar-list { display:grid; grid-template-columns:repeat(auto-fill,minmax(112px,1fr)); gap:7px; max-height:280px; overflow:auto; padding:9px; border:1px solid var(--line); border-radius:8px; background:#fbfcfe; }
.polar-option { display:flex; gap:7px; align-items:center; padding:6px 7px; border-radius:6px; background:white; border:1px solid #edf0f5; }
.polar-option small { display:block; color:var(--muted); margin-top:2px; line-height:1.25; }
.good { color:var(--green); } .warn { color:var(--amber); } .bad { color:var(--red); }
.options-row { display:grid; grid-template-columns:minmax(230px,1fr) minmax(230px,1fr); gap:12px; margin-top:12px; }
.panel { background:white; border:1px solid var(--line); border-radius:12px; padding:16px; margin-top:16px; }
.status { white-space:pre-wrap; font-family:Consolas,"Courier New",monospace; background:#101828; color:#e5e7eb; border-radius:8px; padding:13px; max-height:300px; overflow:auto; font-size:12px; }
.progress-wrap { margin:10px 0 12px; }
.progress-head { display:flex; justify-content:space-between; gap:12px; font-size:12px; font-weight:700; color:#475467; margin-bottom:6px; }
.progress-track { height:18px; border-radius:999px; overflow:hidden; background:#e7ecf3; border:1px solid #d4dce8; }
.progress-bar { height:100%; width:0%; background:linear-gradient(90deg,#2563eb,#4f83f1); color:white; font-size:11px; line-height:16px; text-align:center; transition:width .25s ease; min-width:0; }
.progress-plan { margin-top:7px; color:var(--muted); font-size:12px; }
.result { display:none; margin-top:14px; padding:14px; background:#ecfdf3; border:1px solid #abefc6; border-radius:9px; }
.result a { color:#175cd3; font-weight:700; }
.hidden { display:none !important; }
.file-input { display:none; }
@media (max-width:850px) { .grid,.options-row { grid-template-columns:1fr; } .path-row { grid-template-columns:1fr; } }
</style>
</head>
<body>
<header>
  <h1>JAMAL · CFD post-processing</h1>
  <p>Select a JAMAL base folder, discover ADF POLARs quickly, and process only the selected runs.</p>
</header>
<main>
  <div class="panel" style="margin-top:0;margin-bottom:16px">
    <h3 style="margin-top:0">Generation status</h3>
    <div class="progress-wrap">
      <div class="progress-head"><span id="progressPhase">Ready</span><span id="progressPercent">0%</span></div>
      <div class="progress-track"><div id="progressBar" class="progress-bar"></div></div>
      <div id="progressPlan" class="progress-plan">Incremental cache is checked only after Generate / Update dashboard.</div>
    </div>
    <div id="status" class="status">Ready. Phase 1 performs shallow discovery only.</div>
    <div id="result" class="result">
      <div id="resultPath"></div>
      <div style="margin-top:10px;display:flex;gap:8px;flex-wrap:wrap">
        <button class="primary" onclick="window.open('/report','_blank')">Open dashboard</button>
        <button onclick="openOutputFolder()">Open output folder</button>
      </div>
    </div>
  </div>

  <div class="toolbar">
    <button onclick="addConfiguration()">+ Add configuration</button>
    <button onclick="saveSetup()">Save setup</button>
    <button onclick="document.getElementById('setupFile').click()">Load setup</button>
    <input id="setupFile" class="file-input" type="file" accept="application/json,.json" onchange="loadSetupFile(event)">
    <div style="flex:1"></div>
    <label style="font-size:12px;font-weight:700;color:#475467"><input id="forceFullRebuild" type="checkbox"> Force full rebuild</label>
    <button onclick="clearDashboardCache()">Clear cache</button>
    <button id="generateButton" class="primary" onclick="generateDashboard()">Generate / Update dashboard</button>
  </div>

  <div id="configs"></div>

</main>
<script>
const MAX_CONFIGS = 5;
let nextId = 1;
const configs = new Map();

function esc(value) {
  return String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}
function setStatus(text) { document.getElementById('status').textContent = text; }
function setProgress(percent=0, phase='Ready', plan='') {
  const p=Math.max(0,Math.min(100,Number(percent)||0));
  document.getElementById('progressBar').style.width=`${p}%`;
  document.getElementById('progressBar').textContent=p>=12?`${Math.round(p)}%`:'';
  document.getElementById('progressPercent').textContent=`${Math.round(p)}%`;
  document.getElementById('progressPhase').textContent=phase||'Working';
  if(plan!==undefined) document.getElementById('progressPlan').textContent=plan||'';
}
function configCard(id, data={}) {
  const label = data.label || (id===1 ? 'Baseline' : `Config_${id}`);
  return `<section class="config-card" id="card-${id}">
    <div class="config-head"><h2>Configuration ${id}</h2><button class="danger" onclick="removeConfiguration(${id})">Remove</button></div>
    <div class="grid">
      <div class="field"><label>Label</label><input id="label-${id}" value="${esc(label)}"></div>
      <div class="field"><label>JAMAL base directory</label>
        <div class="path-row"><input id="path-${id}" value="${esc(data.base_directory||data.directory||'')}" placeholder="Z:\\...\\CFD">
        <button onclick="browseDirectory(${id})">Browse</button><button onclick="scanDirectory(${id})">Scan</button></div>
      </div>
    </div>
    <div class="scan-meta" id="scan-${id}">
      <div class="meta-grid" id="meta-${id}"></div>
      <div class="polar-toolbar"><strong>Available POLARs</strong><button onclick="setAllPolars(${id},true)">Select all</button><button onclick="setAllPolars(${id},false)">Clear</button><span id="count-${id}"></span></div>
      <div class="polar-list" id="polars-${id}"></div>
      <div class="options-row">
        <div class="field"><label>Optional drag-rise directory</label><select id="drag-${id}"><option value="">None</option></select></div>
        <div class="field"><label>Scan summary</label><div id="summary-${id}" class="meta-item"></div></div>
      </div>
    </div>
  </section>`;
}
function addConfiguration(data={}) {
  if (configs.size >= MAX_CONFIGS) { alert(`Maximum ${MAX_CONFIGS} configurations.`); return; }
  const id = nextId++;
  configs.set(id, {scan:null});
  document.getElementById('configs').insertAdjacentHTML('beforeend', configCard(id,data));
  if (data.base_directory || data.directory) scanDirectory(id, data.polars || null, data.drag_rise_dir || '');
  persistDraft();
}
function removeConfiguration(id) {
  if (configs.size <= 1) { alert('Keep at least one configuration.'); return; }
  configs.delete(id); document.getElementById(`card-${id}`)?.remove(); persistDraft();
}
async function postJson(url, payload={}) {
  const response = await fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});
  const data = await response.json();
  if (!response.ok || data.error) throw new Error(data.error || `HTTP ${response.status}`);
  return data;
}
async function browseDirectory(id) {
  try {
    setStatus('Opening folder selector...');
    const data = await postJson('/api/browse',{initial_directory:document.getElementById(`path-${id}`).value});
    if (data.path) { document.getElementById(`path-${id}`).value=data.path; await scanDirectory(id); }
    else setStatus('Folder selection cancelled.');
  } catch (error) { setStatus(`Browse error: ${error.message}`); }
}
async function scanDirectory(id, restorePolars=null, restoreDrag='') {
  const path = document.getElementById(`path-${id}`).value.trim();
  if (!path) { alert('Enter or browse to a directory.'); return; }
  try {
    setStatus(`Phase 1 · Fast shallow discovery for configuration ${id}...`);
    const data = await postJson('/api/scan',{base_directory:path});
    configs.get(id).scan=data;
    document.getElementById(`path-${id}`).value=data.base_directory;
    document.getElementById(`scan-${id}`).style.display='block';
    document.getElementById(`meta-${id}`).innerHTML=[
      ['Base',data.base_directory],['ADF',data.adf_directory],['02-RUNS',data.runs_directory],['DRAG-RISE',data.drag_rise_root],['Output',data.dashboard_directory]
    ].map(v=>`<div class="meta-item"><strong>${v[0]}:</strong> ${esc(v[1])}</div>`).join('');
    const wanted = restorePolars ? new Set(restorePolars.map(Number)) : new Set(data.polars.map(p=>p.number));
    document.getElementById(`polars-${id}`).innerHTML=data.polars.map(p=>{
      return `<label class="polar-option"><input type="checkbox" class="polar-check-${id}" value="${p.number}" ${wanted.has(p.number)?'checked':''} onchange="updateCount(${id})"><strong>${p.name}</strong></label>`;
    }).join('') || '<div>No POLAR-XXX.adf files found directly inside 03-RESULTS/ADF.</div>';
    const drag=document.getElementById(`drag-${id}`);
    drag.innerHTML='<option value="">None</option>'+data.drag_rise_directories.map(name=>`<option value="${esc(name)}">${esc(name)}</option>`).join('');
    if (restoreDrag && data.drag_rise_directories.includes(restoreDrag)) drag.value=restoreDrag;
    document.getElementById(`summary-${id}`).textContent=`Fast discovery: ${data.polars.length} root ADF POLAR files · ${data.drag_rise_directories.length} drag-rise folders · 02-RUNS checked only during generation`;
    updateCount(id); setStatus(`Phase 1 complete: ${data.polars.length} POLAR files found in ${data.adf_directory}. No recursive scan and no run validation performed.`); persistDraft();
  } catch (error) { setStatus(`Scan error: ${error.message}`); }
}
function setAllPolars(id, checked) { document.querySelectorAll(`.polar-check-${id}`).forEach(el=>el.checked=checked); updateCount(id); }
function selectedPolars(id) { return [...document.querySelectorAll(`.polar-check-${id}:checked`)].map(el=>Number(el.value)); }
function updateCount(id) { const n=selectedPolars(id).length; document.getElementById(`count-${id}`).textContent=`${n} selected`; persistDraft(); }
function serializeSetup() {
  return {version:'v25',configurations:[...configs.keys()].map(id=>({
    label:document.getElementById(`label-${id}`).value.trim(),base_directory:document.getElementById(`path-${id}`).value.trim(),polars:selectedPolars(id),drag_rise_dir:document.getElementById(`drag-${id}`)?.value||''
  }))};
}
function saveSetup() {
  const payload=serializeSetup(); const blob=new Blob([JSON.stringify(payload,null,2)],{type:'application/json'}); const a=document.createElement('a');
  a.href=URL.createObjectURL(blob); a.download='jamal_dashboard_setup.json'; a.click(); URL.revokeObjectURL(a.href);
}
async function loadSetupFile(event) {
  const file=event.target.files[0]; if(!file)return;
  try { const payload=JSON.parse(await file.text()); await applySetup(payload); } catch(error){ alert(`Invalid setup: ${error.message}`); }
  event.target.value='';
}
async function applySetup(payload) {
  document.getElementById('configs').innerHTML=''; configs.clear(); nextId=1;
  const list=(payload.configurations||[]).slice(0,MAX_CONFIGS); if(!list.length) list.push({label:'Baseline'});
  for (const item of list) addConfiguration(item);
  persistDraft();
}
function persistDraft() { try { localStorage.setItem('JAMAL_v25_source_setup',JSON.stringify(serializeSetup())); } catch (_) {} }
function restoreDraft() { try { const text=localStorage.getItem('JAMAL_v25_source_setup'); if(text) return JSON.parse(text); } catch (_) {} return null; }
let activeJobId = null;
let progressTimer = null;
let pollFailures = 0;
function planText(plan={}) {
  const parts=[];
  for(const key of ['new','modified','cached','rebuild']) if(Number(plan[key]||0)>0) parts.push(`${key}: ${plan[key]}`);
  return parts.length ? `Update plan · ${parts.join(' · ')}` : 'Preparing update plan...';
}
async function pollGenerationJob() {
  if(!activeJobId) return;
  try {
    const response=await fetch(`/api/progress?job_id=${encodeURIComponent(activeJobId)}`,{cache:'no-store'});
    if(response.status===404){
      activeJobId=null; progressTimer=null;
      document.getElementById('generateButton').disabled=false;
      setProgress(0,'Job unavailable','The launcher may have restarted. You can generate again.');
      setStatus('The generation job is no longer available on this launcher. Existing dashboard files are kept.');
      return;
    }
    const data=await response.json();
    if(!response.ok || data.error) throw new Error(data.error||`HTTP ${response.status}`);
    pollFailures=0;
    setProgress(data.percent||0,data.phase||data.status,planText(data.plan||{}));
    setStatus(`${data.phase||'Working'} · ${data.message||''}\n${data.log||''}`);
    if(data.status==='complete') {
      clearInterval(progressTimer); progressTimer=null; activeJobId=null;
      const r=data.result||{};
      setProgress(100,'Complete',`POLARs: parsed ${r.parsed_polars||0} · reused ${r.reused_polars||0} | Distribution files: parsed ${r.distributions?.parsed||0} · reused ${r.distributions?.cached||0} | ${(r.elapsed_seconds||0).toFixed(1)} s`);
      document.getElementById('resultPath').innerHTML=`Standalone dashboard: <strong>${esc(r.report_path||'')}</strong><br>Data snapshot: <strong>${esc(r.json_path||'')}</strong>`;
      document.getElementById('result').style.display='block';
      document.getElementById('generateButton').disabled=false; persistDraft(); return;
    }
    if(data.status==='error') {
      clearInterval(progressTimer); progressTimer=null; activeJobId=null;
      document.getElementById('generateButton').disabled=false;
      setProgress(data.percent||0,'Failed',data.error||data.message||'Generation failed.'); return;
    }
  } catch(error) {
    pollFailures++;
    setStatus(`Progress connection interrupted: ${error.message}. Retrying; the generation job may still be running.`);
  }
  if(activeJobId) progressTimer=setTimeout(pollGenerationJob,Math.min(10000,600*Math.pow(2,Math.min(pollFailures,4))));
}
async function generateDashboard() {
  const payload=serializeSetup();
  payload.force_full_rebuild=document.getElementById('forceFullRebuild').checked;
  for(const cfg of payload.configurations){ if(!cfg.base_directory){alert(`${cfg.label||'Configuration'}: base directory is empty.`);return;} if(!cfg.polars.length){alert(`${cfg.label||'Configuration'}: select at least one POLAR.`);return;} }
  const button=document.getElementById('generateButton'); button.disabled=true; document.getElementById('result').style.display='none';
  setProgress(1,'Starting','Submitting incremental update job...'); setStatus('Starting dashboard generation...');
  try {
    const data=await postJson('/api/generate',payload); activeJobId=data.job_id;
    await pollGenerationJob();
  } catch(error){ setStatus(`Generation error: ${error.message}`); setProgress(0,'Failed',error.message); button.disabled=false; }
}
async function clearDashboardCache() {
  const payload=serializeSetup();
  if(!payload.configurations.length || !payload.configurations[0].base_directory){ alert('Select the first configuration base folder before clearing cache.'); return; }
  if(!confirm('Clear the persistent JAMAL dashboard cache? Generated dashboard files will be kept.')) return;
  try { const data=await postJson('/api/clear-cache',payload); setStatus(`Cache cleared: ${data.cache_directory}`); setProgress(0,'Cache cleared','The next generation will parse selected files again.'); }
  catch(error){ setStatus(`Clear cache error: ${error.message}`); }
}

async function openOutputFolder(){ try{await postJson('/api/open-output');}catch(error){alert(error.message);} }
window.addEventListener('beforeunload',persistDraft);
const draft=restoreDraft(); if(draft) applySetup(draft); else addConfiguration({label:'Baseline'});
fetch('/api/state',{cache:'no-store'}).then(r=>r.json()).then(state=>{
  if(state.active_job_id && !activeJobId){activeJobId=state.active_job_id;document.getElementById('generateButton').disabled=true;pollGenerationJob();}
}).catch(()=>{});
</script>
</body>
</html>'''


class JamalRequestHandler(BaseHTTPRequestHandler):
    server_version = f"JAMALDashboard/{APP_VERSION}"

    def log_message(self, fmt: str, *args) -> None:
        print(f"[HTTP] {self.address_string()} - {fmt % args}")

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, payload: Dict[str, Any], status: int = 200) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self._send(status, body, "application/json; charset=utf-8")

    def _read_json(self) -> Dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length else b"{}"
        return json.loads(raw.decode("utf-8"))

    def do_GET(self) -> None:
        path = urlparse(self.path).path
        if path == "/":
            self._send(200, INDEX_HTML.encode("utf-8"), "text/html; charset=utf-8")
            return
        if path == "/report":
            report = STATE.get("last_report")
            if not report or not Path(report).exists():
                self._send(404, b"No generated report is available.", "text/plain; charset=utf-8")
                return
            with REPORT_LOCK:
                report_bytes = Path(report).read_bytes()
            self._send(200, report_bytes, "text/html; charset=utf-8")
            return
        if path == "/api/state":
            with JOB_LOCK:
                active_job_id = next((key for key, value in STATE['jobs'].items()
                                      if value.get('status') in {'queued', 'running'}), None)
            self._json({
                "active_job_id": active_job_id,
                "last_report": str(STATE["last_report"]) if STATE.get("last_report") else None,
                "last_log": STATE.get("last_log", ""),
                "last_error": STATE.get("last_error"),
            })
            return
        if path == "/api/progress":
            from urllib.parse import parse_qs
            query = parse_qs(urlparse(self.path).query)
            job_id = (query.get("job_id") or [""])[0]
            try:
                self._json(job_snapshot(job_id))
            except KeyError as exc:
                self._json({"error": str(exc)}, status=404)
            return
        if path == "/favicon.ico":
            self._send(204, b"", "image/x-icon")
            return
        self._send(404, b"Not found", "text/plain; charset=utf-8")

    def do_POST(self) -> None:
        path = urlparse(self.path).path
        try:
            payload = self._read_json()
            if path == "/api/browse":
                selected = choose_directory(payload.get("initial_directory"))
                self._json({"path": selected})
                return
            if path == "/api/scan":
                self._json(scan_base_directory(str(payload.get("base_directory") or payload.get("directory") or "")))
                return
            if path == "/api/generate":
                self._json({"job_id": start_generation_job(payload)})
                return
            if path == "/api/clear-cache":
                self._json(clear_cache(payload))
                return
            if path == "/api/open-output":
                report = STATE.get("last_report")
                if not report:
                    raise RuntimeError("No report has been generated yet.")
                open_in_file_manager(Path(report))
                self._json({"ok": True})
                return
            self._json({"error": "Unknown endpoint."}, status=404)
        except Exception as exc:
            traceback.print_exc()
            self._json({"error": str(exc)}, status=400)


def find_free_port(host: str, preferred: int) -> int:
    if preferred > 0:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            try:
                sock.bind((host, preferred))
                return preferred
            except OSError:
                pass
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((host, 0))
        return int(sock.getsockname()[1])


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the local JAMAL two-phase data-source dashboard.")
    parser.add_argument("--host", default="127.0.0.1", help="Local bind address.")
    parser.add_argument("--port", type=int, default=8765, help="Preferred local port.")
    parser.add_argument("--no-browser", action="store_true", help="Do not open the browser automatically.")
    args = parser.parse_args()

    port = find_free_port(args.host, args.port)
    url = f"http://{args.host}:{port}/"
    server = ThreadingHTTPServer((args.host, port), JamalRequestHandler)

    print(f"JAMAL local dashboard {APP_VERSION}")
    print(f"Engine: {ENGINE_FILENAME}")
    print(f"Open: {url}")
    print("Press Ctrl+C to stop the server.")

    if not args.no_browser:
        threading.Timer(0.7, lambda: webbrowser.open(url)).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping JAMAL dashboard server.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
