"""Additive DISTCLCP reader. Physics tested with synthetic fixtures; production validation pending."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

VERSION = '0.2 synthetic-verified'
STATION_TOLERANCE = 1e-6  # metres
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?'
STATION_PATTERN = re.compile(rf'^(section|cp_dist)_state([1-9]\d*)_station({NUMBER})$', re.I)
FORCE_PATTERN = re.compile(r'^total_force_state([1-9]\d*)$', re.I)


def finite_number(token):
    value = float(token.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(value):
        raise ValueError('Non-finite numeric value')
    return value


def parse_curve(path):
    """Preserve surface order; skip all numeric metadata before the axis header."""
    points, started, ended = [], False, False
    with Path(path).open(encoding='utf-8-sig') as source:
        for line_number, line in enumerate(source, 1):
            if not started:
                started = 'ABSCISSA' in line.upper() and 'ORDINATE' in line.upper()
                continue
            text = line.strip()
            if text.upper().startswith('*END'):
                ended = True
                break
            if not text or text.startswith('$'):
                continue
            columns = text.split()
            if len(columns) != 2:
                raise ValueError(f'{path.name}:{line_number}: expected two curve columns')
            points.append([finite_number(value) for value in columns])
    if not started or not ended or len(points) < 2:
        raise ValueError(f'{path.name}: missing curve header, *END, or curve points')
    return points


def parse_force(path):
    rows = []
    with Path(path).open(encoding='utf-8-sig') as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip():
                continue
            columns = line.split()
            if len(columns) != 4:
                raise ValueError(f'{path.name}:{line_number}: expected Y FX FY FZ')
            rows.append([finite_number(value) for value in columns])
    if not rows:
        raise ValueError(f'{path.name}: empty force table')
    return rows


def fingerprint(path):
    stat = path.stat()
    return {'path': str(path.resolve()), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def cached_parse(path, kind, cache_dir, force, counts):
    """Cache raw geometry and per-state inputs independently from aerodynamic caches."""
    before = fingerprint(path)
    key = hashlib.sha256(str(path.resolve()).encode('utf-8')).hexdigest()
    target = cache_dir / kind / (key + '.json')
    if not force:
        try:
            cached = json.loads(target.read_text(encoding='utf-8'))
            if (cached['version'] == VERSION and cached['fingerprint'] == before
                    and isinstance(cached['data'], list) and cached['data']
                    and all(isinstance(row, list) and len(row) == (4 if kind == 'forces' else 2)
                            and all(isinstance(v, (int, float)) and math.isfinite(v) for v in row)
                            for row in cached['data'])):
                counts['cached'] += 1
                return cached['data'], before
        except (OSError, ValueError, KeyError, TypeError):
            pass
    data = parse_force(path) if kind == 'forces' else parse_curve(path)
    if fingerprint(path) != before:
        raise ValueError(f'{path.name}: source changed during parsing; retry generation')
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=key, suffix='.tmp', dir=target.parent)
    temp = Path(name)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as output:
            json.dump({'version': VERSION, 'fingerprint': before, 'data': data}, output, allow_nan=False)
        temp.replace(target)
    finally:
        temp.unlink(missing_ok=True)
    counts['parsed'] += 1
    return data, before


def read_distributions(configurations, summaries, parse_infout, cache_dir, force=False, progress=None):
    result = {'version': VERSION, 'station_tolerance': STATION_TOLERANCE,
              'series': [], 'issues': [], 'sources': [], 'counts': {'parsed': 0, 'cached': 0},
              'validation': 'Synthetic fixture verified; production validation pending'}

    def issue(cfg, polar, component, state, message):
        result['issues'].append({'severity': 'WARNING', 'configuration': cfg['label'],
                                 'polar': polar, 'check': 'Distributions',
                                 'details': f'{component} Â· state {state}: {message}'})

    for cfg in configurations:
        for number in cfg['polars']:
            polar = f'POLAR-{number:03d}'
            root = Path(cfg['base_directory']) / '03-RESULTS' / 'DISTCLCP' / polar
            if not root.is_dir():
                continue  # Optional module; no DISTCLCP is not an error.
            info = parse_infout(Path(cfg['runs_directory']) / polar / 'infout')
            meta, cases = info['meta'], info['cases']
            components = sorted((p for p in root.iterdir() if p.is_dir() and not p.is_symlink()), key=lambda p: p.name)
            for component in components:
                if progress:
                    progress(f"{cfg['label']} Â· {polar} Â· {component.name}")
                geometry, forces, cps = [], {}, {}
                outlines = {}
                span_reference = meta.get('bref')
                reference_source = 'infout BREF'
                reference_path = component / 'component_reference.json'
                if reference_path.is_file():
                    try:
                        reference = json.loads(reference_path.read_text(encoding='utf-8'))
                        span_reference = float(reference['span_reference_m'])
                        if not math.isfinite(span_reference) or span_reference <= 0:
                            raise ValueError('span_reference_m must be positive and finite')
                        reference_source = 'component_reference.json'
                        result['sources'].append({**fingerprint(reference_path), 'configuration': cfg['label'],
                                                  'polar': polar, 'component': component.name, 'kind': 'reference', 'state': None})
                    except (OSError, ValueError, TypeError, KeyError) as exc:
                        span_reference = None
                        issue(cfg, polar, component.name, 'all', f'invalid component span reference: {exc}')
                for path in sorted(component.iterdir()):
                    if not path.is_file() or path.is_symlink():
                        continue
                    station_match, force_match = STATION_PATTERN.fullmatch(path.name), FORCE_PATTERN.fullmatch(path.name)
                    if not station_match and not force_match:
                        continue
                    state = int((station_match or force_match).group(2 if station_match else 1))
                    kind = ('geometry' if station_match.group(1).lower() == 'section' else 'cp') if station_match else 'forces'
                    if kind == 'geometry' and state != 1:
                        continue
                    try:
                        rows, source = cached_parse(path, kind, cache_dir, force, result['counts'])
                        result['sources'].append({**source, 'configuration': cfg['label'], 'polar': polar,
                                                  'component': component.name, 'kind': kind, 'state': state})
                        if kind == 'geometry':
                            y = finite_number(station_match.group(3))
                            xs = [row[0] for row in rows]
                            chord = max(xs) - min(xs)
                            if chord <= 0:
                                raise ValueError('non-positive chord')
                            geometry.append({'y': y, 'xmin': min(xs), 'xmax': max(xs), 'chord': chord})
                            outlines[y] = rows
                        elif kind == 'cp':
                            cps.setdefault(state, []).append((finite_number(station_match.group(3)), rows))
                        else:
                            if state in forces:
                                raise ValueError('duplicate force state')
                            forces[state] = rows
                    except (OSError, ValueError, UnicodeError) as exc:
                        issue(cfg, polar, component.name, state, f'{path.name}: {exc}')

                geometry.sort(key=lambda g: g['y'])
                duplicate_geometry = {g['y'] for g in geometry
                                      if sum(abs(g['y'] - other['y']) <= STATION_TOLERANCE for other in geometry) > 1}
                if duplicate_geometry:
                    issue(cfg, polar, component.name, 1, 'duplicate geometry stations; ambiguous stations omitted')
                    geometry = [g for g in geometry if g['y'] not in duplicate_geometry]

                def match(y):
                    matches = [g for g in geometry if abs(g['y'] - y) <= STATION_TOLERANCE]
                    return matches[0] if len(matches) == 1 else None

                for unexpected in (set(forces) | set(cps)) - set(range(1, len(cases) + 1)):
                    issue(cfg, polar, component.name, unexpected, 'state has no matching infout case; omitted')
                for state, case in enumerate(cases, 1):
                    series = {'configuration': cfg['label'], 'polar': polar, 'component': component.name,
                              'state': state, 'case': case['case'], 'alpha': case['alpha'], 'beta': case['beta'],
                              'mach': case['mach'], 'reynolds': case['reynolds'], 'bref': meta.get('bref'),
                              'span_reference': span_reference, 'span_reference_source': reference_source,
                              'qdin': meta.get('qdin'), 'span': [], 'cp': []}
                    q = meta.get('qdin')
                    valid_q = q is not None and math.isfinite(q) and q > 0
                    # A single flow-condition qdin cannot be assigned across changing Mach/Re cases.
                    same_flow = all(case.get(k) is not None and cases[0].get(k) is not None
                                    and math.isclose(case[k], cases[0][k], rel_tol=1e-6, abs_tol=1e-9)
                                    for k in ('mach', 'reynolds'))
                    angles_valid = all(case.get(k) is not None and math.isfinite(case[k]) for k in ('alpha', 'beta'))
                    valid_lift = valid_q and same_flow and angles_valid and abs(case['beta']) <= 1e-9
                    if not valid_lift:
                        issue(cfg, polar, component.name, state, 'sectional cl unavailable: requires positive qdin, constant flow and beta=0')
                    if state not in forces:
                        issue(cfg, polar, component.name, state, 'missing force state')
                    force_rows = forces.get(state, [])
                    force_stations = [row[0] for row in force_rows]
                    used = set()
                    for y, fx, fy, fz in force_rows:
                        g = match(y)
                        if not g:
                            issue(cfg, polar, component.name, state, f'force station {y:g} has no unique geometry match')
                            continue
                        if sum(abs(y-other) <= 2*STATION_TOLERANCE for other in force_stations) > 1:
                            issue(cfg, polar, component.name, state, f'duplicate force station {y:g}; omitted')
                            continue
                        used.add(g['y'])
                        alpha = math.radians(case['alpha']) if angles_valid else 0
                        lift = fx * math.sin(alpha) - fz * math.cos(alpha) if angles_valid else None
                        series['span'].append({**g, 'force_y': y, 'fx': fx, 'fy': fy, 'fz': fz,
                                               'lift': lift, 'cl': lift / (q * g['chord']) if valid_lift else None})
                    for g in geometry:
                        if g['y'] not in used:
                            issue(cfg, polar, component.name, state, f"geometry station {g['y']:g} has no valid force station")
                    series['span'].sort(key=lambda r: r['y'])
                    if state not in cps:
                        issue(cfg, polar, component.name, state, 'missing Cp state')
                    cp_stations = [y for y, _ in cps.get(state, [])]
                    cp_used = set()
                    for y, points in cps.get(state, []):
                        g = match(y)
                        if not g:
                            issue(cfg, polar, component.name, state, f'Cp station {y:g} has no unique geometry match')
                            continue
                        if sum(abs(y-other) <= 2*STATION_TOLERANCE for other in cp_stations) > 1:
                            issue(cfg, polar, component.name, state, f'duplicate Cp station {y:g}; omitted')
                            continue
                        if any(x < g['xmin']-STATION_TOLERANCE or x > g['xmax']+STATION_TOLERANCE for x, _ in points):
                            issue(cfg, polar, component.name, state, f'Cp station {y:g} extends beyond geometry chord')
                        # Keep raw X as well as candidate x/c. The UI exposes the LE convention.
                        cp_used.add(g['y'])
                        series['cp'].append({**g, 'x': [p[0] for p in points], 'values': [p[1] for p in points],
                                             'xc': [(p[0]-g['xmin'])/g['chord'] for p in points],
                                             'airfoil': {'x': [r[0] for r in outlines[g['y']]],
                                                         'ordinate': [r[1] for r in outlines[g['y']]]}})
                    series['cp'].sort(key=lambda r: r['y'])
                    for g in geometry:
                        if g['y'] not in cp_used:
                            issue(cfg, polar, component.name, state, f"geometry station {g['y']:g} has no valid Cp station")
                    result['series'].append(series)
    return result


def add_to_html(html, data):
    """Read UI assets at generation time and embed them for a portable single-file report."""
    root = Path(__file__).parent
    panel = (root / 'distributions.html').read_text(encoding='utf-8')
    script = (root / 'distributions.js').read_text(encoding='utf-8')
    payload = json.dumps(data or {'series': [], 'issues': [], 'sources': []}, allow_nan=False).replace('<', '\\u003c')
    nav = '<button data-section="distributions" onclick="showSection(\'distributions\', this)">Distributions</button>'
    html = html.replace('<button data-section="convergence"', nav + '\n  <button data-section="convergence"', 1)
    return html.replace('</body>', panel + '\n<script>\nconst distributionData = ' + payload + ';\n' + script + '\n</script>\n</body>')
