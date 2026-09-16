"""Synthetic v25 regression checks; never reads the configured production paths.

Run: python -m unittest discover -s jamal_dashboard -p test_baseline.py -v
Node.js is required for the embedded JavaScript checks.
"""
import contextlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import jamal_dashboard_launcher_v25 as launcher

engine = launcher.ENGINE


def make_fixture(base):
    run = base / '02-RUNS' / 'POLAR-001'
    adf = base / '03-RESULTS' / 'ADF'
    run.mkdir(parents=True)
    adf.mkdir(parents=True)
    (run / 'infout').write_text(
        'POLAR: POLAR-001\n[REFERENCE_DATA]\n'
        'SREF[m2]: 40 CREF[m]: 4.629 BREF[m]: 12\n'
        'XREF[m]: 11.128 YREF[m]: 0 ZREF[m]: 0\n'
        '[CASES]\nCASE MACH REYNOLDS ALPHA BETA NAME ITERS\n'
        '0001 0.2 2E7 0 0 a 22\n'
        '0002 0.2 2E7 2 0 b 22\n'
        '0003 0.2 2E7 4 0 c 22\n', encoding='utf-8')
    header = ' '.join(engine.HISTORY_COLUMNS)
    lines = []
    for alpha, cl in [(0, .2), (2, .4), (4, .6)]:
        for i in range(1, 23):
            if (i - 1) % 11 == 0:
                lines.append(header)
            lines.append(' '.join(map(str, [i] + [1e-5] * 6 +
                                          [.01, 1, 0, -.1 * cl, 0, cl, 0, .02, .1])))
    (run / 'FLUENT_LOG').write_text('\n'.join(lines), encoding='utf-8')
    table = 'MACH REYNOLDS ALPHA BETA CDS CLS CMS25\n'
    table += '\n'.join(f'.2 2E7 {a} 0 .02 {cl} {-.1*cl}'
                       for a, cl in [(0, .2), (2, .4), (4, .6)])
    (adf / 'POLAR-001.adf').write_text(table, encoding='utf-8')
    return {'configurations': [{'label': 'Baseline', 'base_directory': str(base),
                                'polars': [1]}]}


class BaselineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='jamal_baseline_')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.payload = make_fixture(self.base)
        self.run = self.base / '02-RUNS' / 'POLAR-001'

    def generate(self, payload=None):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            launcher._run_generation_job('baseline-test', payload or self.payload)
        job = launcher.job_snapshot('baseline-test')
        self.assertEqual(job['status'], 'complete', job.get('error'))
        return job['result']

    def test_shallow_discovery_without_reading_source_contents(self):
        nested = self.base / '03-RESULTS' / 'ADF' / 'ADF_COMP'
        nested.mkdir()
        (nested / 'POLAR-999.adf').write_text('ignored')
        with patch.object(Path, 'read_text', side_effect=AssertionError('Discovery read contents')):
            scan = launcher.scan_base_directory(str(self.base))
        self.assertEqual([p['number'] for p in scan['polars']], [1])

    def test_reference_units_and_authoritative_case_order(self):
        parsed = engine.parse_infout(self.run / 'infout')
        self.assertEqual(parsed['meta']['cref'], 4.629)
        self.assertEqual(parsed['meta']['sref'], 40)
        self.assertEqual([r['case'] for r in parsed['cases']], ['0001', '0002', '0003'])

    def test_fluent_repeated_headers_and_incomplete_final_case(self):
        parsed = engine.parse_fluent_log_all_rows(self.run / 'FLUENT_LOG')
        cases = engine.parse_infout(self.run / 'infout')['cases']
        self.assertEqual(parsed['diagnostics']['history_headers_found'], 6)
        blocks, _ = engine.split_history_rows_by_infout(parsed['all_rows'], cases)
        self.assertEqual([len(b) for b in blocks], [22, 22, 22])
        blocks, _ = engine.split_history_rows_by_infout(parsed['all_rows'][:-7], cases)
        self.assertEqual([len(b) for b in blocks], [22, 22, 15])

    def test_static_margin_known_linear_slope(self):
        path = self.base / '03-RESULTS' / 'ADF' / 'POLAR-001.adf'
        polar = engine.PolarData('Baseline', 1, path, engine.read_polar_file(path), 'ALPHA')
        sm = engine.compute_static_margin_dataframe(polar)
        for value in sm['STATIC_MARGIN_PERCENT']:
            self.assertAlmostEqual(value, 10)

    def test_drag_rise_uses_lowest_mach_and_negative_filename(self):
        path = self.base / 'drag_rise_clm0p25.dat'
        path.write_text('MACH REYNOLDS ALPHA BETA CDS CLS\n'
                        '.8 2E7 0 0 .04 -.25\n.2 2E7 0 0 .02 -.25\n'
                        '.5 2E7 0 0 .03 -.25\n')
        self.assertEqual(engine.parse_drag_rise_cls_from_filename(path)[1], -.25)
        for actual, expected in zip(engine.read_drag_rise_file(path)['DELTA_CDS'], [0, .01, .02]):
            self.assertAlmostEqual(actual, expected)

    def test_generation_cache_reuse_invalidation_and_force(self):
        self.assertEqual(self.generate()['parsed_polars'], 1)
        with patch.object(engine, 'read_polar_file', side_effect=AssertionError('Reparsed ADF')), \
             patch.object(engine, 'process_polar_convergence', side_effect=AssertionError('Reparsed log')):
            self.assertEqual(self.generate()['reused_polars'], 1)
        log = self.run / 'FLUENT_LOG'
        stat = log.stat()
        os.utime(log, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        changed = self.generate()
        self.assertEqual(changed['plan']['modified'], 1)
        forced = self.generate({**self.payload, 'force_full_rebuild': True})
        self.assertEqual(forced['plan']['rebuild'], 1)
        output = self.base / '03-RESULTS' / 'DASHBOARD'
        self.assertTrue((output / 'dashboard.html').is_file())
        self.assertEqual(len(json.loads((output / 'dashboard.json').read_text())['results']), 3)
        self.assertFalse(list(output.glob('*.csv')))

    def test_python_and_embedded_javascript_syntax(self):
        for source in Path(launcher.__file__).parent.glob('*.py'):
            compile(source.read_text(encoding='utf-8'), str(source), 'exec')
        node = shutil.which('node')
        self.assertIsNotNone(node, 'Install Node.js to check the generated JavaScript')
        for name, html in [('launcher', launcher.INDEX_HTML), ('dashboard', self.html())]:
            scripts = re.findall(r'<script\b[^>]*>(.*?)</script>', html, re.S | re.I)
            for index, script in enumerate(scripts):
                path = self.base / f'{name}_{index}.js'
                path.write_text(script, encoding='utf-8')
                result = subprocess.run([node, '--check', str(path)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)

    def html(self):
        return engine.make_html([], [], {}, {'curves': [], 'static_margin': []},
                                {'curves': []}, {}, [])

    def test_embedded_moment_transfer_and_delta_numeric_regressions(self):
        html = self.html()
        # Execute the actual generated functions, without copying their implementation.
        def section(start, end):
            return html[html.index('function ' + start):html.index('function ' + end)]
        js = "const assert = require('node:assert/strict');\n"
        js += section('bodyToStabilityVector', 'currentMomentReferenceMode')
        js += section('interpolationPoints', 'setComparisonMetric')
        js += section('comparisonSeries', 'drawComparisonPlot')
        js += r'''
const close = (a,b) => assert.ok(Math.abs(a-b)<1e-10, `${a} != ${b}`);
let inputs = {mode:'percent', xPercent:0, yOffset:0, zOffset:0};
function currentMomentInputs(){return inputs;}
const meta = {xref:11.128,yref:0,zref:0,cref:4.629,bref:12};
close(newReferencePoint(meta, inputs).xnew,9.97075);
const row = {ALPHA:0,BETA:0,CDB:.02,CYB:0,CLB:.8,CRB25:0,CMB25:-.05,CNB25:0};
const shifted = shiftedRow(row,meta);
close(shifted.CMB25,-.25); close(shifted.CMS25,-.25); close(shifted.CMW25,-.25);
close(shifted.CLB,.8); close(shifted.CDB,.02); close(row.CMB25,-.05);
inputs.xPercent=25; close(shiftedRow(row,meta).CMB25,-.05);
inputs={mode:'absolute', xAbs:11.128, yAbs:1, zAbs:2};
const lateral=shiftedRow({...row,CYB:.1},meta);
close(lateral.CRB25,.05); close(lateral.CMB25,-.05-.04/4.629);
close(lateral.CNB25,-.02/12);
const pts=[{CLS:0,CDS:.02},{CLS:1,CDS:.05},{CLS:1,CDS:.07}];
close(interpolateAt(pts,'CLS',r=>r.CDS,.5),.04);
assert.equal(interpolateAt(pts,'CLS',r=>r.CDS,2),null);
const ref={rows:[{ALPHA:0,CLS:.25,CDS:.01},{ALPHA:2,CLS:.5,CDS:.02}]};
const cmp={rows:[{ALPHA:0,CLS:0,CDS:.02},{ALPHA:2,CLS:1,CDS:.06}]};
const cd=comparisonSeries(ref,cmp,'CD',{});
assert.deepEqual(cd.xs,[.25,.5]); cd.ys.forEach(v=>close(v,.02));
const cl=comparisonSeries(ref,cmp,'CL',{});
assert.deepEqual(cl.xs,[0,2]); assert.deepEqual(cl.ys,[-.25,.5]);
'''
        path = self.base / 'numeric_regressions.js'
        path.write_text(js, encoding='utf-8')
        node = shutil.which('node')
        self.assertIsNotNone(node, 'Install Node.js for numerical JavaScript checks')
        result = subprocess.run([node, str(path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_coefficient_grid_axis_abscissa_and_ratio(self):
        html = self.html()
        js = "const assert=require('node:assert/strict');\n"
        js += html[html.index('function coefficientKey('):html.index('function coefficientDisplayName(')]
        js += 'function coefficientDisplayName(axis,base){return base+axis;}\n'
        js += html[html.index('function coefficientPlotSpecs('):html.index('function drawStandardAeroPlots(')]
        js += r'''
for(const axis of ['B','S','W']) {
  for(const abscissa of ['ALPHA','CL']) {
    const specs=coefficientPlotSpecs(axis,abscissa);
    assert.equal(specs.length,abscissa==='ALPHA'?7:6);
    assert.equal(specs.some(s=>s.base==='CL'),abscissa==='ALPHA');
    specs.forEach(s=>assert.equal(s.xKey,abscissa==='ALPHA'?'ALPHA':'CL'+axis));
    assert.equal(specs.find(s=>s.base==='CM').yKey,'CM'+axis+'25');
    assert.equal(specs.find(s=>s.base==='CR').yKey,'CR'+axis+'25');
    assert.equal(specs.find(s=>s.base==='CN').yKey,'CN'+axis+'25');
    const rows=[{ALPHA:2,['CL'+axis]:.6,['CD'+axis]:.03},
                {ALPHA:0,['CL'+axis]:.2,['CD'+axis]:0},
                {ALPHA:1,['CL'+axis]:.4,['CD'+axis]:null}];
    const pts=coefficientPlotPoints(rows,specs.find(s=>s.base==='LD'));
    assert.deepEqual(pts.map(p=>p.y),[null,null,20]);
    assert.deepEqual(pts.map(p=>p.x),abscissa==='ALPHA'?[0,1,2]:[.2,.4,.6]);
    assert.equal(pts[2].row,rows[0]);
  }
}
'''
        path = self.base / 'coefficient_grid_checks.js'
        path.write_text(js, encoding='utf-8')
        result = subprocess.run([shutil.which('node'), str(path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
