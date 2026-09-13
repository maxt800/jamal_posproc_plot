"""Numerical and fault-injection tests using the user-supplied synthetic dataset."""
import ast
import contextlib
import io
import json
import math
from pathlib import Path
import shutil
import tempfile
import threading
import unittest
from unittest.mock import patch

import jamal_dashboard_launcher_v25 as launcher
import jamal_distributions as dist

ROOT = Path(__file__).parent
FIXTURE = ROOT / 'JAMAL_SYNTHETIC_CFD'
engine = launcher.ENGINE


class DistributionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='jamal_distributions_')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name) / 'CFD'
        # Keep the original WING-only regression workload stable as optional fixtures grow.
        shutil.copytree(FIXTURE, self.base, ignore=shutil.ignore_patterns('DASHBOARD', '__pycache__', 'VTAIL'))
        self.payload = {'configurations': [{'label': 'Synthetic', 'base_directory': str(self.base),
                                            'polars': [1, 2], 'drag_rise_dir': 'BASELINE'}]}
        self.normalized = launcher._normalize_configurations(self.payload)
        self.cache = Path(self.temp.name) / 'cache'
        self.component = self.base / '03-RESULTS/DISTCLCP/POLAR-001/WING'

    def load(self, force=False):
        return dist.read_distributions(self.normalized, [], engine.parse_infout, self.cache, force)

    def test_expected_values_every_state_station_and_pressure_branch(self):
        data = self.load()
        self.assertEqual(data['issues'], [])
        self.assertEqual(len(data['series']), 10)
        self.assertEqual(data['counts'], {'parsed': 118, 'cached': 0})
        expected = json.loads((FIXTURE / 'expected_values.json').read_text())
        for series in data['series']:
            state = expected['polars'][series['polar']]['states'][series['state']-1]
            self.assertEqual(series['alpha'], state['alpha'])
            self.assertEqual(len(series['span']), 9)
            self.assertEqual(len(series['cp']), 9)
            for span, cp, reference in zip(series['span'], series['cp'], state['sections']):
                self.assertAlmostEqual(span['cl'], reference['cl'], places=8)
                self.assertAlmostEqual(span['chord'], reference['chord'], places=8)
                self.assertAlmostEqual(span['y'], reference['Y'], places=8)
                self.assertAlmostEqual(2*span['y']/series['bref'], reference['Y']/5, places=8)
                self.assertEqual(len(cp['xc']), 81)
                self.assertAlmostEqual(cp['xc'][0], 1)
                self.assertAlmostEqual(cp['xc'][40], 0)
                self.assertAlmostEqual(cp['xc'][-1], 1)
                self.assertAlmostEqual(min(cp['values'][:41]), reference['Cp_upper_min'], places=8)
            values = series['span']
            integral = sum((a['lift']+b['lift'])*(b['y']-a['y'])/2 for a,b in zip(values,values[1:]))
            self.assertAlmostEqual(integral / (series['qdin']*15), state['CLS'], places=8)

    def test_separate_cache_reuse_change_and_force(self):
        self.load()
        with patch.object(dist, 'parse_curve', side_effect=AssertionError('reparsed')), \
             patch.object(dist, 'parse_force', side_effect=AssertionError('reparsed')):
            self.assertEqual(self.load()['counts'], {'parsed': 0, 'cached': 118})
        path = self.component / 'cp_dist_state3_station0.000'
        path.write_text(path.read_text() + '\n')
        self.assertEqual(self.load()['counts'], {'parsed': 1, 'cached': 117})
        self.assertEqual(self.load(force=True)['counts'], {'parsed': 118, 'cached': 0})

    def test_multiple_configurations_share_raw_cache_without_merging_series(self):
        self.normalized.append({**self.normalized[0], 'label': 'Comparison'})
        data = self.load()
        self.assertEqual(data['counts'], {'parsed': 118, 'cached': 118})
        self.assertEqual(len(data['series']), 20)
        self.assertEqual({s['configuration'] for s in data['series']}, {'Synthetic', 'Comparison'})
        self.assertEqual(data['issues'], [])

    def test_missing_geometry_force_and_unmapped_state_are_reported(self):
        (self.component / 'section_state1_station0.000').unlink()
        (self.component / 'total_force_state2').unlink()
        shutil.copy2(self.component / 'total_force_state1', self.component / 'total_force_state9')
        data = self.load()
        reasons = '\n'.join(i['details'] for i in data['issues'])
        self.assertIn('no unique geometry match', reasons)
        self.assertIn('missing force state', reasons)
        self.assertIn('no matching infout case', reasons)
        self.assertEqual(len(data['series']), 10)
        state2 = next(s for s in data['series'] if s['polar']=='POLAR-001' and s['state']==2)
        self.assertEqual(state2['span'], [])

    def test_dynamic_component_names_and_invalid_chords(self):
        path = self.component / 'section_state1_station0.000'
        path.write_text('$ ABSCISSA ORDINATE\n2 0\n2 1\n*END\n')
        self.component.rename(self.component.with_name('TEST_COMPONENT'))
        data = self.load()
        self.assertIn('TEST_COMPONENT', {s['component'] for s in data['series']})
        self.assertTrue(any('non-positive chord' in i['details'] for i in data['issues']))

    def test_nonzero_beta_and_missing_q_do_not_fabricate_cl(self):
        path = self.base / '02-RUNS/POLAR-001/infout'
        path.write_text(path.read_text().replace('qdin[Pa]', 'qimp_unused[Pa]'))
        data = self.load()
        self.assertTrue(all(p['cl'] is None for s in data['series'] if s['polar']=='POLAR-001' for p in s['span']))
        original = engine.parse_infout
        def nonzero_beta(path):
            value = original(path)
            for case in value['cases']:
                case['beta'] = 3
            return value
        with patch.object(engine, 'parse_infout', side_effect=nonzero_beta):
            self.assertTrue(all(p['cl'] is None for s in self.load()['series'] for p in s['span']))

    def test_curve_metadata_and_malformed_numeric_rows(self):
        self.assertEqual(len(dist.parse_curve(self.component / 'section_state1_station0.000')), 81)
        path = Path(self.temp.name) / 'invalid_curve'
        for text in ('$ ABSCISSA ORDINATE\n0 1\n1 nan\n*END',
                     '$ ABSCISSA ORDINATE\n0 1\n1 2',
                     '$ ABSCISSA ORDINATE\n0 1 2\n*END'):
            path.write_text(text)
            with self.assertRaises(ValueError):
                dist.parse_curve(path)

    def test_full_generation_reuse_and_writer_failure_preserve_report(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            launcher._run_generation_job('synthetic-first', self.payload)
        self.assertEqual(launcher.job_snapshot('synthetic-first')['status'], 'complete')
        output = self.base / '03-RESULTS/DASHBOARD'
        report = json.loads((output/'dashboard.json').read_text())
        self.assertEqual([r['actual_iters'] for r in report['results']], [240]*10)
        self.assertEqual(len(report['distributions']['series']), 10)
        self.assertEqual(len(report['drag_rise']['curves']), 2)
        for row in report['adf']['static_margin']:
            self.assertAlmostEqual(row['STATIC_MARGIN_PERCENT'], 10, places=6)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            launcher._run_generation_job('synthetic-reuse', self.payload)
        result = launcher.job_snapshot('synthetic-reuse')['result']
        self.assertEqual(result['reused_polars'], 2)
        self.assertEqual(result['distributions'], {'parsed': 0, 'cached': 118})
        before = {name: (output/name).read_bytes() for name in ('dashboard.html', 'dashboard.json')}
        with patch.object(engine, 'write_html', side_effect=OSError('injected writer failure')), \
             contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            launcher._run_generation_job('synthetic-failure', self.payload)
        self.assertEqual(launcher.job_snapshot('synthetic-failure')['status'], 'error')
        for name, content in before.items():
            self.assertEqual((output/name).read_bytes(), content)


class SafeguardTests(unittest.TestCase):
    def test_exponent_reference_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'infout'
            path.write_text('XREF[m]: 1.1128E+01 YREF[m]: -2D-1 ZREF[m]: +3e-2\nqdin[Pa]: 2.7744e3')
            meta = engine.parse_infout(path)['meta']
            self.assertEqual([meta[k] for k in ('xref','yref','zref','qdin')], [11.128,-.2,.03,2774.4])

    def test_duplicate_labels_rejected(self):
        cfg = {'label':'Same','base_directory':str(FIXTURE),'polars':[1]}
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            launcher._normalize_configurations({'configurations':[cfg, {**cfg, 'label':' same '}]})

    def test_second_generation_and_clear_cache_blocked_until_worker_finishes(self):
        started, finish = threading.Event(), threading.Event()
        def worker(*args):
            started.set()
            finish.wait(5)
        with patch.object(launcher, '_run_generation_job', side_effect=worker):
            try:
                launcher.start_generation_job({})
                self.assertTrue(started.wait(2))
                with self.assertRaisesRegex(RuntimeError, 'already active'):
                    launcher.start_generation_job({})
                with self.assertRaisesRegex(RuntimeError, 'during generation'):
                    launcher.clear_cache({})
            finally:
                finish.set()
                self.assertTrue(launcher.GENERATION_LOCK.acquire(timeout=2))
                launcher.GENERATION_LOCK.release()

    def test_report_pair_rolls_back_when_second_replace_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); stage = root/'stage'; stage.mkdir()
            for name in ('dashboard.json','dashboard.html'):
                (root/name).write_text('old'); (stage/name).write_text('new')
            original = Path.replace
            def replacement(path, destination):
                if path == stage/'dashboard.html':
                    raise OSError('injected replacement failure')
                return original(path, destination)
            with patch.object(Path, 'replace', replacement), self.assertRaises(OSError):
                launcher._publish_report_pair(stage, root)
            self.assertEqual((root/'dashboard.json').read_text(), 'old')
            self.assertEqual((root/'dashboard.html').read_text(), 'old')

    def test_json_script_terminators_are_escaped(self):
        html = engine.make_html([{'case_label':'</script><b>unsafe</b>'}],[],{},
                                {'curves':[], 'static_margin':[]},{'curves':[]},{},[])
        self.assertNotIn('</script><b>unsafe</b>', html)

    def test_protected_modules_unchanged_from_saved_v25(self):
        name = 'jamal_polar_convergence_dashboard_v25.py'
        before = (ROOT/'baselines/v25'/name).read_text(encoding='utf-8')
        after = (ROOT/name).read_text(encoding='utf-8')
        def functions(source):
            return {f.name: ast.dump(f) for f in ast.parse(source).body if isinstance(f, ast.FunctionDef)}
        old, new = functions(before), functions(after)
        for function in ('parse_fluent_log_all_rows','split_history_rows_by_infout',
                         'compute_static_margin_dataframe','read_drag_rise_file'):
            self.assertEqual(old[function], new[function], function)
        for start, end in [('bodyToStabilityVector','currentMomentReferenceMode'),
                           ('derivativeNonuniform','currentStaticMarginRows'),
                           ('interpolationPoints','setComparisonMetric'),
                           ('comparisonSeries','drawComparisonPlot')]:
            self.assertEqual(before[before.index('function '+start):before.index('function '+end)],
                             after[after.index('function '+start):after.index('function '+end)], start)


if __name__ == '__main__':
    unittest.main()
