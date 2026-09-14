"""Paired-tail fixture, component span references, and station panel grouping."""
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import jamal_dashboard_launcher_v25 as launcher
import jamal_distributions as dist

ROOT=Path(__file__).parent
FIXTURE=ROOT/'JAMAL_SYNTHETIC_CFD'


class VtailTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='jamal_vtail_')
        self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name)/'CFD'
        shutil.copytree(FIXTURE,self.base,ignore=shutil.ignore_patterns('DASHBOARD','__pycache__'))
        self.config=launcher._normalize_configurations({'configurations':[
            {'label':'Synthetic','base_directory':str(self.base),'polars':[1,2]}]})

    def load(self):
        return dist.read_distributions(self.config,[],launcher.ENGINE.parse_infout,Path(self.temp.name)/'cache')

    def test_paired_tail_dimensions_loads_and_pressure_consistency(self):
        data=self.load()
        self.assertEqual(data['issues'],[])
        self.assertEqual(len(data['series']),20)
        self.assertEqual(data['counts'],{'parsed':236,'cached':0})
        tail=[s for s in data['series'] if s['component']=='VTAIL']
        expected=json.loads((self.base/'expected_vtail_values.json').read_text())
        for s in tail:
            self.assertEqual(s['bref'],10)
            self.assertEqual(s['span_reference'],5)
            self.assertEqual([s['span'][0]['y'],s['span'][-1]['y']],[-2.5,2.5])
            reference=expected['polars'][s['polar']]['states'][s['state']-1]
            for span,cp,target in zip(s['span'],s['cp'],reference['sections']):
                self.assertAlmostEqual(span['cl'],target['cl'],places=8)
                self.assertAlmostEqual(span['fy'],target['FY'],places=5)
                self.assertAlmostEqual(span['chord'],target['chord'],places=8)
                self.assertAlmostEqual(min(cp['values'][:41]),target['Cp_upper_min'],places=8)
                self.assertEqual(len(cp['values']),81)
                upper=list(reversed(cp['values'][:41]));lower=[cp['values'][40]]+cp['values'][41:]
                xs=list(reversed(cp['xc'][:41]))
                delta=[b-a for a,b in zip(upper,lower)]
                integral=sum((a+b)*(d-c)/2 for a,b,c,d in zip(delta,delta[1:],xs,xs[1:]))
                self.assertAlmostEqual(integral,-span['fz']/(s['qdin']*span['chord']),places=8)
                if span['y']:
                    self.assertAlmostEqual(abs(span['fy']/span['fz']),1,places=8)
            fy=[p['fy'] for p in s['span']];ys=[p['y'] for p in s['span']]
            self.assertAlmostEqual(sum((a+b)*(d-c)/2 for a,b,c,d in zip(fy,fy[1:],ys,ys[1:])),0,places=6)
        first=next(s for s in tail if s['polar']=='POLAR-001' and s['state']==1)
        self.assertAlmostEqual(first['span'][4]['cl'],.18,places=9)
        self.assertEqual(self.load()['counts'],{'parsed':0,'cached':236})

    def test_component_reference_change_and_invalid_reference(self):
        self.load()
        path=self.base/'03-RESULTS/DISTCLCP/POLAR-001/VTAIL/component_reference.json'
        metadata=json.loads(path.read_text());metadata['span_reference_m']=6
        path.write_text(json.dumps(metadata))
        data=self.load()
        self.assertEqual(data['counts']['parsed'],0)
        self.assertTrue(all(s['span_reference']==6 for s in data['series'] if s['polar']=='POLAR-001' and s['component']=='VTAIL'))
        metadata['span_reference_m']=0;path.write_text(json.dumps(metadata))
        data=self.load()
        self.assertTrue(any('invalid component span reference' in i['details'] for i in data['issues']))
        self.assertTrue(all(s['span_reference'] is None for s in data['series'] if s['polar']=='POLAR-001' and s['component']=='VTAIL'))

    def test_airfoil_coordinates_preserve_section_order_and_cached_values(self):
        data = self.load()
        for series in data['series']:
            for cp in series['cp']:
                path = (self.base / '03-RESULTS/DISTCLCP' / series['polar'] /
                        series['component'] / f"section_state1_station{cp['y']:.3f}")
                points = dist.parse_curve(path)
                self.assertEqual(cp['airfoil']['x'], [p[0] for p in points])
                self.assertEqual(cp['airfoil']['ordinate'], [p[1] for p in points])
                normalized = [(x-cp['xmin'])/cp['chord'] for x in cp['airfoil']['x']]
                self.assertAlmostEqual(min(normalized), 0)
                self.assertAlmostEqual(max(normalized), 1)
        cached = self.load()
        self.assertEqual(cached['counts'], {'parsed': 0, 'cached': 236})
        self.assertEqual(cached['series'], data['series'])

    def test_original_and_additional_fixture_checksums(self):
        for manifest in ('SHA256SUMS.txt','VTAIL_SHA256SUMS.txt'):
            for line in (self.base/manifest).read_text().splitlines():
                checksum,name=line.split('  ',1)
                self.assertEqual(hashlib.sha256((self.base/name).read_bytes()).hexdigest(),checksum,name)

    def test_station_percentage_and_overlay_grouping_in_actual_javascript(self):
        js=(ROOT/'distributions.js').read_text(encoding='utf-8').split('(() => {',1)[0]
        js += r'''
const assert=require('node:assert/strict');
const wing={component:'WING',bref:10,span_reference:10};
const tail={component:'VTAIL',bref:10,span_reference:5};
assert.equal(distributionStation(wing,{y:5}).label,'WING · 2Y/BREF = +100%');
assert.equal(distributionStation(tail,{y:2.5}).label,'VTAIL · 2Y/BREF = +50%');
assert.equal(distributionStation(tail,{y:-1.25},'yb').label,'VTAIL · Y/BREF = -12.5%');
assert.equal(distributionStation(tail,{y:2.5},'y').label,'VTAIL · Y = 2.500 m');
assert.equal(distributionStation(wing,{y:2.5}).key,
             distributionStation({...wing,bref:20},{y:5}).key);
assert.notEqual(distributionStation(wing,{y:2.5}).key,distributionStation(tail,{y:2.5}).key);
assert.equal(distributionStation(tail,{y:2.5}).key,distributionStation(tail,{y:2.5},'yb').key);
assert.equal(distributionSpanReference({...tail,span_reference:null}),10);
assert.ok(distributionStation({...tail,bref:null},{y:2.5}).label.includes('BREF unavailable'));
assert.equal(distributionCoordinate(tail,2.5,'eta'),.5);
assert.equal(distributionCoordinate(tail,2.5,'yb'),.25);
assert.equal(distributionCoordinate(tail,2.5,'y'),2.5);
assert.equal(distributionChordLoad({cl:.2,chord:1.5}),.2*1.5);
assert.equal(distributionChordLoad({cl:0,chord:1.5}),0);
assert.equal(distributionChordLoad({cl:null,chord:1.5}),null);
assert.equal(distributionChordLoad({cl:.2,chord:0}),null);
const sample={...tail,configuration:'Test',polar:'POLAR-001',state:1,alpha:-1,beta:0,
 span:[{y:2.5,chord:.6,cl:.2},{y:0,chord:1.2,cl:null}]};
const rows=distributionLoadsText([sample], 'eta').trim().split('\r\n');
assert.equal(rows.length,6);
assert.deepEqual(rows[4].split('\t'),['Test','POLAR-001','VTAIL','1','-1','0','2.5','10','0.5','0.6','0.2','0.12']);
assert.deepEqual(rows[5].split('\t').slice(-2),['NA','NA']);
assert.equal(distributionLoadsText([sample], 'y').split('\r\n')[4].split('\t')[8],'2.5');
'''
        path=Path(self.temp.name)/'station_checks.js';path.write_text(js,encoding='utf-8')
        node=shutil.which('node');self.assertIsNotNone(node)
        result=subprocess.run([node,str(path)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)


if __name__=='__main__':
    unittest.main()
