"""Validate this synthetic fixture using Python 3's standard library.

Run: python validate_test_data.py [path/to/JAMAL_SYNTHETIC_CFD]
This checks fixture integrity, not compatibility with an installed dashboard.
"""
from pathlib import Path
import json, math, re, sys

BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
checks=0
def check(ok, message):
    global checks
    checks+=1
    if not ok: raise AssertionError(message)
def near(a,b,label,tol=2e-7):
    check(math.isclose(a,b,rel_tol=tol,abs_tol=tol),f'{label}: {a} != {b}')
def trap(v,x):
    return sum((a+b)*(d-c)/2 for a,b,c,d in zip(v,v[1:],x,x[1:]))
def table(path):
    lines=path.read_text().splitlines()
    idx=next(i for i,l in enumerate(lines) if 'MACH' in l.split() and 'CDS' in l.split())
    keys=lines[idx].split()
    return [dict(zip(keys,map(float,l.split()))) for l in lines[idx+1:] if l.strip() and not l.startswith('#')]
def curve(path):
    lines=path.read_text().splitlines()
    idx=next(i for i,l in enumerate(lines) if 'ABSCISSA' in l and 'ORDINATE' in l)
    end=next(i for i in range(idx+1,len(lines)) if lines[i].startswith('*END'))
    pts=[tuple(map(float,l.split())) for l in lines[idx+1:end] if l.strip()]
    check(len(pts)==81,f'{path.name}: expected 81 points, excluding metadata')
    check(all(len(p)==2 and all(math.isfinite(x) for x in p) for p in pts),'finite two-column curves')
    return pts

expected=json.loads((BASE/'expected_values.json').read_text())
adf=BASE/'03-RESULTS/ADF'
check(sorted(p.stem for p in adf.glob('POLAR-*.adf'))==['POLAR-001','POLAR-002'],'shallow discovery')
check(len(list(adf.rglob('POLAR-*.adf')))==3,'nested decoy exists')
for polar,pe in expected['polars'].items():
    run=BASE/'02-RUNS'/polar
    info=(run/'infout').read_text()
    number=r'([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)'
    def metadata(label): return float(re.search(re.escape(label)+r':\s*'+number,info)[1])
    q=metadata('qdin[Pa]'); bref=metadata('BREF[m]'); sref=metadata('SREF[m2]')
    near(q,pe['qdin'],'qdin'); near(metadata('Reynolds'),pe['Reynolds'],'Reynolds')
    near(q,.5*metadata('rho[kg/m3]')*metadata('V[m/s]')**2,'dynamic pressure')
    near(metadata('Reynolds'),metadata('rho[kg/m3]')*metadata('V[m/s]')*metadata('CREF[m]')/metadata('mu[Pa.s]'),'Re from flow')
    for label in ['SREF','CREF','BREF','XREF','YREF','ZREF']:
        near(metadata(label+('[m2]' if label=='SREF' else '[m]')),expected['geometry'][label],label)
    cases=[l.split() for l in info.split('[CASES]')[1].splitlines() if re.match(r'^\d{4}\s',l)]
    rows=table(adf/(polar+'.adf'))
    check(len(cases)==len(rows)==5,'five cases and ADF rows')
    # The angle banner delimits states; repeated iter headers never do.
    # Simpler line-based scanner keeps solver chatter and repeated headers out.
    histories=[]; header_counts=[]; solves=[]
    for line in (run/'FLUENT_LOG').read_text().splitlines():
        if re.search(r';\*+\s+ALPHA:',line):
            histories.append([]);header_counts.append(0);solves.append(0)
        elif histories and line.strip().startswith('iter '): header_counts[-1]+=1
        elif histories and 'solve  iterate' in line: solves[-1]+=1
        elif histories and re.match(r'^\s*\d+\s+[+-]?\d+\.\d+[Ee]',line):
            t=line.split(); histories[-1].append([float(x) for x in t[:15]])
    check(len(histories)==5,'five log states')
    wing=BASE/'03-RESULTS/DISTCLCP'/polar/'WING'
    check(len(list(wing.glob('section_state1_station*')))==9,'nine geometry files')
    check(len(list(wing.glob('cp_dist_state*_station*')))==45,'45 Cp files')
    check(len(list(wing.glob('total_force_state*')))==5,'five force files')
    for case,r,se,h,hcount,solve in zip(cases,rows,pe['states'],histories,header_counts,solves):
        s=se['state']; alpha=float(case[3]);a=math.radians(alpha)
        near(int(case[0]),s,'state order');near(alpha,se['alpha'],'alpha');near(r['ALPHA'],alpha,'ADF angle')
        near(float(case[1]),r['MACH'],'case Mach');near(float(case[2]),r['REYNOLDS'],'case Re')
        near(float(case[4]),0,'beta');near(int(case[6]),len(h),'ITERS')
        check(len(h)==240 and hcount==12 and solve==2,'concatenate monitor subblocks')
        check([int(x[0]) for x in h]==list(range(se['first_iteration'],se['last_iteration']+1)),'continuous iteration sequence')
        near(h[-1][1],se['final_continuity'],'final residual');near(h[-1][8],1.05,'final cp-max')
        near(max(x[8] for x in h),3.6,'early cp spike')
        for idx,key in [(9,'CNB25'),(10,'CMB25'),(11,'CRB25'),(12,'CLB'),(13,'CYB'),(14,'CDB')]: near(h[-1][idx],r[key],'final monitor '+key)
        for key in ['CLS','CDS','CMS25']: near(r[key],se[key],key)
        near(r['CDS'],r['CDB']*math.cos(a)+r['CLB']*math.sin(a),'drag rotation')
        near(r['CLS'],-r['CDB']*math.sin(a)+r['CLB']*math.cos(a),'lift rotation')
        for stab,wind in [('CDS','CDW'),('CYS','CYW'),('CLS','CLW'),('CRS25','CRW25'),('CMS25','CMW25'),('CNS25','CNW25')]:near(r[stab],r[wind],'beta-zero frames')
        near(r['CMS25'],-.02-.1*r['CLS'],'10% static margin law')
        force=[list(map(float,l.split())) for l in (wing/f'total_force_state{s}').read_text().splitlines()]
        check(len(force)==9 and all(len(f)==4 for f in force),'headerless Y FX FY FZ')
        lift=[];drag=[];station=[]
        for f,sect in zip(force,se['sections']):
            y=sect['Y']; near(f[0],y,'station tolerance',1e-6)
            geo=curve(wing/f'section_state1_station{y:.3f}')
            xmin=min(x for x,z in geo);chord=max(x for x,z in geo)-xmin
            near(chord,sect['chord'],'chord');near(xmin,sect['Xmin'],'Xmin');check(chord>0,'positive chord')
            lp=f[1]*math.sin(a)-f[3]*math.cos(a)
            dp=-f[1]*math.cos(a)-f[3]*math.sin(a)
            near(lp/(q*chord),sect['cl'],'section cl');check(lp>0,'positive lift')
            near(2*y/bref,expected['geometry']['eta'][len(station)],'eta')
            lift.append(lp);drag.append(dp);station.append(y)
            cp=curve(wing/f'cp_dist_state{s}_station{y:.3f}')
            xs=[(x-xmin)/chord for x,v in cp]
            check(min(xs)>-1e-8 and max(xs)<1+1e-8,'Cp x/c range')
            upper=list(reversed(cp[:41]));lower=cp[40:]
            check(min(v for x,v in upper)<0,'upper suction')
            integral=trap([lo[1]-up[1] for up,lo in zip(upper,lower)],[(x-xmin)/chord for x,v in upper])
            near(integral,-f[3]/(q*chord),'Cp normal coefficient integral')
        near(trap(lift,station)/(q*sref),r['CLS'],'integrated span lift vs ADF')
        near(trap(drag,station)/(q*sref),r['CDS'],'integrated span drag vs ADF')
for name,entries in expected['drag_rise'].items():
    rows=table(BASE/'03-RESULTS/DRAG-RISE/BASELINE'/name)
    check(len(rows)==7,'seven Mach points')
    check(all(b['MACH']>a['MACH'] and b['CDS']>a['CDS'] for a,b in zip(rows,rows[1:])),'increasing Mach and CDS')
    for r,e in zip(rows,entries):
        near(r['CLS'],entries[0]['CLS'],'constant CLS')
        near(r['CDS']-rows[0]['CDS'],e['Delta_CD'],'drag rise')
print(f'PASS: {checks} checks; 2 POLARs, 10 states, 18 geometry curves, 90 Cp curves, 10 force tables, 2 drag-rise sweeps.')
