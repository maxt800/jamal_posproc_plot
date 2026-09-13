"""Deterministic additive paired-V-tail fixture; original wing/run/ADF files are untouched."""
from pathlib import Path
import hashlib
import json
import math
import sys


def trapezoid(values, coordinates):
    return sum((a+b)*(d-c)/2 for a,b,c,d in zip(values,values[1:],coordinates,coordinates[1:]))


def write_curve(path, points, title):
    header = ('*KEYWORD\n*DEFINE_CURVE_TITLE\n' + title + '\n'
              '$ LCID SIDR SFA SFO OFFA OFFO DATTYP\n'
              ' 1 0 1.0 1.0 0 0 0\n$ ABSCISSA ORDINATE\n')
    path.write_text(header + ''.join(f'{x: .12E} {y: .12E}\n' for x,y in points) + '*END\n', encoding='utf-8')


def generate(base):
    base = Path(base)
    wing = json.loads((base/'expected_values.json').read_text(encoding='utf-8'))
    span = wing['geometry']['BREF']/2
    cant = math.radians(45)
    stations = [span*(i-4)/8 for i in range(9)]
    x_ascending = [(1-math.cos(math.pi*i/40))/2 for i in range(41)]
    shape = [math.sqrt(x)*(1-x) for x in x_ascending]
    shape_integral = trapezoid(shape, x_ascending)
    expected = {'synthetic': True, 'component':'VTAIL', 'projected_span':span,
                'wing_span':wing['geometry']['BREF'], 'cant_degrees':45,
                'projected_area':4.5, 'actual_panel_span_total':span/math.cos(cant), 'polars':{}}
    written = []
    for polar, polar_info in wing['polars'].items():
        root = base/'03-RESULTS/DISTCLCP'/polar/'VTAIL'
        root.mkdir(parents=True, exist_ok=True)
        reference = {'schema_version':1, 'span_reference_m':span, 'synthetic':True,
                     'geometry':'paired V-tail', 'cant_degrees':45,
                     'station_coordinate':'projected Y [m]; centerline 0',
                     'panel_center_z':'abs(Y)*tan(45 degrees)',
                     'curve_ordinate':'local section-normal coordinate [m]',
                     'force_density':'body forces per projected Y metre [N/m]',
                     'notes':'Demonstration component only; not included in the existing whole-aircraft ADF coefficients.'}
        path=root/'component_reference.json'
        path.write_text(json.dumps(reference,indent=2)+'\n',encoding='utf-8');written.append(path)
        expected['polars'][polar]={'qdin':polar_info['qdin'],'states':[]}
        for y in stations:
            eta=2*y/span
            chord=1.2-.6*abs(eta)
            xmin=7+.25*abs(y)
            points=[]
            for indices, sign in [(range(40,-1,-1),1),(range(1,41),-1)]:
                for i in indices:
                    x=x_ascending[i]
                    thickness=5*.10*(.2969*math.sqrt(x)-.126*x-.3516*x*x+.2843*x**3-.1036*x**4)
                    points.append((xmin+chord*x, sign*chord*thickness))
            path=root/f'section_state1_station{y:.3f}'
            write_curve(path,points,'SYNTHETIC VTAIL local section geometry');written.append(path)
        for state in polar_info['states']:
            alpha=math.radians(state['alpha']);q=polar_info['qdin']
            sections=[];force_rows=[]
            for y in stations:
                eta=2*y/span;chord=1.2-.6*abs(eta);xmin=7+.25*abs(y)
                cn=(.18+.035*state['alpha'])*(1+.1*(polar_info['Mach']-.2))*(1-.4*eta*eta)
                cd=.012+.025*cn*cn
                # Panel-span force -> projected-Y density includes ds/dY = 1/cos(cant).
                normal_per_y=q*chord*cn/math.cos(cant)
                fx=-q*chord*cd/math.cos(cant)
                fy=0 if y==0 else -math.copysign(1,y)*normal_per_y*math.sin(cant)
                fz=-normal_per_y*math.cos(cant)
                cl=(fx*math.sin(alpha)-fz*math.cos(alpha))/(q*chord)
                upper=[];lower=[]
                for x,weight in zip(x_ascending,shape):
                    recovery=math.exp(-60*x)*(1-x)
                    upper.append(recovery-.5*cn*weight/shape_integral)
                    lower.append(recovery+.5*cn*weight/shape_integral)
                cp=[(xmin+chord*x_ascending[i],upper[i]) for i in range(40,-1,-1)]
                cp += [(xmin+chord*x_ascending[i],lower[i]) for i in range(1,41)]
                path=root/f"cp_dist_state{state['state']}_station{y:.3f}"
                write_curve(path,cp,'SYNTHETIC VTAIL Cp');written.append(path)
                stored_y=y-2e-7 if math.isclose(y,span/4) else y
                force_rows.append((stored_y,fx,fy,fz))
                sections.append({'Y':y,'station_percent':100*y/span,'Z_center':abs(y)*math.tan(cant),
                                 'chord':chord,'Xmin':xmin,'normal_coefficient':cn,'cl':cl,
                                 'FX':fx,'FY':fy,'FZ':fz,'Cp_upper_min':min(upper),
                                 'Cp_normal_integral':trapezoid([b-a for a,b in zip(upper,lower)],x_ascending)})
            path=root/f"total_force_state{state['state']}"
            path.write_text(''.join(' '.join(f'{v:.12E}' for v in row)+'\n' for row in force_rows),encoding='utf-8');written.append(path)
            expected['polars'][polar]['states'].append({'state':state['state'],'alpha':state['alpha'],
                                                       'beta':state['beta'],'sections':sections})
    path=base/'expected_vtail_values.json'
    path.write_text(json.dumps(expected,indent=2)+'\n',encoding='utf-8');written.append(path)
    path=base/'VTAIL_SHA256SUMS.txt'
    path.write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(base).as_posix()}\n'
                            for p in sorted(written)),encoding='ascii')
    return expected


if __name__=='__main__':
    base=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'JAMAL_SYNTHETIC_CFD'
    result=generate(base)
    print(f"Created paired VTAIL: {result['projected_span']:g} m projected span, 45 degree cant, 2 POLARs, 10 states, 9 stations/state.")
