/* Additive UI: existing aerodynamic calculations remain in the v25 engine. */
function distributionSpanReference(series) {
  const value = series.bref;
  return Number.isFinite(value) && value > 0 ? value : null;
}
function distributionCoordinate(series, y, mode) {
  const reference = distributionSpanReference(series);
  return mode === 'y' ? y : reference === null ? null : y / reference * (mode === 'eta' ? 2 : 1);
}
function distributionChordLoad(point) {
  return Number.isFinite(point.cl) && Number.isFinite(point.chord) && point.chord > 0 ? point.cl * point.chord : null;
}
function distributionStation(series, station, mode='eta') {
  const value = distributionCoordinate(series, station.y, mode);
  const position = mode === 'y' ? `Y = ${station.y.toFixed(3)} m`
    : value === null ? `Y = ${station.y.toFixed(3)} m · BREF unavailable`
    : `${mode === 'eta' ? '2Y/BREF' : 'Y/BREF'} = ${value > 0 ? '+' : ''}${Number((100*value).toFixed(2))}%`;
  // eta and yb use the same selection identity; dimensional stations group by Y.
  const fraction = distributionCoordinate(series, station.y, 'yb');
  const dimensional = mode === 'y' || fraction === null;
  return {key:JSON.stringify([series.component, dimensional ? 'Y' : 'fraction',
    (dimensional ? station.y : fraction).toFixed(7)]), component:series.component,
    order:dimensional ? station.y : fraction, label:`${series.component} · ${position}`};
}
function distributionLoadsText(series, mode) {
  const field = value => value === null || value === undefined || (typeof value === 'number' && !Number.isFinite(value))
    ? 'NA' : String(value).replace(/[\t\r\n]/g,' ');
  const rows = ['# JAMAL spanwise loads; tab-separated; missing values = NA',
    '# Normalization: infout BREF. cl.c = cl * local chord [m]. All stations of displayed series.',
    `# Selected coordinate: ${mode === 'eta' ? '2Y/BREF' : mode === 'yb' ? 'Y/BREF' : 'Y [m]'} (normalized coordinates are fractions)`,
    'configuration\tpolar\tcomponent\tstate\talpha_deg\tbeta_deg\tY_m\tBREF_m\tcoordinate\tchord_m\tcl\tcl.c_m'];
  series.forEach(s => s.span.forEach(p => rows.push([s.configuration,s.polar,s.component,s.state,s.alpha,s.beta,
    p.y,distributionSpanReference(s),distributionCoordinate(s,p.y,mode),p.chord,p.cl,distributionChordLoad(p)].map(field).join('\t'))));
  return rows.join('\r\n')+'\r\n';
}
(() => {
  const byId = id => document.getElementById(id);
  const data = distributionData.series || [];
  const overlays = new Set();
  let ready = false;
  const selectedStations = new Set();
  const storageKey = 'JAMAL_distributions_v3:' + JSON.stringify(distinctSources());
  const panels = new Map();
  function distinctSources() {return [...new Set((distributionData.sources||[]).map(s => s.path.split(/DISTCLCP/i)[0]))].sort();}
  const key = s => JSON.stringify([s.configuration, s.polar, s.component, s.state]);
  const text = s => `${s.configuration} · ${s.polar} · ${s.component} · α=${s.alpha}° β=${s.beta}°`;
  const options = (id, values, label = x => x) => {
    const select = byId(id), old = select.value;
    select.replaceChildren(...values.map(value => new Option(label(value), String(value))));
    if (values.some(v => String(v) === old)) select.value = old;
  };
  const distinct = list => [...new Set(list)];
  const candidates = level => data.filter(s =>
    (level < 1 || s.configuration === byId('distConfig').value) &&
    (level < 2 || s.polar === byId('distPolar').value) &&
    (level < 3 || s.component === byId('distComponent').value));
  const selected = () => candidates(3).find(s => s.state === Number(byId('distState').value));
  const shown = () => data.filter(s => overlays.has(key(s)) || s === selected());
  const station = (s,p) => distributionStation(s,p,byId('distCoordinate').value);
  const stationCatalog = () => {
    const entries = new Map();
    shown().forEach(s => s.cp.forEach(p => {const item=station(s,p);entries.set(item.key,item);}));
    return [...entries.values()].sort((a,b)=>a.component.localeCompare(b.component)||a.order-b.order);
  };
  function remember() {
    if(!ready) return;
    try {localStorage.setItem(storageKey,JSON.stringify({
      configuration:byId('distConfig').value,polar:byId('distPolar').value,component:byId('distComponent').value,
      state:byId('distState').value,coordinate:byId('distCoordinate').value,leadingEdge:byId('distLeadingEdge').value,
      overlays:[...overlays],stations:[...selectedStations],single:byId('distSingle').checked
    }));} catch (_) {}
  }
  function updateControls(level = 0) {
    if (level <= 0) options('distConfig', distinct(data.map(s => s.configuration)));
    if (level <= 1) options('distPolar', distinct(candidates(1).map(s => s.polar)));
    if (level <= 2) options('distComponent', distinct(candidates(2).map(s => s.component)));
    options('distState', candidates(3).map(s => s.state), state => {
      const s = candidates(3).find(s => s.state === state);
      return `α=${s.alpha}° · β=${s.beta}°`;
    });
    updateStations(true); draw();
  }
  function updateStations(autoSelect=false) {
    const catalog = stationCatalog();
    const current = selected();
    const currentKeys = current ? current.cp.map(p=>station(current,p).key) : [];
    if ((!ready && !selectedStations.size) || (autoSelect && !currentKeys.some(k=>selectedStations.has(k)))) {
      if(byId('distSingle').checked) selectedStations.clear();
      currentKeys.forEach((k,i)=>{if(!byId('distSingle').checked || i===0) selectedStations.add(k);});
    }
    byId('distStations').replaceChildren(...catalog.map(item => {
      const label = document.createElement('label'), checkbox = document.createElement('input');
      checkbox.type = 'checkbox'; checkbox.value = item.key; checkbox.checked = selectedStations.has(item.key);
      checkbox.addEventListener('change', () => {
        if (byId('distSingle').checked) selectedStations.clear();
        if (checkbox.checked) selectedStations.add(item.key); else selectedStations.delete(item.key);
        updateStations(); draw();
      });
      label.append(checkbox, document.createTextNode(` ${item.label}`));
      return label;
    }));
  }
  function draw() {
    if (!ready || !byId('section-distributions').classList.contains('active')) return;
    const series = shown(), mode = byId('distCoordinate').value, le = byId('distLeadingEdge').value;
    const span = [], chordSpan = [], groups = new Map();
    let missingBref = false;
    const safe = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
    series.forEach((s, i) => {
      const style = traceStyle(s.configuration, s.polar);
      const name = safe(text(s)), shortName = safe(`${s.configuration} · ${s.polar.replace('POLAR-','P')} · ${s.component} · α=${s.alpha}° β=${s.beta}°`), bref = distributionSpanReference(s);
      const normalized = bref !== null;
      if (!normalized && mode !== 'y') missingBref = true;
      const x = y => mode === 'y' ? y : normalized ? y / bref * (mode === 'eta' ? 2 : 1) : null;
      span.push({x:s.span.map(p => x(p.y)), y:s.span.map(p => p.cl), mode:'lines+markers', name:shortName,
        line:{color:style.color,dash:POLAR_DASHES[i % POLAR_DASHES.length]},
        marker:{symbol:MARKER_SYMBOLS[i % MARKER_SYMBOLS.length]}, connectgaps:false,
        hovertemplate:`${name}<br>Span: %{x:.4f}<br>cl: %{y:.5f}<extra></extra>`});
      chordSpan.push({...span[span.length-1], y:s.span.map(distributionChordLoad),
        hovertemplate:`${name}<br>Span: %{x:.4f}<br>cl.c: %{y:.5f} m<extra></extra>`});
      s.cp.filter(p => selectedStations.has(station(s,p).key)).forEach(p => {
        const item = station(s,p);
        if(!groups.has(item.key)) groups.set(item.key,{...item,traces:[]});
        const eta = normalized ? (2*p.y/bref).toFixed(4) : 'unavailable';
        groups.get(item.key).traces.push({x:le === 'raw' ? p.x : le === 'xmax' ? p.xc.map(v => 1-v) : p.xc,
          y:p.values, mode:'lines', name:safe(`${s.configuration} · ${s.polar.replace('POLAR-','P')} · α=${s.alpha}° β=${s.beta}°`),
          line:{color:style.color,dash:POLAR_DASHES[i % POLAR_DASHES.length]},
          hoverlabel:{bgcolor:'rgba(0,0,0,0)',bordercolor:'rgba(0,0,0,0)',font:{color:style.color}},
          hovertemplate:`${name}<br>${safe(item.label)}<br>Y=${p.y.toFixed(3)} m · 2Y/BREF=${eta}<br>${le === 'raw' ? 'X' : 'x/c'}: %{x:.4f}<br>Cp: %{y:.5f}<extra></extra>`});
      });
    });
    const layout = (id,title,count) => {
      const columns = Math.max(1,Math.floor((byId(id).clientWidth-90)/280));
      const bottom = 90+Math.ceil(count/columns)*28;
      byId(id).style.height = `${360+bottom}px`;
      return {title,height:360+bottom,margin:{l:65,r:25,t:55,b:bottom},uirevision:'distributions',
        showlegend:true,legend:{orientation:'h',y:-.28,yanchor:'top',x:0,entrywidth:260}, font:{family:'Arial, Helvetica, sans-serif'}};
    };
    Plotly.react('distSpanPlot', span, mixedTypographyLayout({...layout('distSpanPlot','Sectional cl',span.length),
      xaxis:{title:mode==='eta'?'2Y/BREF':mode==='yb'?'Y/BREF':'Y [m]'},yaxis:{title:'cl'},uirevision:mode}), {responsive:true});
    Plotly.react('distSpanChordPlot', chordSpan, mixedTypographyLayout({...layout('distSpanChordPlot','Sectional cl.c',chordSpan.length),
      xaxis:{title:mode==='eta'?'2Y/BREF':mode==='yb'?'Y/BREF':'Y [m]'},yaxis:{title:'cl.c [m]'},uirevision:mode}), {responsive:true});
    byId('distExport').disabled = !series.some(s=>s.span.length);
    // Reuse one Plotly surface per station; release removed panels and their listeners.
    for(const [key,panel] of panels) if(!groups.has(key)) {
      Plotly.purge(panel.plot);panel.card.remove();panels.delete(key);
    }
    const ordered = [...groups.values()].sort((a,b)=>a.component.localeCompare(b.component)||a.order-b.order);
    let low=Infinity, high=-Infinity, curveCount=0;
    ordered.forEach(g=>g.traces.forEach(t=>{curveCount++;t.y.forEach(v=>{if(Number.isFinite(v)){low=Math.min(low,v);high=Math.max(high,v);}});}));
    const padding=Number.isFinite(low)?Math.max(.05,(high-low)*.08):.05;
    const cpRange=Number.isFinite(low)?[high+padding,low-padding]:[1,-1];
    ordered.forEach(g=>{
      let panel=panels.get(g.key);
      if(!panel){
        const card=document.createElement('article'), heading=document.createElement('h4'), plot=document.createElement('div');
        card.className='dist-cp-panel';plot.className='plot';card.append(heading,plot);
        panel={card,heading,plot};panels.set(g.key,panel);
      }
      panel.heading.textContent=g.label;
      byId('distCpGrid').append(panel.card);
      const bottom=85+g.traces.length*22, height=300+bottom;
      panel.plot.style.height=`${height}px`;
      Plotly.react(panel.plot,g.traces,mixedTypographyLayout({height,margin:{l:45,r:12,t:15,b:bottom},
        xaxis:{title:le==='raw'?'X [m]':'x/c',range:le==='raw'?undefined:[0,1]},
        yaxis:{title:'Cp',range:cpRange,autorange:false},
        showlegend:true,legend:{orientation:'h',y:-.28,yanchor:'top',font:{size:10}},
        uirevision:JSON.stringify([le,g.key,cpRange])}),{responsive:true,displaylogo:false});
    });
    byId('distOverlays').textContent = overlays.size ? `Pinned overlays: ${data.filter(s => overlays.has(key(s))).map(text).join(' | ')}` : 'Current state shown. Add overlay to retain it while selecting another state or configuration.';
    byId('distStatus').textContent = !data.length ? 'No distribution data found for the selected POLARs. Expected: 03-RESULTS/DISTCLCP/POLAR-XXX/<component>.' :
      `${series.length} state(s) · ${groups.size} station plot(s) · ${curveCount} Cp curve(s). ${missingBref ? 'Missing positive span reference: choose dimensional Y. ' : ''}${curveCount ? '' : 'Select a station to show Cp. '}${(distributionData.issues||[]).length} distribution integrity warning(s).`;
    byId('distReference').textContent = (mode==='y' ? 'Station labels use dimensional Y [m]. ' :
      `Station percentage = ${mode==='eta'?'200':'100'} × Y / infout BREF. `)+
      distinct(series.map(s=>`${s.configuration} · ${s.polar}: BREF = ${distributionSpanReference(s) ?? 'unavailable'} m`)).join(' · ');
    remember();
  }
  function init() {
    updateControls();
    try {
      const saved = JSON.parse(localStorage.getItem(storageKey) || 'null');
      if(saved) {
        const restore = (id, value) => {if([...byId(id).options].some(o=>o.value===String(value))) byId(id).value=String(value);};
        restore('distConfig',saved.configuration); updateControls(1);
        restore('distPolar',saved.polar); updateControls(2);
        restore('distComponent',saved.component); updateControls(3);
        restore('distState',saved.state); restore('distCoordinate',saved.coordinate); restore('distLeadingEdge',saved.leadingEdge);
        (saved.overlays||[]).filter(k=>data.some(s=>key(s)===k)).forEach(k=>overlays.add(k));
        selectedStations.clear(); (saved.stations||[]).forEach(y=>selectedStations.add(y));
        byId('distSingle').checked=Boolean(saved.single);
      }
    } catch (_) {}
    ready = true; updateStations();
    byId('distIntegrity').textContent = [
      ...(distributionData.issues||[]).map(i => `${i.configuration} · ${i.polar} · ${i.details}`),
      '\nSources (size and modification timestamp retained in dashboard.json):',
      ...(distributionData.sources||[]).map(s => s.path)
    ].join('\n');
    ['distConfig','distPolar','distComponent'].forEach((id,i) => byId(id).addEventListener('change', () => updateControls(i+1)));
    byId('distState').addEventListener('change', () => {updateStations(); draw();});
    let previousMode = byId('distCoordinate').value;
    byId('distCoordinate').addEventListener('change', () => {
      const mode = byId('distCoordinate').value, retained = new Set();
      shown().forEach(s=>s.cp.forEach(p=>{
        if(selectedStations.has(distributionStation(s,p,previousMode).key)) retained.add(distributionStation(s,p,mode).key);
      }));
      selectedStations.clear(); retained.forEach(k=>selectedStations.add(k));
      previousMode=mode; updateStations(); draw();
    });
    byId('distLeadingEdge').addEventListener('change', draw);
    byId('distExport').onclick = () => {
      const series=shown(), blob=new Blob([distributionLoadsText(series,byId('distCoordinate').value)],{type:'text/plain;charset=utf-8'});
      const url=URL.createObjectURL(blob), link=document.createElement('a');
      link.href=url;link.download='JAMAL_spanwise_loads.txt';document.body.append(link);link.click();link.remove();
      setTimeout(()=>URL.revokeObjectURL(url),10000);
      byId('distExportStatus').textContent=`Download requested: ${series.reduce((n,s)=>n+s.span.length,0)} stations from ${series.length} displayed series.`;
    };
    byId('distAdd').onclick = () => {if(selected()) overlays.add(key(selected())); updateStations(); draw();};
    byId('distClear').onclick = () => {overlays.clear(); updateStations(); draw();};
    function pick(n) {
      selectedStations.clear();
      [...byId('distStations').querySelectorAll('input')].forEach((el,i) => {
        if (n && i%n===0 && (!byId('distSingle').checked || !selectedStations.size)) selectedStations.add(el.value);
      });
      updateStations(); draw();
    }
    byId('distAll').onclick = () => pick(1);
    byId('distNone').onclick = () => pick(0);
    byId('distEvery').onclick = () => pick(Math.max(1, Math.floor(Number(byId('distNth').value)||1)));
    byId('distSingle').onchange = () => {if(byId('distSingle').checked && selectedStations.size>1) {
      const first = selectedStations.values().next().value; selectedStations.clear(); selectedStations.add(first); updateStations(); draw();
    }};
  }
  const originalShow = showSection;
  showSection = function(section, button) {
    originalShow(section, button);
    byId('analysisToolbar').style.display = section === 'distributions' ? 'none' : '';
    byId('advancedPanel').style.display = section === 'distributions' ? 'none' : '';
    if(section === 'distributions') {if(!ready) init(); draw();}
  };
})();
