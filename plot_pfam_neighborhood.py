#!/usr/bin/env python3
"""
Plot PFAM neighborhood around a reference PFAM domain.

Usage:
    python plot_pfam_neighborhood.py <tab_file> <pfam_id> <window_bp> <ratio_cutoff>

Example:
    python plot_pfam_neighborhood.py data.tab PF05014.22 1000 10
"""

import re, json, colorsys, sys, os
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def parse_args():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    tab_file   = sys.argv[1]
    pfam_id    = sys.argv[2]   # e.g. PF05014.22  (prefix match used)
    window     = int(sys.argv[3])
    ratio_cut  = float(sys.argv[4])
    return tab_file, pfam_id, window, ratio_cut


def load_data(tab_file):
    data = defaultdict(list)
    with open(tab_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 19:
                continue
            contig_full = parts[0]
            m = re.match(r'^(.+?_ratio_[\d.]+)-\d+[FR]$', contig_full)
            if not m:
                continue
            base = m.group(1)
            data[base].append({
                'pfam_name': parts[3],
                'pfam_id':   parts[4],
                'start':     int(parts[17]),
                'end':       int(parts[18]),
            })
    return data


def get_ratio(base):
    m = re.search(r'ratio_([\d.]+)', base)
    return float(m.group(1)) if m else 0.0


def get_neighbors(anns, ref_pfam_prefix, window):
    ref = next((a for a in anns if a['pfam_id'].startswith(ref_pfam_prefix)), None)
    if not ref:
        return []
    return [a for a in anns
            if a['end'] >= ref['start'] - window and a['start'] <= ref['end'] + window]


def get_ref_mid(neighbors, ref_pfam_prefix):
    for a in neighbors:
        if a['pfam_id'].startswith(ref_pfam_prefix):
            return (a['start'] + a['end']) // 2
    return 0


def assign_colors(pfam_ids):
    pfam_ids = sorted(pfam_ids)
    n = max(len(pfam_ids), 1)
    colors = {}
    for i, pid in enumerate(pfam_ids):
        rc, gc, bc = colorsys.hsv_to_rgb(i / n, 0.70, 0.85)
        colors[pid] = f'rgb({int(rc*255)},{int(gc*255)},{int(bc*255)})'
    return colors


def build_js_data(data, ref_pfam_id, window, ratio_cut):
    # Use prefix up to the dot for matching (e.g. "PF05014" from "PF05014.22")
    ref_prefix = ref_pfam_id.split('.')[0]

    pf_contigs = {b: anns for b, anns in data.items()
                  if any(a['pfam_id'].startswith(ref_prefix) for a in anns)}

    if not pf_contigs:
        print(f"No contigs found with PFAM ID starting with '{ref_prefix}'.", file=sys.stderr)
        sys.exit(1)

    rows = []
    high_pfams = set()

    for base, anns in pf_contigs.items():
        ratio = get_ratio(base)
        neighbors = get_neighbors(anns, ref_prefix, window)
        if ratio > ratio_cut:
            for a in neighbors:
                high_pfams.add(a['pfam_id'])
        rows.append({'base': base, 'ratio': ratio, 'neighbors': neighbors})

    rows.sort(key=lambda x: -x['ratio'])

    above = [r for r in rows if r['ratio'] >  ratio_cut]
    below = [r for r in rows if r['ratio'] <= ratio_cut]
    n_above, n_below = len(above), len(below)

    all_pfam_ids = sorted({a['pfam_id'] for r in rows for a in r['neighbors']})
    colors = assign_colors(all_pfam_ids)

    pfam_name_map = {a['pfam_id']: a['pfam_name']
                     for r in rows for a in r['neighbors']}

    # Fisher's exact test: for each PFAM, presence/absence in above vs below
    def pfam_presence(row_list, pid):
        return sum(1 for r in row_list if any(a['pfam_id'] == pid for a in r['neighbors']))

    candidate_pfams = sorted(high_pfams)  # only test PFAMs that appear in high-ratio contigs
    pvals_raw = []
    for pid in candidate_pfams:
        a_yes = pfam_presence(above, pid)
        b_yes = pfam_presence(below, pid)
        a_no  = n_above - a_yes
        b_no  = n_below - b_yes
        _, p  = fisher_exact([[a_yes, a_no], [b_yes, b_no]], alternative='greater')
        pvals_raw.append(p)

    n_tests = len(candidate_pfams)
    # BH/FDR correction
    if n_tests > 0:
        reject, pvals_fdr, _, _ = multipletests(pvals_raw, alpha=0.05, method='fdr_bh')
    else:
        reject, pvals_fdr = [], []

    enriched_pfams  = {}  # pid -> {'p_raw', 'p_fdr', 'sig_fdr', 'sig_nominal'}
    for pid, p_raw, p_fdr, rej in zip(candidate_pfams, pvals_raw, pvals_fdr, reject):
        enriched_pfams[pid] = {
            'p_raw': p_raw, 'p_fdr': p_fdr,
            'sig_fdr':     bool(rej),
            'sig_nominal': p_raw < 0.05 and not rej,
        }

    sig_fdr     = {pid for pid, v in enriched_pfams.items() if v['sig_fdr']}
    sig_nominal = {pid for pid, v in enriched_pfams.items() if v['sig_nominal']}
    sig_pfams   = sig_fdr | sig_nominal   # used for track highlight

    print(f"PFAMs tested   : {n_tests}  (FDR q<0.05: {len(sig_fdr)},  nominal p<0.05: {len(sig_nominal)})")
    for pid in candidate_pfams:
        v = enriched_pfams[pid]
        if v['sig_fdr'] or v['sig_nominal']:
            tag = '★ FDR' if v['sig_fdr'] else '◆ nominal'
            print(f"  {pid:20s}  {pfam_name_map.get(pid,'?'):30s}  p={v['p_raw']:.2e}  FDR={v['p_fdr']:.2e}  {tag}")

    def _pval_sort_key(pid):
        v = enriched_pfams.get(pid, {'p_raw': 1.0, 'p_fdr': 1.0, 'sig_fdr': False, 'sig_nominal': False})
        tier = 0 if v['sig_fdr'] else (1 if v['sig_nominal'] else 2)
        return (tier, v['p_raw'])

    legend_entries = [
        {'id':          pid,
         'name':        pfam_name_map.get(pid, pid),
         'color':       colors[pid],
         'p_raw':       round(enriched_pfams.get(pid, {}).get('p_raw', 1.0), 6),
         'p_fdr':       round(enriched_pfams.get(pid, {}).get('p_fdr', 1.0), 6),
         'sig_fdr':     pid in sig_fdr,
         'sig_nominal': pid in sig_nominal,
         'enriched':    pid in sig_pfams}
        for pid in sorted(high_pfams)
    ]
    legend_entries.sort(key=lambda x: (
        not x['sig_fdr'], not x['sig_nominal'],
        not x['id'].startswith(ref_prefix),
        x['p_raw'], x['name']
    ))

    def build_row(r):
        ref_mid = get_ref_mid(r['neighbors'], ref_prefix)
        return {
            'base':  r['base'],
            'ratio': round(r['ratio'], 4),
            'anns':  [{'pfam_id': a['pfam_id'], 'pfam_name': a['pfam_name'],
                       'start': a['start'] - ref_mid, 'end': a['end'] - ref_mid}
                      for a in r['neighbors']],
        }

    print(f"Reference PFAM : {ref_pfam_id}  (prefix: {ref_prefix})")
    print(f"Window         : ±{window} bp")
    print(f"Ratio cutoff   : {ratio_cut}")
    print(f"Contigs found  : {len(pf_contigs)}  (above cutoff: {len(above)}, below: {len(below)})")

    return {
        'above':      [build_row(r) for r in above],
        'below':      [build_row(r) for r in below],
        'high_pfams':   list(high_pfams),
        'sig_fdr':      list(sig_fdr),
        'sig_nominal':  list(sig_nominal),
        'sig_pfams':    list(sig_pfams),
        'colors':     colors,
        'legend':     legend_entries,
        'ref_pfam':   ref_pfam_id,
        'window':     window,
        'ratio_cut': ratio_cut,
    }


def render_html(js_data):
    ref_pfam   = js_data['ref_pfam']
    window     = js_data['window']
    ratio_cut  = js_data['ratio_cut']
    n_above    = len(js_data['above'])
    n_below    = len(js_data['below'])
    data_json  = json.dumps(js_data)
    ref_prefix = ref_pfam.split('.')[0]

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{ref_pfam} Neighborhood</title>
<style>
  body {{ font-family: monospace; font-size: 11px; background: #fff; margin: 0; padding: 10px; }}
  h1   {{ font-size: 14px; margin: 6px 0; }}
  h2   {{ font-size: 12px; margin: 8px 0 4px 0; }}
  .section {{ border: 1px solid #888; padding: 6px; margin-bottom: 16px; }}
  .row {{ display: flex; align-items: center; margin: 2px 0; height: 16px; }}
  .label {{ width: 420px; min-width: 420px; text-align: right; padding-right: 6px;
             white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
             font-size: 9px; color: #333; }}
  svg.track {{ display: block; }}
  #tooltip {{
    position: fixed; background: rgba(0,0,0,0.82); color: #fff;
    padding: 6px 10px; border-radius: 4px; font-size: 11px;
    pointer-events: none; display: none; max-width: 340px; z-index: 999; line-height: 1.5;
  }}
  #legend {{ border: 1px solid #ccc; padding: 8px 12px; margin-bottom: 12px; background: #fafafa; }}
  #legend h2 {{ margin-top: 0; }}
  #legend-grid {{ display: flex; flex-wrap: wrap; gap: 4px 14px; }}
  .legend-item {{ display: flex; align-items: center; gap: 5px; cursor: default;
                  font-size: 10px; white-space: nowrap; }}
  .legend-swatch {{ width: 14px; height: 10px; border-radius: 2px; border: 1px solid #555; flex-shrink: 0; }}
  .legend-item.highlight .legend-label {{ font-weight: bold; text-decoration: underline; }}
  .legend-item.enriched .legend-label  {{ color: #c00; }}
  .legend-item.enriched .legend-swatch {{ box-shadow: 0 0 0 2px #222; }}
  .enriched-star {{ color: #c00; margin-left: 3px; font-size: 11px; }}
  .pval-badge {{ font-size: 9px; color: #666; margin-left: 4px; }}
  #dl-bar {{ position:fixed; top:10px; right:14px; display:flex; gap:6px; z-index:200; }}
  #dl-bar button {{
    padding:4px 10px; font-size:11px; font-family:monospace; cursor:pointer;
    border:1px solid #888; border-radius:4px; background:#f5f5f5;
  }}
  #dl-bar button:hover {{ background:#e0e0e0; }}
  @media print {{
    #dl-bar, #tooltip {{ display:none !important; }}
    .label {{ display:none !important; }}
    body {{ padding:0; }}
  }}
</style>
</head>
<body>
<div id="dl-bar">
  <button onclick="downloadPNG()">&#8659; PNG</button>
  <button onclick="window.print()">&#8659; PDF</button>
</div>
<div id="tooltip"></div>
<h1>{ref_pfam} Neighborhood &plusmn;{window} bp</h1>
<p style="font-size:10px;color:#555;margin:2px 0 8px">
  Colored fill = PFAM found in ratio&gt;{ratio_cut} contigs. Outline only = not found in high-ratio contigs.
  <b>Hover</b> on any box to see name and overlapping annotations.
</p>
<div id="legend">
  <h2>Legend &mdash; PFAMs present in ratio &gt; {ratio_cut} contigs
    &nbsp;<span style="font-size:10px;font-weight:normal;color:#555">
    <span style="color:#c00">&#9733;</span> FDR q&lt;0.05 (BH) &nbsp;|&nbsp;
    <span style="color:#e07000">&#9670;</span> nominal p&lt;0.05 (dashed border on plot)
    </span>
  </h2>
  <div id="legend-grid"></div>
</div>
<div id="app"></div>

<script>
const DATA = {data_json};

const PLOT_WIDTH = 600;
const DOMAIN     = {window} + 300;   // a little margin beyond the window
const SCALE      = PLOT_WIDTH / (2 * DOMAIN);
const MID        = PLOT_WIDTH / 2;
const ROW_H      = 14;
const ANN_H      = 10;
const ANN_Y      = 2;
const REF_PREFIX = {json.dumps(ref_prefix)};

function bpToPx(bp) {{ return MID + bp * SCALE; }}

const highPfams   = new Set(DATA.high_pfams);
const sigFdr      = new Set(DATA.sig_fdr);
const sigNominal  = new Set(DATA.sig_nominal);
const sigPfams    = new Set(DATA.sig_pfams);
const colors      = DATA.colors;

// Legend
const legendGrid = document.getElementById('legend-grid');
DATA.legend.forEach(entry => {{
  const item   = document.createElement('div');
  item.className = 'legend-item';
  item.dataset.pfamId = entry.id;

  const swatch = document.createElement('div');
  swatch.className = 'legend-swatch';
  swatch.style.background   = entry.id.startsWith(REF_PREFIX) ? '#e63946' : entry.color;
  swatch.style.borderColor  = entry.id.startsWith(REF_PREFIX) ? '#c1121f' : '#555';

  if (entry.enriched) item.classList.add('enriched');

  const pfamBase = entry.id.split('.')[0];
  const lbl = document.createElement('a');
  lbl.className = 'legend-label';
  lbl.textContent = entry.name + ' (' + entry.id + ')';
  lbl.href = 'https://www.ebi.ac.uk/interpro/entry/pfam/' + pfamBase + '/';
  lbl.target = '_blank';
  lbl.style.cssText = 'color:inherit;text-decoration:none;';
  lbl.addEventListener('mouseenter', () => lbl.style.textDecoration = 'underline');
  lbl.addEventListener('mouseleave', () => lbl.style.textDecoration = 'none');

  item.appendChild(swatch);
  item.appendChild(lbl);

  if (entry.sig_fdr || entry.sig_nominal) {{
    const star = document.createElement('span');
    star.className = 'enriched-star';
    if (entry.sig_fdr) {{
      star.textContent = '★';
      star.title = 'Significantly enriched (BH/FDR q < 0.05)';
      star.style.color = '#c00';
    }} else {{
      star.textContent = '◆';
      star.title = 'Nominally enriched (p < 0.05, not FDR-corrected)';
      star.style.color = '#e07000';
    }}
    item.appendChild(star);
    const badge = document.createElement('span');
    badge.className = 'pval-badge';
    badge.textContent = (entry.sig_fdr ? 'FDR=' + entry.p_fdr.toExponential(1) : 'p=' + entry.p_raw.toExponential(1));
    item.appendChild(badge);
  }}

  legendGrid.appendChild(item);
}});

// Tooltip
const tooltip = document.getElementById('tooltip');
document.addEventListener('mousemove', e => {{
  tooltip.style.left = (e.clientX + 14) + 'px';
  tooltip.style.top  = (e.clientY - 10) + 'px';
}});

function makeSVGEl(tag, attrs) {{
  const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  return el;
}}

function buildTrack(row, inAbove) {{
  const svg = makeSVGEl('svg', {{
    viewBox: `0 0 ${{PLOT_WIDTH}} ${{ROW_H}}`,
    width: PLOT_WIDTH, height: ROW_H
  }});
  svg.appendChild(makeSVGEl('line', {{x1:0, y1:ROW_H/2, x2:PLOT_WIDTH, y2:ROW_H/2, stroke:'#ddd', 'stroke-width':'0.5'}}));
  svg.appendChild(makeSVGEl('line', {{x1:MID, y1:0, x2:MID, y2:ROW_H, stroke:'#bbb', 'stroke-width':'0.8', 'stroke-dasharray':'2,2'}}));

  row.anns.forEach(a => {{
    const x1 = bpToPx(a.start), x2 = bpToPx(a.end);
    const w  = Math.max(x2 - x1, 2);
    const isRef  = a.pfam_id.startsWith(REF_PREFIX);
    const isHigh = highPfams.has(a.pfam_id);
    const fill   = isRef  ? '#e63946' : (isHigh ? (colors[a.pfam_id] || '#aaa') : 'none');
    const stroke = isRef  ? '#c1121f' : (isHigh ? (colors[a.pfam_id] || '#888') : '#888');

    const isSigFdr  = sigFdr.has(a.pfam_id);
    const isSigNom  = sigNominal.has(a.pfam_id);
    const rect = makeSVGEl('rect', {{
      x: x1, y: ANN_Y, width: w, height: ANN_H,
      fill, stroke,
      'stroke-width':    isRef ? '1.2' : (isSigFdr ? '2.5' : (isSigNom ? '1.8' : '0.8')),
      'stroke-dasharray': isSigNom && !isRef ? '3,1.5' : 'none',
      rx: '1', ry: '1', style: 'cursor:pointer',
      opacity: (isSigFdr || isSigNom || isRef) ? '1.0' : '0.6'
    }});

    const overlapping = row.anns.filter(b => b !== a && b.start < a.end && b.end > a.start);
    const overlapStr  = overlapping.length
      ? overlapping.map(b => `&bull; ${{b.pfam_name}} (${{b.pfam_id}})`).join('<br>')
      : 'none';

    rect.addEventListener('mouseenter', () => {{
      const enrichTag = isSigFdr ? ' <span style="color:#f44">★ FDR-enriched</span>'
                      : isSigNom ? ' <span style="color:#f90">◆ nominally enriched</span>' : '';
      tooltip.innerHTML = `<b>${{a.pfam_name}}</b> (${{a.pfam_id}})${{enrichTag}}<br>pos: ${{a.start}} &ndash; ${{a.end}}<br><span style="color:#ccc">overlapping:</span><br>${{overlapStr}}`;
      tooltip.style.display = 'block';
      document.querySelectorAll('.legend-item').forEach(el =>
        el.classList.toggle('highlight', el.dataset.pfamId === a.pfam_id));
    }});
    rect.addEventListener('mouseleave', () => {{
      tooltip.style.display = 'none';
      document.querySelectorAll('.legend-item').forEach(el => el.classList.remove('highlight'));
    }});
    svg.appendChild(rect);
  }});
  return svg;
}}

function buildAxisSVG() {{
  const svg = makeSVGEl('svg', {{viewBox:`0 0 ${{PLOT_WIDTH}} 20`, width:PLOT_WIDTH, height:20, style:'display:block'}});
  const ticks = [];
  for (let t = -DATA.window; t <= DATA.window; t += Math.ceil(DATA.window / 4 / 100) * 100)
    ticks.push(t);
  if (!ticks.includes(0)) ticks.push(0);
  ticks.forEach(t => {{
    const x = bpToPx(t);
    svg.appendChild(makeSVGEl('line', {{x1:x, y1:0, x2:x, y2:6, stroke:'#555', 'stroke-width':'1'}}));
    const lbl = makeSVGEl('text', {{x, y:16, 'text-anchor':'middle', 'font-size':'8', fill:'#333'}});
    lbl.textContent = t === 0 ? DATA.ref_pfam : `${{t>0?'+':''}}${{t}}`;
    svg.appendChild(lbl);
  }});
  svg.appendChild(makeSVGEl('line', {{x1:0, y1:0, x2:PLOT_WIDTH, y2:0, stroke:'#555', 'stroke-width':'0.8'}}));
  return svg;
}}

const BAR_W = 80;   // px width of the ratio bar column

function buildBarSVG(ratio, maxRatio) {{
  const svg = makeSVGEl('svg', {{viewBox:`0 0 ${{BAR_W}} ${{ROW_H}}`, width:BAR_W, height:ROW_H, style:'display:block'}});
  const barW = Math.min((ratio / maxRatio) * (BAR_W - 22), BAR_W - 22);
  svg.appendChild(makeSVGEl('rect', {{x:0, y:ANN_Y, width: Math.max(barW,1), height:ANN_H, fill:'#4a90d9', rx:'1'}}));
  const lbl = makeSVGEl('text', {{x: Math.max(barW,1) + 2, y: ANN_Y + ANN_H - 1,
    'font-size':'7', fill:'#333', 'dominant-baseline':'auto'}});
  lbl.textContent = ratio >= 1000 ? ratio.toExponential(1) : ratio.toFixed(1);
  svg.appendChild(lbl);
  return svg;
}}

function buildSection(rows, title, inAbove) {{
  const section = document.createElement('div');
  section.className = 'section';
  const h = document.createElement('h2');
  h.textContent = title;
  section.appendChild(h);

  const maxRatio = Math.max(...DATA.above.concat(DATA.below).map(r => r.ratio), 1);

  rows.forEach(row => {{
    const div  = document.createElement('div');
    div.className = 'row';
    const lbl  = document.createElement('div');
    lbl.className = 'label';
    lbl.textContent = row.base;
    lbl.title = row.base;
    const wrap = document.createElement('div');
    wrap.style.cssText = 'flex:1;position:relative';
    wrap.appendChild(buildTrack(row, inAbove));
    const barWrap = document.createElement('div');
    barWrap.style.cssText = `width:${{BAR_W}}px;min-width:${{BAR_W}}px;`;
    barWrap.appendChild(buildBarSVG(row.ratio, maxRatio));
    div.appendChild(lbl);
    div.appendChild(wrap);
    div.appendChild(barWrap);
    section.appendChild(div);
  }});

  const axisRow = document.createElement('div');
  axisRow.className = 'row';
  const axisLbl = document.createElement('div');
  axisLbl.className = 'label';
  axisLbl.textContent = 'bp from ' + DATA.ref_pfam + ' center';
  const axisWrap = document.createElement('div');
  axisWrap.style.cssText = 'flex:1';
  axisWrap.appendChild(buildAxisSVG());
  axisRow.appendChild(axisLbl);
  axisRow.appendChild(axisWrap);
  // bar column axis label
  const barAxisWrap = document.createElement('div');
  barAxisWrap.style.cssText = `width:${{BAR_W}}px;min-width:${{BAR_W}}px;font-size:8px;color:#555;padding-top:2px;`;
  barAxisWrap.textContent = 'ratio';
  axisRow.appendChild(barAxisWrap);
  section.appendChild(axisRow);
  return section;
}}

const app = document.getElementById('app');
app.appendChild(buildSection(DATA.above, `Ratio > ${{DATA.ratio_cut}}  (n=${{DATA.above.length}})`, true));
app.appendChild(buildSection(DATA.below, `Ratio ≤ ${{DATA.ratio_cut}}  (n=${{DATA.below.length}})`, false));

// ── PNG download ──────────────────────────────────────────────────────────────
async function downloadPNG() {{
  const rows = [...document.querySelectorAll('.row')];

  const GAP    = 2;
  const rowH   = 16;
  const trackW = 600;
  const barW   = 80;
  const totalW = trackW + barW + 20;

  // Legend layout: 3 columns
  const legEntries = DATA.legend;
  const LEG_COLS   = 3;
  const LEG_ROW_H  = 14;
  const SWATCH_W   = 12, SWATCH_H = 9;
  const legRows    = Math.ceil(legEntries.length / LEG_COLS);
  const legendH    = legRows * LEG_ROW_H + 28;

  const plotH  = rows.length * (rowH + GAP) + 4;
  const totalH = 40 + plotH + legendH + 10;

  // Render at a high pixel density for publication-quality output, capped so
  // we never exceed browser canvas size limits.
  const DESIRED_SCALE = 4;
  const MAX_DIM = 8000;
  const SCALE = Math.max(1, Math.min(DESIRED_SCALE, MAX_DIM / totalW, MAX_DIM / totalH));

  const canvas = document.createElement('canvas');
  canvas.width  = Math.round(totalW * SCALE);
  canvas.height = Math.round(totalH * SCALE);
  const ctx = canvas.getContext('2d');
  ctx.scale(SCALE, SCALE);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, totalW, totalH);

  // Title
  ctx.fillStyle = '#111';
  ctx.font = 'bold 13px monospace';
  ctx.fillText(document.querySelector('h1').textContent, 4, 22);

  async function svgToImage(svgEl) {{
    const clone = svgEl.cloneNode(true);
    const str = new XMLSerializer().serializeToString(clone);
    const blob = new Blob([str], {{type:'image/svg+xml'}});
    const url  = URL.createObjectURL(blob);
    return new Promise(resolve => {{
      const img = new Image();
      img.onload  = () => {{ URL.revokeObjectURL(url); resolve(img); }};
      img.onerror = () => {{ URL.revokeObjectURL(url); resolve(null); }};
      img.src = url;
    }});
  }}

  let y = 40;
  for (const row of rows) {{
    const svgEls   = row.querySelectorAll('svg');
    const trackSvg = svgEls[0] || null;
    const barSvg   = svgEls[1] || null;
    if (trackSvg) {{
      const img = await svgToImage(trackSvg);
      if (img) ctx.drawImage(img, 0, y, trackW, rowH);
    }}
    if (barSvg) {{
      const img = await svgToImage(barSvg);
      if (img) ctx.drawImage(img, trackW, y, barW, rowH);
    }}
    y += rowH + GAP;
  }}

  // ── Legend table ──────────────────────────────────────────────────────────
  y += 8;
  ctx.fillStyle = '#111';
  ctx.font = 'bold 10px monospace';
  ctx.fillText('Legend', 4, y);
  y += 14;

  const colW = Math.floor(totalW / LEG_COLS);
  const REF_PREFIX_PNG = DATA.ref_pfam.split('.')[0];
  legEntries.forEach((entry, i) => {{
    const col = i % LEG_COLS;
    const row = Math.floor(i / LEG_COLS);
    const lx = col * colW + 4;
    const ly = y + row * LEG_ROW_H;

    // color swatch
    ctx.fillStyle = entry.id.startsWith(REF_PREFIX_PNG) ? '#e63946' : entry.color;
    ctx.fillRect(lx, ly, SWATCH_W, SWATCH_H);
    ctx.strokeStyle = '#555';
    ctx.lineWidth = 0.5;
    ctx.strokeRect(lx, ly, SWATCH_W, SWATCH_H);

    // enrichment marker
    let marker = '';
    if (entry.sig_fdr)     marker = ' ★';
    else if (entry.sig_nominal) marker = ' ◆';

    ctx.fillStyle = entry.enriched ? '#c00' : '#222';
    ctx.font = '8px monospace';
    const label = entry.name + ' (' + entry.id + ')' + marker;
    ctx.fillText(label.slice(0, 38), lx + SWATCH_W + 3, ly + SWATCH_H - 1);
  }});

  const a = document.createElement('a');
  a.download = 'neighborhood_{ref_pfam}.png'.replace(/[.]/g,'_');
  a.href = canvas.toDataURL('image/png');
  a.click();
}}
</script>
</body>
</html>'''


def main():
    tab_file, pfam_id, window, ratio_cut = parse_args()
    data    = load_data(tab_file)
    js_data = build_js_data(data, pfam_id, window, ratio_cut)

    base_name = os.path.splitext(os.path.basename(tab_file))[0]
    out_file  = f"{base_name}_{pfam_id.replace('.', '_')}_w{window}_r{ratio_cut}.html"
    out_path  = os.path.join(os.path.dirname(os.path.abspath(tab_file)), out_file)

    with open(out_path, 'w') as f:
        f.write(render_html(js_data))
    print(f"Output         : {out_path}")


if __name__ == '__main__':
    main()
