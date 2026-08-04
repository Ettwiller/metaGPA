#!/usr/bin/env python3
"""
Recursive PFAM co-enrichment network.

Starting from a seed PFAM, finds all PFAMs significantly enriched (BH FDR q<0.05)
in high-ratio contigs, then recurses on each newly found PFAM until no new
significant associations are found. Outputs an interactive D3 force-directed network.

--hmmer_evalue_max caps hits by the domain's independent E-value (i-Evalue,
column 13 of the tab file): only hits with i-Evalue < this threshold are kept.

Usage:
    python pfam_coenrichment_network.py --hmmer_output <tab_file> --hmm_id <pfam_id> --window <window_bp> --ratio <ratio_cutoff> [--max_depth N] [--hmmer_evalue_max VALUE]

Example:
    python pfam_coenrichment_network.py --hmmer_output data.tab --hmm_id PF05014.22 --window 1000 --ratio 10 --max_depth 3 --hmmer_evalue_max 1e-5
"""

import re, json, sys, os, math
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def parse_args():
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--hmmer_output',     required=True, metavar='TAB_FILE',      help='HMMER-format PFAM annotation tab file')
    p.add_argument('--hmm_id',           required=True, metavar='PFAM_ID',       help='Seed PFAM domain ID (e.g. PF05014.22); version suffix ignored', type=lambda x: x.split('.')[0])
    p.add_argument('--window',           required=True, metavar='WINDOW_BP',     type=int,   help='Window in bp around seed domain')
    p.add_argument('--ratio',            required=True, metavar='RATIO_CUTOFF',  type=float, help='Enrichment ratio cutoff')
    p.add_argument('--max_depth',        default=3,     metavar='N',             type=int,   help='Recursion depth (default: 3)')
    p.add_argument('--hmmer_evalue_max', default=None,  metavar='VALUE',         type=float, help='Max i-Evalue for HMMER hits')
    a = p.parse_args()
    return a.hmmer_output, a.hmm_id, a.window, a.ratio, a.max_depth, a.hmmer_evalue_max


def load_data(tab_file, hmmer_evalue_max=None):
    if not os.path.exists(tab_file):
        print(f"Error: file not found: {tab_file}", file=sys.stderr)
        sys.exit(1)
    data = defaultdict(list)
    with open(tab_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 19:
                continue
            if hmmer_evalue_max is not None:
                try:
                    i_evalue = float(parts[12])
                except ValueError:
                    continue
                if i_evalue >= hmmer_evalue_max:
                    continue
            m = re.match(r'^(.+?_ratio_[\d.]+)-\d+[FR]$', parts[0])
            if not m:
                continue
            data[m.group(1)].append({
                'pfam_id':   parts[4],
                'pfam_name': parts[3],
                'start':     int(parts[17]),
                'end':       int(parts[18]),
            })
    return data


def get_ratio(base):
    m = re.search(r'ratio_([\d.]+)', base)
    return float(m.group(1)) if m else 0.0


def get_neighbors(anns, ref_prefix, window):
    ref = next((a for a in anns if a['pfam_id'].startswith(ref_prefix)), None)
    if not ref:
        return []
    return [a for a in anns if a['end'] >= ref['start'] - window and a['start'] <= ref['end'] + window]


def find_enriched(data, ref_prefix, window, ratio_cut):
    """Return list of (pfam_id, pfam_name, p_fdr) enriched FDR q<0.05 in high-ratio contigs."""
    pf_contigs = {b: anns for b, anns in data.items()
                  if any(a['pfam_id'].startswith(ref_prefix) for a in anns)}
    if not pf_contigs:
        return []

    pfam_name_map = {}
    rows = []
    for base, anns in pf_contigs.items():
        neighbors = get_neighbors(anns, ref_prefix, window)
        for a in neighbors:
            pfam_name_map[a['pfam_id']] = a['pfam_name']
        rows.append({'ratio': get_ratio(base), 'neighbors': neighbors})

    above = [r for r in rows if r['ratio'] >= ratio_cut]
    below = [r for r in rows if r['ratio'] <  ratio_cut]
    if not above:
        return []

    # Candidates: PFAMs present in at least one above-cutoff contig, excluding self
    candidates = sorted({a['pfam_id'] for r in above for a in r['neighbors']
                         if not a['pfam_id'].startswith(ref_prefix)})
    if not candidates:
        return []

    pvals = []
    for pid in candidates:
        ay = sum(1 for r in above if any(a['pfam_id'] == pid for a in r['neighbors']))
        by = sum(1 for r in below if any(a['pfam_id'] == pid for a in r['neighbors']))
        _, p = fisher_exact([[ay, len(above) - ay], [by, len(below) - by]], alternative='greater')
        pvals.append(p)

    reject, pvals_fdr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    return [(pid, pfam_name_map.get(pid, pid), float(p_fdr))
            for pid, p_fdr, rej in zip(candidates, pvals_fdr, reject) if rej]


def build_network(data, seed_pfam_id, window, ratio_cut, max_depth=3):
    seed_prefix = seed_pfam_id.split('.')[0]

    # Resolve full ID and name for seed
    pfam_name_map = {}
    pfam_full_id  = {}  # prefix -> first full id seen
    for anns in data.values():
        for a in anns:
            pfam_name_map[a['pfam_id']] = a['pfam_name']
            pfam_full_id.setdefault(a['pfam_id'].split('.')[0], a['pfam_id'])

    nodes = {}   # prefix -> {id, name, depth}
    edges = {}   # (src_prefix, tgt_prefix) -> min p_fdr (keep strongest)

    nodes[seed_prefix] = {
        'id':    pfam_full_id.get(seed_prefix, seed_pfam_id),
        'name':  pfam_name_map.get(pfam_full_id.get(seed_prefix, seed_pfam_id), seed_pfam_id),
        'depth': 0,
    }

    queue    = [seed_prefix]
    visited  = {seed_prefix}
    depth    = 0

    while queue:
        next_queue = []
        depth += 1
        for current_prefix in queue:
            enriched = find_enriched(data, current_prefix, window, ratio_cut)
            for pid, pname, p_fdr in enriched:
                prefix = pid.split('.')[0]
                # Add node
                if prefix not in nodes:
                    nodes[prefix] = {
                        'id':    pfam_full_id.get(prefix, pid),
                        'name':  pname,
                        'depth': depth,
                    }
                # Add / update edge (keep minimum p_fdr between any pair)
                key = tuple(sorted([current_prefix, prefix]))
                if key not in edges or p_fdr < edges[key]:
                    edges[key] = p_fdr
                # Queue for further recursion
                if prefix not in visited:
                    visited.add(prefix)
                    next_queue.append(prefix)

        n_new = len(next_queue)
        print(f"  depth {depth}: {n_new} new nodes, {len(nodes)} total nodes, {len(edges)} edges")
        if not n_new or depth >= max_depth:
            break
        queue = next_queue

    # Build JS-ready structures
    prefix_list = list(nodes.keys())
    node_list = []
    for prefix in prefix_list:
        n = nodes[prefix]
        node_list.append({
            'id':      prefix,
            'full_id': n['id'],
            'name':    n['name'],
            'depth':   n['depth'],
            'is_seed': prefix == seed_prefix,
        })

    edge_list = []
    for (src, tgt), p_fdr in edges.items():
        edge_list.append({
            'source': src,
            'target': tgt,
            'p_fdr':  p_fdr,
            'weight': round(-math.log10(max(p_fdr, 1e-300)), 3),
        })

    return node_list, edge_list


def load_d3():
    """Return inline D3 JS, downloading if needed."""
    import os, requests
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'd3.min.js')
    if not os.path.exists(cache):
        print("Downloading D3.js...", end=' ', flush=True)
        r = requests.get('https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js', timeout=15)
        r.raise_for_status()
        with open(cache, 'wb') as f:
            f.write(r.content)
        print("done.")
    with open(cache) as f:
        return f.read()


def render_html(node_list, edge_list, seed_pfam_id, window, ratio_cut):
    d3_js  = load_d3()
    js_data = json.dumps({
        'nodes': node_list,
        'edges': edge_list,
        'seed':  seed_pfam_id.split('.')[0],
        'seed_full': seed_pfam_id,
        'window': window,
        'ratio_cut': ratio_cut,
    })

    max_weight = max((e['weight'] for e in edge_list), default=1)
    n_nodes = len(node_list)
    n_edges = len(edge_list)

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>PFAM co-enrichment: {seed_pfam_id}</title>
<style>
  body   {{ margin:0; background:#fff; font-family:monospace; overflow:hidden; }}
  #info  {{ position:fixed; top:10px; left:10px; color:#222; font-size:11px; z-index:10;
            background:rgba(240,240,245,0.85); padding:8px 12px; border-radius:6px; line-height:1.6;
            border:1px solid #ccc; }}
  #tooltip {{
    position:fixed; background:rgba(30,30,30,0.92); color:#fff;
    padding:7px 11px; border-radius:5px; font-size:11px; pointer-events:none;
    display:none; max-width:280px; z-index:999; line-height:1.6;
  }}
  #legend {{ position:fixed; bottom:14px; left:14px; color:#333; font-size:10px;
             background:rgba(240,240,245,0.85); padding:8px 12px; border-radius:6px;
             border:1px solid #ccc; }}
  .legend-row {{ display:flex; align-items:center; gap:8px; margin:3px 0; }}
  .depth-dot  {{ width:10px; height:10px; border-radius:50%; flex-shrink:0; }}
  #dl-bar {{ position:fixed; top:10px; right:14px; display:flex; gap:6px; z-index:200; }}
  #dl-bar button {{
    padding:4px 10px; font-size:11px; font-family:monospace; cursor:pointer;
    border:1px solid #888; border-radius:4px; background:#f0f0f0;
  }}
  #dl-bar button:hover {{ background:#ddd; }}
  @media print {{
    #dl-bar, #info, #legend, #tooltip {{ display:none !important; }}
  }}
</style>
</head>
<body>
<div id="dl-bar">
  <button onclick="downloadPNG()">&#8659; PNG</button>
  <button onclick="window.print()">&#8659; PDF</button>
</div>
<div id="info">
  <b>{seed_pfam_id}</b> co-enrichment network &nbsp;|&nbsp; ±{window} bp &nbsp;|&nbsp; ratio &ge; {ratio_cut}<br>
  {n_nodes} nodes &nbsp;|&nbsp; {n_edges} edges &nbsp;|&nbsp; edge thickness ∝ −log₁₀(FDR q)<br>
  <span style="color:#888">scroll=zoom &nbsp;·&nbsp; drag=pan &nbsp;·&nbsp; drag node=move &nbsp;·&nbsp; click=EBI</span>
</div>
<div id="tooltip"></div>
<div id="legend"></div>
<svg id="svg" style="width:100vw;height:100vh;"></svg>

<script>{d3_js}</script>
<script>
const DATA = {js_data};
const MAX_WEIGHT = {max_weight};

const DEPTH_COLORS = ['#e63946','#1d6fa4','#e07b00','#2a7d5e',
                      '#7b3fbf','#b5830a','#1e7d45','#c2185b','#5c6bc0','#6d4c41'];
function depthColor(d) {{ return DEPTH_COLORS[Math.min(d, DEPTH_COLORS.length-1)]; }}

// Legend
const leg = document.getElementById('legend');
leg.innerHTML = '<b style="font-size:11px">Depth from seed</b>';
const depths = [...new Set(DATA.nodes.map(n => n.depth))].sort((a,b) => a-b);
depths.forEach(d => {{
  const row = document.createElement('div'); row.className = 'legend-row';
  const dot = document.createElement('div'); dot.className = 'depth-dot';
  dot.style.background = depthColor(d);
  const lbl = document.createElement('span');
  lbl.textContent = d === 0 ? '0 – seed' : String(d);
  row.appendChild(dot); row.appendChild(lbl); leg.appendChild(row);
}});
leg.innerHTML += `<div style="margin-top:4px;color:#888;font-size:9px">Edge thickness = −log₁₀(FDR q), max={max_weight:.1f}</div>`;

// Pre-compute degree
const degree = {{}};
DATA.nodes.forEach(n => degree[n.id] = 0);
DATA.edges.forEach(e => {{ degree[e.source] = (degree[e.source]||0)+1; degree[e.target] = (degree[e.target]||0)+1; }});
const maxDeg = Math.max(...Object.values(degree), 1);

const edgeWidthScale = d3.scaleLinear().domain([0, MAX_WEIGHT]).range([0.4, 6]);
const nodeRadiusScale = d3.scaleSqrt().domain([0, maxDeg]).range([4, 18]);

const W = window.innerWidth, H = window.innerHeight;
const svg = d3.select('#svg');
const g   = svg.append('g');

svg.call(d3.zoom().scaleExtent([0.05, 8]).on('zoom', e => g.attr('transform', e.transform)));

const simulation = d3.forceSimulation(DATA.nodes)
  .force('link',      d3.forceLink(DATA.edges).id(d => d.id)
                        .distance(d => 60 + (MAX_WEIGHT - d.weight) * 5).strength(0.6))
  .force('charge',    d3.forceManyBody().strength(-220))
  .force('center',    d3.forceCenter(W/2, H/2))
  .force('collision', d3.forceCollide().radius(d => nodeRadiusScale(degree[d.id]||0) + 3));

const link = g.append('g')
  .selectAll('line').data(DATA.edges).join('line')
  .attr('stroke', '#aab')
  .attr('stroke-width',   d => edgeWidthScale(d.weight))
  .attr('stroke-opacity', d => 0.3 + 0.55 * (d.weight / MAX_WEIGHT));

const node = g.append('g')
  .selectAll('circle').data(DATA.nodes).join('circle')
  .attr('r',            d => nodeRadiusScale(degree[d.id]||0))
  .attr('fill',         d => depthColor(d.depth))
  .attr('stroke',       d => d.is_seed ? '#333' : 'none')
  .attr('stroke-width', d => d.is_seed ? 2 : 0)
  .style('cursor', 'pointer')
  .call(d3.drag()
    .on('start', (ev, d) => {{ if (!ev.active) simulation.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; }})
    .on('drag',  (ev, d) => {{ d.fx=ev.x; d.fy=ev.y; }})
    .on('end',   (ev, d) => {{ if (!ev.active) simulation.alphaTarget(0); d.fx=null; d.fy=null; }}));

const label = g.append('g')
  .selectAll('text').data(DATA.nodes).join('text')
  .text(d => d.name)
  .attr('font-size',      d => d.is_seed ? '11px' : '8px')
  .attr('font-weight',    d => d.is_seed ? 'bold' : 'normal')
  .attr('fill',           d => d.is_seed ? '#111' : '#444')
  .attr('pointer-events', 'none')
  .attr('dx', d => nodeRadiusScale(degree[d.id]||0) + 2)
  .attr('dy', '0.35em');

const tooltip = document.getElementById('tooltip');
document.addEventListener('mousemove', e => {{
  tooltip.style.left = (e.clientX + 14) + 'px';
  tooltip.style.top  = (e.clientY - 10) + 'px';
}});

node
  .on('mouseenter', (ev, d) => {{
    const deg = degree[d.id] || 0;
    const nbrs = DATA.edges
      .filter(e => e.source.id===d.id || e.target.id===d.id)
      .sort((a,b) => a.p_fdr - b.p_fdr)
      .map(e => {{
        const o = e.source.id===d.id ? e.target : e.source;
        return `&bull; ${{o.name}} (${{o.full_id}}) q=${{e.p_fdr.toExponential(1)}}`;
      }}).join('<br>');
    tooltip.innerHTML =
      `<b>${{d.name}}</b> (${{d.full_id}})<br>`+
      `depth: ${{d.depth}} &nbsp;|&nbsp; connections: ${{deg}}<br>`+
      `<span style="color:#bbb">connected to:</span><br>${{nbrs||'none'}}`;
    tooltip.style.display = 'block';
    link
      .attr('stroke',         e => (e.source.id===d.id||e.target.id===d.id) ? '#222' : '#dde')
      .attr('stroke-opacity', e => (e.source.id===d.id||e.target.id===d.id) ? 1.0 : 0.12);
  }})
  .on('mouseleave', () => {{
    tooltip.style.display = 'none';
    link
      .attr('stroke', '#aab')
      .attr('stroke-opacity', d => 0.3 + 0.55*(d.weight/MAX_WEIGHT));
  }})
  .on('click', (ev, d) => {{
    window.open('https://www.ebi.ac.uk/interpro/entry/pfam/'+d.full_id.split('.')[0]+'/', '_blank');
  }});

simulation.on('tick', () => {{
  link.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
  node .attr('cx', d => d.x).attr('cy', d => d.y);
  label.attr('x',  d => d.x).attr('y',  d => d.y);
}});

function downloadPNG() {{
  const svgEl = document.getElementById('svg');
  const W = svgEl.clientWidth, H = svgEl.clientHeight;
  // Render at a high pixel density for publication-quality output, capped so
  // we never exceed browser canvas size limits.
  const DESIRED_SCALE = 4;
  const MAX_DIM = 8000;
  const SCALE = Math.max(1, Math.min(DESIRED_SCALE, MAX_DIM / W, MAX_DIM / H));
  const serializer = new XMLSerializer();
  const svgStr = serializer.serializeToString(svgEl);
  const blob = new Blob([svgStr], {{type:'image/svg+xml'}});
  const url  = URL.createObjectURL(blob);
  const img  = new Image();
  img.onload = () => {{
    const canvas = document.createElement('canvas');
    canvas.width = Math.round(W * SCALE); canvas.height = Math.round(H * SCALE);
    const ctx = canvas.getContext('2d');
    ctx.scale(SCALE, SCALE);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, W, H);
    ctx.drawImage(img, 0, 0, W, H);
    URL.revokeObjectURL(url);
    const a = document.createElement('a');
    a.download = 'network_{seed_pfam_id}.png'.replace(/[.]/g,'_');
    a.href = canvas.toDataURL('image/png');
    a.click();
  }};
  img.src = url;
}}
</script>
</body>
</html>'''


def main():
    tab_file, pfam_id, window, ratio_cut, max_depth, hmmer_evalue_max = parse_args()
    data = load_data(tab_file, hmmer_evalue_max)

    print(f"Seed PFAM      : {pfam_id}")
    print(f"Window         : ±{window} bp")
    print(f"Ratio cutoff   : {ratio_cut}")
    print(f"Max depth      : {max_depth}")
    if hmmer_evalue_max is not None:
        print(f"E-value cutoff : i-Evalue < {hmmer_evalue_max}")
    print("Building co-enrichment network...")

    node_list, edge_list = build_network(data, pfam_id, window, ratio_cut, max_depth)

    print(f"Network        : {len(node_list)} nodes, {len(edge_list)} edges")

    base_name = os.path.splitext(os.path.basename(tab_file))[0]
    out_file  = f"{base_name}_{pfam_id.replace('.','_')}_w{window}_r{ratio_cut}_d{max_depth}_network.html"
    out_path  = os.path.join(os.path.dirname(os.path.abspath(tab_file)), out_file)

    with open(out_path, 'w') as f:
        f.write(render_html(node_list, edge_list, pfam_id, window, ratio_cut))
    print(f"Output         : {out_path}")


if __name__ == '__main__':
    main()
