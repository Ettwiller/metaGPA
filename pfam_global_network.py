#!/usr/bin/env python3
"""
Global PFAM co-enrichment network (no seed required).

Step 1 – Find all PFAMs significantly enriched in high-ratio contigs
         vs low-ratio contigs (Fisher's exact test, BH/FDR correction).
Step 2 – Among enriched PFAMs, test pairwise co-occurrence on high-ratio
         contigs (Fisher's exact, BH/FDR). Build network from significant pairs.
Step 3 – Render an interactive D3 force-directed HTML network.

Clicking a node opens its neighborhood plot (like create_pfam_network.py);
shift+click opens the EBI InterPro page instead.

Usage:
    python pfam_global_network.py --hmmer_output <tab_file> --ratio <ratio_cutoff> [--fdr <fdr_threshold>] [--min_contigs N] [--window <window_bp>]

Arguments:
    tab_file        HMMER-format PFAM annotation file
    ratio_cutoff    Contigs with ratio >= this are "high"; others are "low"
    fdr_threshold   BH/FDR q-value cutoff (default: 0.05)
    min_contigs     Min number of high-ratio contigs a PFAM must appear in
                    to be tested (default: 3)
    window_bp       +/- bp window used for the per-node neighborhood plots
                    (default: 1000)
    --hmmer_evalue_max VALUE
                    Caps hits by the domain's independent E-value (i-Evalue,
                    column 13 of the tab file): only hits with i-Evalue <
                    this threshold are kept.

Example:
    python pfam_global_network.py --hmmer_output all_contigs.nd_pfam_with_enrichment.tab --ratio 10
    python pfam_global_network.py --hmmer_output all_contigs.nd_pfam_with_enrichment.tab --ratio 10 --fdr 0.05 --min_contigs 5 --window 1000
    python pfam_global_network.py --hmmer_output all_contigs.nd_pfam_with_enrichment.tab --ratio 10 --hmmer_evalue_max 1e-5
"""

import re, json, sys, os, math, urllib.request, io, contextlib, importlib.util
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))

def _load_sibling(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, f'{name}.py'))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

nb = _load_sibling('plot_pfam_neighborhood')


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args():
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--hmmer_output',     required=True, metavar='TAB_FILE',      help='HMMER-format PFAM annotation tab file')
    p.add_argument('--ratio',            required=True, metavar='RATIO_CUTOFF',  type=float, help='Enrichment ratio cutoff')
    p.add_argument('--fdr',              default=0.05,  metavar='FDR_THRESHOLD', type=float, help='BH/FDR q-value cutoff (default: 0.05)')
    p.add_argument('--min_contigs',      default=3,     metavar='N',             type=int,   help='Min high-ratio contigs per PFAM to test (default: 3)')
    p.add_argument('--window',           default=1000,  metavar='WINDOW_BP',     type=int,   help='Window in bp for neighborhood plots (default: 1000)')
    p.add_argument('--hmmer_evalue_max', default=None,  metavar='VALUE',         type=float, help='Max i-Evalue for HMMER hits')
    a = p.parse_args()
    return a.hmmer_output, a.ratio, a.fdr, a.min_contigs, a.window, a.hmmer_evalue_max


# ── Data loading ───────────────────────────────────────────────────────────────
def load_data(tab_file, hmmer_evalue_max=None):
    """Return dict: contig_base -> list of {pfam_id, pfam_name}."""
    if not os.path.exists(tab_file):
        print(f"Error: file not found: {tab_file}", file=sys.stderr)
        sys.exit(1)
    data = defaultdict(list)
    with open(tab_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 19:
                continue
            m = re.match(r'^(.+?_ratio_[\d.]+)-\d+[FR]$', parts[0])
            if not m:
                continue
            pfam_id   = parts[4]
            pfam_name = parts[3]
            if not pfam_id.startswith('PF'):
                continue
            if hmmer_evalue_max is not None:
                try:
                    i_evalue = float(parts[12])
                except ValueError:
                    continue
                if i_evalue >= hmmer_evalue_max:
                    continue
            base = m.group(1)
            # deduplicate per contig (same PFAM can appear multiple times as different genes)
            data[base].append({'pfam_id': pfam_id, 'pfam_name': pfam_name})
    return data


def get_ratio(base):
    m = re.search(r'ratio_([\d.]+)', base)
    return float(m.group(1)) if m else 0.0


# ── Step 1: global PFAM enrichment ────────────────────────────────────────────
def find_enriched_pfams(data, ratio_cut, fdr_threshold, min_contigs):
    """
    Returns:
        enriched  : dict pfam_id -> {pfam_name, p_raw, p_fdr, n_high, n_low}
        above     : list of contig bases with ratio >= ratio_cut
        below     : list of contig bases with ratio < ratio_cut
        pfam_names: dict pfam_id -> pfam_name
    """
    above = [b for b in data if get_ratio(b) >= ratio_cut]
    below = [b for b in data if get_ratio(b) <  ratio_cut]
    n_above, n_below = len(above), len(below)
    print(f"Contigs above cutoff : {n_above}")
    print(f"Contigs below cutoff : {n_below}")

    # Build presence sets: pfam_id -> set of contig bases
    pfam_names = {}
    above_set  = defaultdict(set)
    below_set  = defaultdict(set)

    for b in above:
        for ann in data[b]:
            pid = ann['pfam_id'].split('.')[0]  # normalise to base ID
            above_set[pid].add(b)
            pfam_names[pid] = ann['pfam_name']

    for b in below:
        for ann in data[b]:
            pid = ann['pfam_id'].split('.')[0]
            below_set[pid].add(b)
            pfam_names.setdefault(pid, ann['pfam_name'])

    # Filter: PFAM must appear in at least min_contigs high-ratio contigs
    candidates = [pid for pid, s in above_set.items() if len(s) >= min_contigs]
    print(f"PFAMs tested         : {len(candidates)}  (≥{min_contigs} high-ratio contigs)")

    pvals_raw = []
    for pid in candidates:
        a_yes = len(above_set[pid])
        b_yes = len(below_set.get(pid, set()))
        a_no  = n_above - a_yes
        b_no  = n_below - b_yes
        _, p  = fisher_exact([[a_yes, a_no], [b_yes, b_no]], alternative='greater')
        pvals_raw.append(p)

    reject, pvals_fdr, _, _ = multipletests(pvals_raw, alpha=fdr_threshold, method='fdr_bh')

    enriched = {}
    for pid, p_raw, p_fdr, rej in zip(candidates, pvals_raw, pvals_fdr, reject):
        if rej:
            enriched[pid] = {
                'pfam_name': pfam_names[pid],
                'p_raw': p_raw,
                'p_fdr': p_fdr,
                'n_high': len(above_set[pid]),
                'n_low':  len(below_set.get(pid, set())),
            }

    print(f"Enriched PFAMs (FDR<{fdr_threshold}) : {len(enriched)}")
    return enriched, above, below, pfam_names, above_set


# ── Step 2: pairwise co-occurrence among enriched PFAMs ───────────────────────
def build_cooccurrence_network(enriched, above, above_set, fdr_threshold):
    """
    For each pair of enriched PFAMs, test whether they co-occur in high-ratio
    contigs more than expected by chance (Fisher's exact, two-sided).
    Returns node_list, edge_list.
    """
    n_above   = len(above)
    above_set_frozen = {pid: frozenset(s) for pid, s in above_set.items() if pid in enriched}
    pids      = sorted(enriched.keys())
    n         = len(pids)

    if n == 0:
        return [], []

    # Pairwise tests
    pairs     = [(pids[i], pids[j]) for i in range(n) for j in range(i+1, n)]
    pvals_raw = []
    for pid_a, pid_b in pairs:
        set_a = above_set_frozen[pid_a]
        set_b = above_set_frozen[pid_b]
        both  = len(set_a & set_b)
        only_a = len(set_a) - both
        only_b = len(set_b) - both
        neither = n_above - both - only_a - only_b
        if neither < 0:
            neither = 0
        _, p = fisher_exact([[both, only_a], [only_b, neither]], alternative='greater')
        pvals_raw.append(p)

    print(f"Pairwise co-occurrence tests : {len(pairs)}")

    if not pvals_raw:
        return [], []

    reject, pvals_fdr, _, _ = multipletests(pvals_raw, alpha=fdr_threshold, method='fdr_bh')

    sig_edges = []
    for (pid_a, pid_b), p_raw, p_fdr, rej in zip(pairs, pvals_raw, pvals_fdr, reject):
        if rej:
            weight = -math.log10(max(p_fdr, 1e-300))
            sig_edges.append({
                'source': pid_a, 'target': pid_b,
                'p_raw': round(p_raw, 8), 'p_fdr': round(p_fdr, 8),
                'weight': round(weight, 4),
            })

    # Only keep nodes that appear in at least one edge
    connected = {e['source'] for e in sig_edges} | {e['target'] for e in sig_edges}
    # If no edges, include all enriched nodes anyway
    if not connected:
        connected = set(pids)

    node_list = [
        {
            'id':        pid,
            'name':      enriched[pid]['pfam_name'],
            'full_id':   pid,
            'p_fdr':     round(enriched[pid]['p_fdr'], 8),
            'n_high':    enriched[pid]['n_high'],
            'n_low':     enriched[pid]['n_low'],
        }
        for pid in pids if pid in connected
    ]

    print(f"Network nodes (connected) : {len(node_list)}")
    print(f"Network edges (FDR<{fdr_threshold})  : {len(sig_edges)}")
    return node_list, sig_edges


# ── D3 helper ─────────────────────────────────────────────────────────────────
def load_d3():
    here    = os.path.dirname(os.path.abspath(__file__))
    d3_path = os.path.join(here, 'd3.min.js')
    if os.path.exists(d3_path):
        with open(d3_path) as f:
            return f.read()
    url = 'https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js'
    print(f"Downloading D3 from {url} ...")
    urllib.request.urlretrieve(url, d3_path)
    with open(d3_path) as f:
        return f.read()


# ── HTML render ───────────────────────────────────────────────────────────────
def render_html(node_list, edge_list, ratio_cut, fdr_threshold, d3_js, tab_file):
    js_data    = json.dumps({'nodes': node_list, 'edges': edge_list,
                             'ratio_cut': ratio_cut, 'fdr': fdr_threshold})
    max_weight = max((e['weight'] for e in edge_list), default=1)
    n_nodes, n_edges = len(node_list), len(edge_list)
    base_name  = os.path.basename(tab_file)

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>PFAM global co-enrichment network (ratio &ge; {ratio_cut})</title>
<style>
  body   {{ margin:0; background:#fff; font-family:monospace; overflow:hidden; }}
  #info  {{ position:fixed; top:10px; left:10px; color:#222; font-size:11px; z-index:10;
            background:rgba(240,240,245,0.9); padding:8px 12px; border-radius:6px;
            line-height:1.6; border:1px solid #ccc; }}
  #tooltip {{
    position:fixed; background:rgba(30,30,30,0.92); color:#fff;
    padding:7px 11px; border-radius:5px; font-size:11px; pointer-events:none;
    display:none; max-width:300px; z-index:999; line-height:1.6;
  }}
  #legend {{ position:fixed; bottom:14px; left:14px; color:#333; font-size:10px;
             background:rgba(240,240,245,0.9); padding:8px 12px; border-radius:6px;
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
  <b>Global PFAM co-enrichment network</b> &nbsp;|&nbsp; ratio &ge; {ratio_cut}<br>
  {n_nodes} nodes &nbsp;|&nbsp; {n_edges} edges &nbsp;|&nbsp; edge thickness ∝ −log₁₀(co-occurrence FDR q)<br>
  <span style="color:#888">
    scroll=zoom &nbsp;·&nbsp; drag=pan &nbsp;·&nbsp; drag node=move<br>
    <b>click</b>=open neighborhood &nbsp;·&nbsp; <b>shift+click</b>=EBI page
  </span>
</div>
<div id="tooltip"></div>
<div id="legend"></div>
<svg id="svg" style="width:100vw;height:100vh;"></svg>

<script>{d3_js}</script>
<script>
const DATA       = {js_data};
const MAX_WEIGHT = {max_weight};

// Node color: blue (low enrichment) → red (high enrichment) based on -log10(p_fdr)
const maxNodeScore = Math.max(...DATA.nodes.map(n => -Math.log10(Math.max(n.p_fdr, 1e-300))), 1);
const colorScale   = d3.scaleSequential(d3.interpolateYlOrRd)
                       .domain([0, maxNodeScore]);
function nodeColor(n) {{ return colorScale(-Math.log10(Math.max(n.p_fdr, 1e-300))); }}

const degree = {{}};
DATA.nodes.forEach(n => degree[n.id] = 0);
DATA.edges.forEach(e => {{ degree[e.source]=(degree[e.source]||0)+1; degree[e.target]=(degree[e.target]||0)+1; }});
const maxDeg = Math.max(...Object.values(degree), 1);

const edgeWidthScale  = d3.scaleLinear().domain([0, MAX_WEIGHT]).range([0.3, 5]);
const nodeRadiusScale = d3.scaleSqrt().domain([0, maxDeg]).range([4, 18]);

const W = window.innerWidth, H = window.innerHeight;
const svg = d3.select('#svg');
const g   = svg.append('g');

// Legend: binned enrichment-strength swatches (same visual language as depth legend)
const leg = document.getElementById('legend');
leg.innerHTML = '<b style="font-size:11px">Enrichment strength (&minus;log&#8321;&#8320; FDR q)</b>';
const N_BINS = 5;
for (let i = N_BINS - 1; i >= 0; i--) {{
  const score = (i / (N_BINS - 1)) * maxNodeScore;
  const row = document.createElement('div'); row.className = 'legend-row';
  const dot = document.createElement('div'); dot.className = 'depth-dot';
  dot.style.background = colorScale(score);
  const lbl = document.createElement('span');
  lbl.textContent = i === N_BINS - 1 ? `high (${{score.toFixed(1)}}+)`
                   : i === 0          ? 'low (~0)'
                   : score.toFixed(1);
  row.appendChild(dot); row.appendChild(lbl); leg.appendChild(row);
}}
leg.innerHTML += `<div style="margin-top:4px;color:#888;font-size:9px">Node size = degree (co-occurrence partners)</div>`;

svg.call(d3.zoom().scaleExtent([0.02, 10]).on('zoom', e => g.attr('transform', e.transform)));

const simulation = d3.forceSimulation(DATA.nodes)
  .force('link',      d3.forceLink(DATA.edges).id(d => d.id)
                        .distance(d => 50 + (MAX_WEIGHT - d.weight) * 4).strength(0.5))
  .force('charge',    d3.forceManyBody().strength(-180))
  .force('center',    d3.forceCenter(W/2, H/2))
  .force('collision', d3.forceCollide().radius(d => nodeRadiusScale(degree[d.id]||0) + 2));

const link = g.append('g')
  .selectAll('line').data(DATA.edges).join('line')
  .attr('stroke', '#bbc')
  .attr('stroke-width',   d => edgeWidthScale(d.weight))
  .attr('stroke-opacity', d => 0.25 + 0.60 * (d.weight / MAX_WEIGHT));

const node = g.append('g')
  .selectAll('circle').data(DATA.nodes).join('circle')
  .attr('r',    d => nodeRadiusScale(degree[d.id]||0))
  .attr('fill', d => nodeColor(d))
  .attr('stroke', '#555').attr('stroke-width', 0.5)
  .style('cursor', 'pointer')
  .call(d3.drag()
    .on('start', (ev, d) => {{ if (!ev.active) simulation.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; }})
    .on('drag',  (ev, d) => {{ d.fx=ev.x; d.fy=ev.y; }})
    .on('end',   (ev, d) => {{ if (!ev.active) simulation.alphaTarget(0); d.fx=null; d.fy=null; }}));

const label = g.append('g')
  .selectAll('text').data(DATA.nodes).join('text')
  .text(d => d.name)
  .attr('font-size', '8px')
  .attr('fill', '#333')
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
      .sort((a,b) => a.p_fdr - b.p_fdr).slice(0, 8)
      .map(e => {{
        const o = e.source.id===d.id ? e.target : e.source;
        return `&bull; ${{o.name}} q=${{e.p_fdr.toExponential(1)}}`;
      }}).join('<br>');
    tooltip.innerHTML =
      `<b>${{d.name}}</b> (${{d.id}})<br>`+
      `enrichment FDR q=${{d.p_fdr.toExponential(2)}}<br>`+
      `in high-ratio: ${{d.n_high}} contigs &nbsp;|&nbsp; low: ${{d.n_low}}<br>`+
      `connections: ${{deg}}<br>`+
      `<span style="color:#aaa">click → neighborhood &nbsp;·&nbsp; shift+click → EBI</span><br>`+
      `<span style="color:#bbb">top co-occurring:</span><br>${{nbrs||'none'}}`;
    tooltip.style.display = 'block';
    link
      .attr('stroke',         e => (e.source.id===d.id||e.target.id===d.id) ? '#222' : '#dde')
      .attr('stroke-opacity', e => (e.source.id===d.id||e.target.id===d.id) ? 1.0  : 0.08);
  }})
  .on('mouseleave', () => {{
    tooltip.style.display = 'none';
    link.attr('stroke', '#bbc')
        .attr('stroke-opacity', d => 0.25 + 0.60*(d.weight/MAX_WEIGHT));
  }})
  .on('click', (ev, d) => {{
    if (ev.shiftKey) {{
      window.open('https://www.ebi.ac.uk/interpro/entry/pfam/' + d.id + '/', '_blank');
    }} else {{
      window.open('neighborhood_' + d.id + '.html', '_blank');
    }}
  }});

simulation.on('tick', () => {{
  link.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
  node .attr('cx', d => d.x).attr('cy', d => d.y);
  label.attr('x',  d => d.x).attr('y',  d => d.y);
}});

// ── PNG download ──────────────────────────────────────────────────────────────
function downloadPNG() {{
  const svgEl = document.getElementById('svg');
  const W = svgEl.clientWidth, H = svgEl.clientHeight;
  const str  = new XMLSerializer().serializeToString(svgEl);
  const blob = new Blob([str], {{type:'image/svg+xml'}});
  const url  = URL.createObjectURL(blob);
  const img  = new Image();
  img.onload = () => {{
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, W, H);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    const a = document.createElement('a');
    a.download = 'global_network_r{ratio_cut}.png';
    a.href = canvas.toDataURL('image/png');
    a.click();
  }};
  img.src = url;
}}
</script>
</body>
</html>'''


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    tab_file, ratio_cut, fdr_threshold, min_contigs, window_bp, hmmer_evalue_max = parse_args()

    if hmmer_evalue_max is not None:
        print(f"E-value cutoff : i-Evalue < {hmmer_evalue_max}")

    print(f"Loading data from {tab_file} ...")
    data = load_data(tab_file, hmmer_evalue_max)
    print(f"Total contigs loaded : {len(data)}")

    print(f"\nStep 1 – Finding enriched PFAMs (ratio >= {ratio_cut}, FDR < {fdr_threshold}) ...")
    enriched, above, below, pfam_names, above_set = find_enriched_pfams(
        data, ratio_cut, fdr_threshold, min_contigs)

    if not enriched:
        print("No enriched PFAMs found. Try lowering the ratio_cutoff or fdr_threshold.")
        sys.exit(0)

    print(f"\nStep 2 – Building pairwise co-occurrence network ...")
    node_list, edge_list = build_cooccurrence_network(enriched, above, above_set, fdr_threshold)

    d3_js = load_d3()

    base    = os.path.splitext(os.path.basename(tab_file))[0]
    out_dir = os.path.join(os.path.dirname(os.path.abspath(tab_file)),
                            f"{base}_global_r{ratio_cut}_fdr{fdr_threshold}_linked")
    os.makedirs(out_dir, exist_ok=True)

    network_html = render_html(node_list, edge_list, ratio_cut, fdr_threshold, d3_js, tab_file)
    network_path = os.path.join(out_dir, 'network.html')
    with open(network_path, 'w') as f:
        f.write(network_html)
    print(f"\nNetwork HTML : {network_path}")

    # ── Per-node neighborhood plots ──────────────────────────────────────────
    print(f"\nStep 3 – Generating {len(node_list)} neighborhood plots (window=±{window_bp}bp) ...")
    nb_data = nb.load_data(tab_file, hmmer_evalue_max)
    n_written = 0
    for i, node in enumerate(node_list):
        pfam_id = node['id']  # e.g. PF10417 (prefix match used internally by nb)
        try:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                js_data = nb.build_js_data(nb_data, pfam_id, window_bp, ratio_cut)
            html = nb.render_html(js_data)
        except SystemExit:
            continue

        out_path = os.path.join(out_dir, f'neighborhood_{pfam_id}.html')
        with open(out_path, 'w') as f:
            f.write(html)
        n_written += 1

        if (i + 1) % 20 == 0 or (i + 1) == len(node_list):
            print(f"  {i+1}/{len(node_list)} done")

    print(f"Neighborhood plots written : {n_written}/{len(node_list)}")
    print(f"\nDone. Open in browser:\n  {network_path}")


if __name__ == '__main__':
    main()
