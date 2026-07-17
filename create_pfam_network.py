#!/usr/bin/env python3
"""
Linked PFAM visualization: network + per-node neighborhood plots.

Generates an output folder containing:
  - network.html          : co-enrichment network (click node → neighborhood)
  - neighborhood_PF*.html : neighborhood plot for every node in the network

Usage:
    python pfam_linked_viz.py <tab_file> <pfam_id> <window_bp> <ratio_cutoff> [max_depth]

Example:
    python pfam_linked_viz.py data.tab PF05014.22 1000 10 3
"""

import sys, os, re, json, importlib.util

# ── Load sibling modules ───────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))

def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, f'{name}.py'))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

nb  = _load('plot_pfam_neighborhood')       # neighborhood script
net = _load('pfam_coenrichment_network')    # network script


# ── Args ───────────────────────────────────────────────────────────────────────
def parse_args():
    if len(sys.argv) not in (5, 6):
        print(__doc__)
        sys.exit(1)
    max_depth = int(sys.argv[5]) if len(sys.argv) == 6 else 3
    return sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), max_depth


# ── Network HTML (click → neighborhood) ───────────────────────────────────────
def render_network_html(node_list, edge_list, seed_pfam_id, window, ratio_cut, d3_js):
    js_data    = json.dumps({'nodes': node_list, 'edges': edge_list,
                             'seed': seed_pfam_id.split('.')[0], 'seed_full': seed_pfam_id,
                             'window': window, 'ratio_cut': ratio_cut})
    max_weight = max((e['weight'] for e in edge_list), default=1)
    n_nodes, n_edges = len(node_list), len(edge_list)

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>PFAM co-enrichment: {seed_pfam_id}</title>
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
  <b>{seed_pfam_id}</b> co-enrichment network &nbsp;|&nbsp; ±{window} bp &nbsp;|&nbsp; ratio &ge; {ratio_cut}<br>
  {n_nodes} nodes &nbsp;|&nbsp; {n_edges} edges &nbsp;|&nbsp; edge thickness ∝ −log₁₀(FDR q)<br>
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
const DATA = {js_data};
const MAX_WEIGHT = {max_weight};

const DEPTH_COLORS = ['#e63946','#1d6fa4','#e07b00','#2a7d5e',
                      '#7b3fbf','#b5830a','#1e7d45','#c2185b','#5c6bc0','#6d4c41'];
function depthColor(d) {{ return DEPTH_COLORS[Math.min(d, DEPTH_COLORS.length-1)]; }}

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

const degree = {{}};
DATA.nodes.forEach(n => degree[n.id] = 0);
DATA.edges.forEach(e => {{ degree[e.source]=(degree[e.source]||0)+1; degree[e.target]=(degree[e.target]||0)+1; }});
const maxDeg = Math.max(...Object.values(degree), 1);

const edgeWidthScale  = d3.scaleLinear().domain([0, MAX_WEIGHT]).range([0.4, 6]);
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
        return `&bull; ${{o.name}} q=${{e.p_fdr.toExponential(1)}}`;
      }}).join('<br>');
    tooltip.innerHTML =
      `<b>${{d.name}}</b> (${{d.full_id}})<br>`+
      `depth: ${{d.depth}} &nbsp;|&nbsp; connections: ${{deg}}<br>`+
      `<span style="color:#aaa">click → neighborhood &nbsp;·&nbsp; shift+click → EBI</span><br>`+
      `<span style="color:#bbb">connected to:</span><br>${{nbrs||'none'}}`;
    tooltip.style.display = 'block';
    link
      .attr('stroke',         e => (e.source.id===d.id||e.target.id===d.id) ? '#222' : '#dde')
      .attr('stroke-opacity', e => (e.source.id===d.id||e.target.id===d.id) ? 1.0  : 0.12);
  }})
  .on('mouseleave', () => {{
    tooltip.style.display = 'none';
    link.attr('stroke', '#aab')
        .attr('stroke-opacity', d => 0.3 + 0.55*(d.weight/MAX_WEIGHT));
  }})
  .on('click', (ev, d) => {{
    if (ev.shiftKey) {{
      // Shift+click → EBI
      window.open('https://www.ebi.ac.uk/interpro/entry/pfam/'+d.full_id.split('.')[0]+'/', '_blank');
    }} else {{
      // Regular click → neighborhood plot
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


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    tab_file, pfam_id, window, ratio_cut, max_depth = parse_args()

    if not os.path.exists(tab_file):
        print(f"Error: file not found: {tab_file}", file=sys.stderr)
        sys.exit(1)

    # Output directory
    base     = os.path.splitext(os.path.basename(tab_file))[0]
    out_dir  = os.path.join(os.path.dirname(os.path.abspath(tab_file)),
                            f"{base}_{pfam_id.replace('.','_')}_w{window}_r{ratio_cut}_d{max_depth}_linked")
    os.makedirs(out_dir, exist_ok=True)
    print(f"Output dir     : {out_dir}")

    # Load data once, shared by both scripts
    print("Loading data...")
    data = net.load_data(tab_file)

    # ── Build network ──────────────────────────────────────────────────────────
    print(f"Building co-enrichment network (depth={max_depth})...")
    node_list, edge_list = net.build_network(data, pfam_id, window, ratio_cut, max_depth)
    print(f"Network        : {len(node_list)} nodes, {len(edge_list)} edges")

    d3_js = net.load_d3()
    network_html = render_network_html(node_list, edge_list, pfam_id, window, ratio_cut, d3_js)
    network_path = os.path.join(out_dir, 'network.html')
    with open(network_path, 'w') as f:
        f.write(network_html)
    print(f"Network HTML   : {network_path}")

    # ── Build one neighborhood HTML per network node ───────────────────────────
    print(f"Generating {len(node_list)} neighborhood plots...")
    for i, node in enumerate(node_list):
        pfam_full = node['full_id']   # e.g. PF05014.22
        pfam_prefix = node['id']      # e.g. PF05014

        try:
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                js_data = nb.build_js_data(data, pfam_full, window, ratio_cut)
            html = nb.render_html(js_data)
        except SystemExit:
            continue

        out_path = os.path.join(out_dir, f'neighborhood_{pfam_prefix}.html')
        with open(out_path, 'w') as f:
            f.write(html)

        if (i+1) % 20 == 0 or (i+1) == len(node_list):
            print(f"  {i+1}/{len(node_list)} done")

    print(f"\nDone. Open in browser:\n  {network_path}")


if __name__ == '__main__':
    main()
