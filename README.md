# metaGPA

Tools for exploring PFAM domain co-enrichment in ratio-labeled contig annotations
(e.g. hmmscan/hmmsearch `domtblout` output from a metaGPA-style enrichment screen).

## Scripts

### `create_pfam_network.py` (main entry point)

Builds an interactive co-enrichment network around a seed PFAM domain and, for
every node in that network, a linked neighborhood plot. Output is a folder of
HTML files — click a node in `network.html` to open its neighborhood plot.

```
python create_pfam_network.py <tab_file> <pfam_id> <window_bp> <ratio_cutoff> [max_depth]
```

Example:

```
python create_pfam_network.py data.tab PF05014.22 1000 10 3
```

- `tab_file` — hmmscan/hmmsearch-style domain table (whitespace-delimited).
- `pfam_id` — seed PFAM accession, e.g. `PF05014.22`.
- `window_bp` — ± window (bp) around the seed domain to look for neighboring domains.
- `ratio_cutoff` — contigs with ratio above this are "high-ratio"; used for enrichment testing.
- `max_depth` — optional, default `3`. How many recursive hops from the seed to follow when
  expanding the network.

Output is written to a folder next to `tab_file` named
`<tab_basename>_<pfam_id>_w<window>_r<ratio>_d<depth>_linked/`, containing:

- `network.html` — the co-enrichment network.
- `neighborhood_PF*.html` — one neighborhood plot per network node.

### `pfam_coenrichment_network.py`

Standalone version of the network builder (used internally by `create_pfam_network.py`).
Can also be run directly to produce just `network.html` (its node clicks link to EBI
instead of local neighborhood pages):

```
python pfam_coenrichment_network.py <tab_file> <pfam_id> <window_bp> <ratio_cutoff> [max_depth]
```

### `plot_pfam_neighborhood.py`

Standalone version of the neighborhood plot for a single PFAM domain:

```
python plot_pfam_neighborhood.py <tab_file> <pfam_id> <window_bp> <ratio_cutoff>
```

## How it works

Contigs are matched by name against `^(.+?_ratio_[\d.]+)-\d+[FR]$`, so each contig's
ratio is parsed straight out of its ID. For a given PFAM domain, contigs are split into
"high-ratio" (`ratio > ratio_cutoff`) and "low-ratio" groups. For every other PFAM found
within `window_bp` of the reference domain, a Fisher's exact test (BH-FDR corrected,
q < 0.05) checks whether it's enriched in the high-ratio group. `create_pfam_network.py`
recurses this process outward from the seed domain up to `max_depth` hops, building a graph
of significant co-enrichment relationships.

## Input format

Whitespace-delimited domain table (e.g. hmmscan `--domtblout`) with at least 19 columns:

- column 1: contig/sequence ID, must match `<name>_ratio_<float>-<index><F|R>`
- column 4: PFAM domain name
- column 5: PFAM accession (e.g. `PF00005.30`)
- columns 18–19: alignment start/end coordinates

## Dependencies

```
pip install scipy statsmodels requests
```

`d3.min.js` is vendored in this repo; if missing, `pfam_coenrichment_network.py` will
download it from cdnjs on first run (requires `requests` and network access).
