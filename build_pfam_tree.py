#!/usr/bin/env python3
"""
Build a MAFFT alignment and FastTree phylogeny for a PFAM domain, plus an
iTOL-style mapping file coloring leaves by enrichment ratio.

Unlike get_proteins.py, every hit of <pfam_id> in <tab_file> is included
regardless of ratio (duplicates are still removed); <ratio_cutoff> is used
only to color each leaf as "selected" (ratio >= cutoff) or "unselected".
Sequence extraction (--domain_only / --full_length / --include_truncated)
works exactly as in get_proteins.py, whose functions this script reuses.

Each retained sequence gets a unique name of the form
"<S|unselected>_<n>_<ratio>" (e.g. "S_1_3.7", "unselected_14_0.8"), which is
used as its ID in the FASTA, alignment, tree, and mapping file.

Usage:
    python build_pfam_tree.py <tab_file> <faa_file> <pfam_id> <ratio_cutoff> [--domain_only|--full_length] [--include_truncated]

Example:
    python build_pfam_tree.py data.tab proteins.faa PF00709.24 10 --domain_only

Requires `mafft` and `FastTree` (or `fasttree`) on PATH.
"""

import sys, os, shutil, subprocess, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, f'{name}.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gp = _load('get_proteins')  # find_hits, extract_sequences, find_orf_bounds


def parse_args():
    args = sys.argv[1:]

    mode = 'domain_only'
    if '--full_length' in args:
        mode = 'full_length'
        args.remove('--full_length')
    if '--domain_only' in args:
        mode = 'domain_only'
        args.remove('--domain_only')

    include_truncated = '--include_truncated' in args
    if include_truncated:
        args.remove('--include_truncated')

    if len(args) != 4:
        print(__doc__)
        sys.exit(1)

    tab_file, faa_file, pfam_id, ratio_cutoff = args[0], args[1], args[2], float(args[3])
    return tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated


def extract_all_sequences(tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated):
    """Return a list of {name, ratio, selected} dicts, one per unique sequence,
    in the order hits were encountered in tab_file. No ratio filtering (all
    hits kept); exact-duplicate sequences are collapsed to their first
    occurrence.
    """
    hits = gp.find_hits(tab_file, pfam_id, None)
    if not hits:
        print("No hits found.", file=sys.stderr)
        sys.exit(1)

    seqs = gp.extract_sequences(faa_file, {h['contig_id'] for h in hits})
    missing = [h['contig_id'] for h in hits if h['contig_id'] not in seqs]
    if missing:
        print(f"Warning: {len(missing)} contig(s) not found in {faa_file}, skipped.",
              file=sys.stderr)

    counters = {'S': 0, 'unselected': 0}
    seen_seqs = set()
    records = []
    n_truncated = 0
    n_skipped_truncated = 0
    n_duplicates = 0

    for h in hits:
        protein = seqs.get(h['contig_id'])
        if protein is None:
            continue

        orf_start, orf_end, flanked = gp.find_orf_bounds(protein, h['start'], h['end'])
        truncated = not flanked
        if truncated and not include_truncated:
            n_skipped_truncated += 1
            continue

        if mode == 'domain_only':
            start = max(h['start'], 1)
            end   = min(h['end'], len(protein))
            seq = protein[start - 1:end]
        else:
            seq = protein[orf_start:orf_end]

        if seq in seen_seqs:
            n_duplicates += 1
            continue
        seen_seqs.add(seq)

        if truncated:
            n_truncated += 1

        selected = h['ratio'] >= ratio_cutoff
        group = 'S' if selected else 'unselected'
        counters[group] += 1
        name = f"{group}_{counters[group]}_{h['ratio']:.1f}"

        records.append({
            'name': name,
            'ratio': h['ratio'],
            'selected': selected,
            'seq': seq,
            'truncated': truncated,
        })

    print(f"Sequences      : {len(records)} unique"
          + (f" ({n_truncated} truncated)" if n_truncated else ""))
    if n_duplicates:
        print(f"Deduplicated   : {n_duplicates} identical sequence(s) dropped")
    if n_skipped_truncated:
        print(f"Skipped        : {n_skipped_truncated} not flanked by stop codons "
              f"(use --include_truncated to keep)")

    return records


def find_tool(*names):
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def run(cmd, stdout_path, step_name):
    print(f"Running {step_name}: {' '.join(cmd)}")
    with open(stdout_path, 'w') as out:
        # stderr is left to inherit the parent's, so progress output
        # (e.g. MAFFT's) streams to the terminal live instead of being
        # captured silently.
        result = subprocess.run(cmd, stdout=out)
    if result.returncode != 0:
        print(f"Error: {step_name} failed (exit {result.returncode}).", file=sys.stderr)
        sys.exit(1)


def main():
    tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated = parse_args()

    print(f"PFAM           : {pfam_id}")
    print(f"Mode           : {mode}")
    print(f"Ratio cutoff   : >= {ratio_cutoff} (selected vs. unselected coloring only)")

    records = extract_all_sequences(tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated)
    if len(records) < 2:
        print("Error: need at least 2 unique sequences to align/build a tree.", file=sys.stderr)
        sys.exit(1)

    mafft_bin = find_tool('mafft')
    if mafft_bin is None:
        print("Error: 'mafft' not found on PATH. Install it (e.g. `conda install -c bioconda mafft`).",
              file=sys.stderr)
        sys.exit(1)

    fasttree_bin = find_tool('FastTree', 'fasttree', 'FastTreeMP')
    if fasttree_bin is None:
        print("Error: FastTree not found on PATH. Install it "
              "(e.g. `conda install -c bioconda fasttree`).", file=sys.stderr)
        sys.exit(1)

    base_name = os.path.splitext(os.path.basename(tab_file))[0]
    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(tab_file)),
        f"{base_name}_{pfam_id.replace('.', '_')}_r{ratio_cutoff}_{mode}_tree"
    )
    os.makedirs(out_dir, exist_ok=True)
    print(f"Output dir     : {out_dir}")

    seq_path     = os.path.join(out_dir, "sequences.fasta")
    aligned_path = os.path.join(out_dir, "aligned.fasta")
    tree_path    = os.path.join(out_dir, "tree.nwk")
    mapping_path = os.path.join(out_dir, "mapping.txt")

    with open(seq_path, 'w') as f:
        for r in records:
            f.write(f">{r['name']}\n{r['seq']}\n")

    run([mafft_bin, '--maxiterate', '2', '--localpair', seq_path], aligned_path, "MAFFT")
    run([fasttree_bin, '-gamma', aligned_path], tree_path, "FastTree")

    with open(mapping_path, 'w') as f:
        f.write("name\tleaf_dot_color\tleaf_label_color\tbar1_height\tbar1_gradient\n")
        for r in records:
            dot_color   = 'bp_green' if r['selected'] else 'k_grey'
            label_color = 'ptm_rose' if r['selected'] else 'ptm_sand'
            f.write(f"{r['name']}\t{dot_color}\t{label_color}\t{r['ratio']:.1f}\tPurples\n")

    print(f"Alignment      : {aligned_path}")
    print(f"Tree           : {tree_path}")
    print(f"Mapping        : {mapping_path}")


if __name__ == '__main__':
    main()
