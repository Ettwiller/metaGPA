#!/usr/bin/env python3
"""
Build a MAFFT alignment and FastTree phylogeny for a PFAM domain, plus an
iTOL-style mapping file coloring leaves by enrichment ratio.

Unlike get_proteins.py, every hit of <pfam_id> in <tab_file> is included
regardless of ratio (duplicates are still removed); <ratio_cutoff> is used
only to color each leaf as "selected" (ratio >= cutoff) or "unselected".
Sequence extraction (--domain_only / --full_length / --include_truncated)
works exactly as in get_proteins.py, whose functions this script reuses.

--hmmer_evalue_max caps hits by the domain's independent E-value (i-Evalue,
column 13 of the tab file): only hits with i-Evalue < this threshold are kept.

By default, MAFFT runs its standard progressive alignment (FFT-NS-2). Pass
--localpair to use L-INS-i instead (more accurate for divergent sequences,
but all-pairs Smith-Waterman makes it much slower on large sequence sets).

Each retained sequence gets a unique name of the form
"<S|unselected>_<n>_<ratio>" (e.g. "S_1_3.7", "unselected_14_0.8"), which is
used as its ID in the FASTA, alignment, tree, and mapping file.

Usage:
    python build_pfam_tree.py --hmmer_output <tab_file> --faa <faa_file> --hmm_id <pfam_id> --ratio <ratio_cutoff> [--domain_only|--full_length] [--include_truncated] [--hmmer_evalue_max VALUE] [--localpair]

Example:
    python build_pfam_tree.py --hmmer_output data.tab --faa proteins.faa --hmm_id PF00709.24 --ratio 10 --domain_only --hmmer_evalue_max 1e-5

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
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--hmmer_output',     required=True, metavar='TAB_FILE',      help='HMMER-format PFAM annotation tab file')
    p.add_argument('--faa',              required=True, metavar='FAA_FILE',      help='Protein FASTA file (frame translations)')
    p.add_argument('--hmm_id',           required=True, metavar='PFAM_ID',       help='PFAM domain ID (e.g. PF00709.24); version suffix ignored', type=lambda x: x.split('.')[0])
    p.add_argument('--ratio',            required=True, metavar='RATIO_CUTOFF',  type=float, help='Enrichment ratio cutoff')
    mode_grp = p.add_mutually_exclusive_group()
    mode_grp.add_argument('--domain_only', action='store_true', help='Extract aligned domain span only (default)')
    mode_grp.add_argument('--full_length', action='store_true', help='Extract full ORF around domain')
    p.add_argument('--include_truncated', action='store_true',  help='Keep sequences without flanking stop codons')
    p.add_argument('--localpair',         action='store_true',  help='Use MAFFT L-INS-i (slower, more accurate)')
    p.add_argument('--hmmer_evalue_max', default=None,  metavar='VALUE',         type=float, help='Max i-Evalue for HMMER hits')
    a = p.parse_args()
    mode = 'full_length' if a.full_length else 'domain_only'
    return a.hmmer_output, a.faa, a.hmm_id, a.ratio, mode, a.include_truncated, a.hmmer_evalue_max, a.localpair


def extract_all_sequences(tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max=None):
    """Return a list of {name, ratio, selected} dicts, one per unique sequence,
    in the order hits were encountered in tab_file. No ratio filtering (all
    hits kept); exact-duplicate sequences are collapsed to their first
    occurrence.
    """
    hits = gp.find_hits(tab_file, pfam_id, None, hmmer_evalue_max)
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
    tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max, localpair = parse_args()

    print(f"PFAM           : {pfam_id}")
    print(f"Mode           : {mode}")
    print(f"Ratio cutoff   : >= {ratio_cutoff} (selected vs. unselected coloring only)")
    if hmmer_evalue_max is not None:
        print(f"E-value cutoff : i-Evalue < {hmmer_evalue_max}")
    print(f"MAFFT mode     : {'L-INS-i (--localpair)' if localpair else 'FFT-NS-2 (default)'}")

    records = extract_all_sequences(tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max)
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

    mafft_cmd = [mafft_bin, '--maxiterate', '2']
    if localpair:
        mafft_cmd.append('--localpair')
    mafft_cmd.append(seq_path)
    run(mafft_cmd, aligned_path, "MAFFT")
    run([fasttree_bin, '-gamma', aligned_path], tree_path, "FastTree")

    with open(mapping_path, 'w') as f:
        f.write("name\tleaf_dot_color\tbranch_color\tbar1_height\tbar1_gradient\n")
        for r in records:
            dot_color   = 'bp_green' if r['selected'] else 'k_grey'
            label_color = 'ptm_rose' if r['selected'] else 'ptm_sand'
            f.write(f"{r['name']}\t{dot_color}\t{label_color}\t{r['ratio']:.1f}\tPurples\n")

    print(f"Alignment      : {aligned_path}")
    print(f"Tree           : {tree_path}")
    print(f"Mapping        : {mapping_path}")


if __name__ == '__main__':
    main()
