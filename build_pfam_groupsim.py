#!/usr/bin/env python3
"""
Build a MAFFT alignment for a PFAM domain and run GroupSim to find
specificity-determining positions between two ratio-defined groups.

Like build_pfam_tree.py, every hit of <pfam_id> in <tab_file> is included
regardless of ratio (duplicates are still removed); <ratio_cutoff> splits
sequences into two groups for GroupSim: group1 = ratio < cutoff, group2 =
ratio >= cutoff. Sequence extraction (--domain_only / --full_length /
--include_truncated) works exactly as in get_proteins.py, whose functions
this script reuses.

--hmmer_evalue_max caps hits by the domain's independent E-value (i-Evalue,
column 13 of the tab file): only hits with i-Evalue < this threshold are kept.

Each retained sequence gets a unique name of the form
"<group1|group2>_<n>_<ratio>" (e.g. "group2_1_17.7", "group1_9_4.3"), which
is used as its ID in the FASTA, alignment, and GroupSim's manual group file.

Usage:
    python build_pfam_groupsim.py <tab_file> <faa_file> <pfam_id> <ratio_cutoff> [--domain_only|--full_length] [--include_truncated] [--hmmer_evalue_max VALUE]

Example:
    python build_pfam_groupsim.py data.tab proteins.faa PF00709.24 10 --domain_only --hmmer_evalue_max 1e-5

Requires `mafft` on PATH, and GroupSim-py3 (https://github.com/jacgonisa/groupsim-py3)
available locally — point GROUPSIM_SCRIPT at its src/groupsim.py, or have a
`groupsim.py` on PATH. GroupSim-py3 itself requires biopython, numpy, pandas,
scipy, matplotlib, and seaborn in the Python environment running this script.
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

    hmmer_evalue_max = None
    if '--hmmer_evalue_max' in args:
        idx = args.index('--hmmer_evalue_max')
        hmmer_evalue_max = float(args[idx + 1])
        del args[idx:idx + 2]

    if len(args) != 4:
        print(__doc__)
        sys.exit(1)

    tab_file, faa_file, pfam_id, ratio_cutoff = args[0], args[1], args[2], float(args[3])
    return tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max


def extract_all_sequences(tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max=None):
    """Return a list of {name, ratio, group, seq} dicts, one per unique
    sequence, in the order hits were encountered in tab_file. No ratio
    filtering (all hits kept); exact-duplicate sequences are collapsed to
    their first occurrence.
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

    counters = {'group1': 0, 'group2': 0}
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

        group = 'group2' if h['ratio'] >= ratio_cutoff else 'group1'
        counters[group] += 1
        name = f"{group}_{counters[group]}_{h['ratio']:.1f}"

        records.append({'name': name, 'ratio': h['ratio'], 'group': group, 'seq': seq})

    print(f"Sequences      : {len(records)} unique"
          + (f" ({n_truncated} truncated)" if n_truncated else ""))
    if n_duplicates:
        print(f"Deduplicated   : {n_duplicates} identical sequence(s) dropped")
    if n_skipped_truncated:
        print(f"Skipped        : {n_skipped_truncated} not flanked by stop codons "
              f"(use --include_truncated to keep)")
    print(f"Group1 (< {ratio_cutoff}) : {counters['group1']} sequences")
    print(f"Group2 (>= {ratio_cutoff}): {counters['group2']} sequences")

    return records


VENDOR_SCRIPT = os.path.join(HERE, 'vendor', 'groupsim-py3', 'src', 'groupsim.py')
VENDOR_PYTHON = os.path.join(HERE, 'vendor', 'groupsim-env', 'bin', 'python3')


def find_groupsim_script():
    env_path = os.environ.get('GROUPSIM_SCRIPT')
    if env_path and os.path.exists(env_path):
        return env_path
    if os.path.exists(VENDOR_SCRIPT):
        return VENDOR_SCRIPT
    return shutil.which('groupsim.py')


def find_groupsim_python():
    env_path = os.environ.get('GROUPSIM_PYTHON')
    if env_path and os.path.exists(env_path):
        return env_path
    if os.path.exists(VENDOR_PYTHON):
        return VENDOR_PYTHON
    return sys.executable


def run(cmd, step_name, stdout_path=None):
    print(f"Running {step_name}: {' '.join(cmd)}")
    out = open(stdout_path, 'w') if stdout_path else None
    try:
        # stderr is left to inherit the parent's, so tool progress streams
        # to the terminal live instead of being captured silently.
        result = subprocess.run(cmd, stdout=out)
    finally:
        if out:
            out.close()
    if result.returncode != 0:
        print(f"Error: {step_name} failed (exit {result.returncode}).", file=sys.stderr)
        sys.exit(1)


def main():
    tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max = parse_args()

    print(f"PFAM           : {pfam_id}")
    print(f"Mode           : {mode}")
    print(f"Ratio cutoff   : {ratio_cutoff} (group1 < cutoff, group2 >= cutoff)")
    if hmmer_evalue_max is not None:
        print(f"E-value cutoff : i-Evalue < {hmmer_evalue_max}")

    records = extract_all_sequences(tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max)

    n_group1 = sum(1 for r in records if r['group'] == 'group1')
    n_group2 = sum(1 for r in records if r['group'] == 'group2')
    if n_group1 == 0 or n_group2 == 0:
        print("Error: GroupSim needs at least one sequence in each group "
              "(ratio < cutoff and ratio >= cutoff).", file=sys.stderr)
        sys.exit(1)
    if len(records) < 4:
        print("Error: need at least 4 unique sequences (2+ per group) to run GroupSim.",
              file=sys.stderr)
        sys.exit(1)

    mafft_bin = shutil.which('mafft')
    if mafft_bin is None:
        print("Error: 'mafft' not found on PATH. Install it (e.g. `conda install -c bioconda mafft`).",
              file=sys.stderr)
        sys.exit(1)

    groupsim_script = find_groupsim_script()
    if groupsim_script is None:
        print("Error: GroupSim-py3 script not found. Clone "
              "https://github.com/jacgonisa/groupsim-py3 into vendor/groupsim-py3 "
              "next to this script (auto-detected), or put its src/groupsim.py on "
              "PATH (as `groupsim.py`), or set the GROUPSIM_SCRIPT environment "
              "variable to its full path. It also requires: pip install biopython "
              "numpy pandas scipy matplotlib seaborn — see README.md for the "
              "recommended vendor/ + venv setup.", file=sys.stderr)
        sys.exit(1)

    groupsim_python = find_groupsim_python()

    base_name = os.path.splitext(os.path.basename(tab_file))[0]
    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(tab_file)),
        f"{base_name}_{pfam_id.replace('.', '_')}_r{ratio_cutoff}_{mode}_groupsim"
    )
    os.makedirs(out_dir, exist_ok=True)
    print(f"Output dir     : {out_dir}")

    seq_path     = os.path.join(out_dir, "sequences.fasta")
    aligned_path = os.path.join(out_dir, "aligned.fasta")
    groups_path  = os.path.join(out_dir, "groups.txt")
    out_prefix   = os.path.join(out_dir, "groupsim")

    with open(seq_path, 'w') as f:
        for r in records:
            f.write(f">{r['name']}\n{r['seq']}\n")

    run([mafft_bin, '--maxiterate', '2', '--localpair', seq_path], "MAFFT", stdout_path=aligned_path)

    with open(groups_path, 'w') as f:
        for group in ('group1', 'group2'):
            names = [r['name'] for r in records if r['group'] == group]
            f.write(f"{group}: {', '.join(names)}\n")

    run([groupsim_python, groupsim_script, '-k', groups_path, '-o', out_prefix, aligned_path],
        "GroupSim")

    print(f"Alignment      : {aligned_path}")
    print(f"Groups         : {groups_path}")
    print(f"Scores         : {out_prefix}.txt")
    print(f"Plot           : {out_prefix}_manhattan_plot.png")


if __name__ == '__main__':
    main()
