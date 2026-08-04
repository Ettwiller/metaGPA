#!/usr/bin/env python3
"""
Extract amino-acid sequences for a PFAM domain, ranked by enrichment ratio.

For every hit of <pfam_id> in <tab_file>, looks up the corresponding raw
frame-translation sequence in <faa_file> (matched by the contig ID in
column 1) and extracts either just the aligned domain span or the
full-length ORF around it. Sequences are written out ranked from highest
to lowest contig ratio; if two or more hits yield the exact same sequence,
only the highest-ratio instance is kept.

Since <faa_file> is a whole-contig frame translation (not called ORFs) and
contains '*' stop codons, the ORF around the domain is located as: start =
first 'M' after the nearest upstream stop codon, end = the nearest
downstream stop codon (exclusive). This check applies in both modes: if
the domain's surrounding ORF isn't flanked by a stop codon on one or both
sides, that record is skipped by default (even in --domain_only mode);
pass --include_truncated to keep it (best-effort boundary) with
"(truncated)" noted in its header. --full_length outputs the full ORF span;
--domain_only still outputs just the aligned domain span.

--hmmer_evalue_max caps hits by the domain's independent E-value (i-Evalue,
column 13 of the tab file): only hits with i-Evalue < this threshold are kept.

Usage:
    python get_proteins.py --hmmer_output <tab_file> --faa <faa_file> --hmm_id <pfam_id> [--ratio <ratio_cutoff>] [--domain_only|--full_length] [--include_truncated] [--hmmer_evalue_max VALUE]

Example:
    python get_proteins.py --hmmer_output data.tab --faa proteins.faa --hmm_id PF00709.24 --ratio 10 --domain_only --hmmer_evalue_max 1e-5
"""

import re, sys, os


def parse_args():
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--hmmer_output',     required=True, metavar='TAB_FILE',      help='HMMER-format PFAM annotation tab file')
    p.add_argument('--faa',              required=True, metavar='FAA_FILE',      help='Protein FASTA file (frame translations)')
    p.add_argument('--hmm_id',           required=True, metavar='PFAM_ID',       help='PFAM domain ID (e.g. PF00709.24)')
    p.add_argument('--ratio',            default=None,  metavar='RATIO_CUTOFF',  type=float, help='Enrichment ratio cutoff (optional)')
    mode_grp = p.add_mutually_exclusive_group()
    mode_grp.add_argument('--domain_only', action='store_true', help='Extract aligned domain span only (default)')
    mode_grp.add_argument('--full_length', action='store_true', help='Extract full ORF around domain')
    p.add_argument('--include_truncated', action='store_true',  help='Keep sequences without flanking stop codons')
    p.add_argument('--hmmer_evalue_max', default=None,  metavar='VALUE',         type=float, help='Max i-Evalue for HMMER hits')
    a = p.parse_args()
    mode = 'full_length' if a.full_length else 'domain_only'
    return a.hmmer_output, a.faa, a.hmm_id, a.ratio, mode, a.include_truncated, a.hmmer_evalue_max


def get_ratio(contig_id):
    m = re.search(r'ratio_([\d.]+)', contig_id)
    return float(m.group(1)) if m else 0.0


def find_hits(tab_file, pfam_id, ratio_cutoff, hmmer_evalue_max=None):
    """Return one entry per domain hit for pfam_id (prefix match, version-agnostic)."""
    if not os.path.exists(tab_file):
        print(f"Error: file not found: {tab_file}", file=sys.stderr)
        sys.exit(1)

    prefix = pfam_id.split('.')[0]
    hits = []
    with open(tab_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 19:
                continue
            if not parts[4].startswith(prefix):
                continue
            if hmmer_evalue_max is not None:
                try:
                    i_evalue = float(parts[12])
                except ValueError:
                    continue
                if i_evalue >= hmmer_evalue_max:
                    continue
            contig_id = parts[0]
            ratio = get_ratio(contig_id)
            if ratio_cutoff is not None and ratio < ratio_cutoff:
                continue
            hits.append({
                'contig_id': contig_id,
                'pfam_full': parts[4],
                'start':     int(parts[17]),
                'end':       int(parts[18]),
                'dom_num':   parts[9],
                'dom_of':    parts[10],
                'ratio':     ratio,
            })
    return hits


def extract_sequences(faa_file, wanted_ids):
    """Stream the .faa file once, returning {id: sequence} for wanted_ids.

    Standard FASTA: records start with '>id', sequence may wrap over any
    number of lines. Header boundaries are detected via '>', not by
    membership in wanted_ids, so unwanted records in between are skipped
    cleanly rather than being appended to the previous sequence.
    """
    if not os.path.exists(faa_file):
        print(f"Error: file not found: {faa_file}", file=sys.stderr)
        sys.exit(1)

    seqs = {}
    remaining = set(wanted_ids)
    current_id, current_seq = None, []

    with open(faa_file) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id in remaining:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id in remaining:
            seqs[current_id] = ''.join(current_seq)

    return seqs


def find_orf_bounds(seq, start, end):
    """Locate the ORF around a 1-based inclusive domain span [start, end]
    within a raw frame-translation sequence containing '*' stop codons.

    Returns (orf_start, orf_end, flanked) as a 0-based Python slice range;
    flanked is False if either boundary had to fall back to a non-stop-codon
    edge (sequence end, or a stop codon with no 'M' before the domain).
    """
    domain_start0 = start - 1
    domain_end0   = end  # exclusive

    up_stop = seq.rfind('*', 0, domain_start0)
    if up_stop == -1:
        upstream_ok = False
        orf_start = 0
    else:
        m_idx = seq.find('M', up_stop + 1, domain_start0)
        if m_idx == -1:
            upstream_ok = False
            orf_start = up_stop + 1
        else:
            upstream_ok = True
            orf_start = m_idx

    down_stop = seq.find('*', domain_end0)
    if down_stop == -1:
        downstream_ok = False
        orf_end = len(seq)
    else:
        downstream_ok = True
        orf_end = down_stop

    return orf_start, orf_end, upstream_ok and downstream_ok


def main():
    tab_file, faa_file, pfam_id, ratio_cutoff, mode, include_truncated, hmmer_evalue_max = parse_args()

    print(f"PFAM           : {pfam_id}")
    print(f"Mode           : {mode}")
    if ratio_cutoff is not None:
        print(f"Ratio cutoff   : >= {ratio_cutoff}")
    if hmmer_evalue_max is not None:
        print(f"E-value cutoff : i-Evalue < {hmmer_evalue_max}")

    hits = find_hits(tab_file, pfam_id, ratio_cutoff, hmmer_evalue_max)
    print(f"Domain hits    : {len(hits)}")
    if not hits:
        print("No hits found.", file=sys.stderr)
        sys.exit(1)

    hits.sort(key=lambda h: -h['ratio'])

    seqs = extract_sequences(faa_file, {h['contig_id'] for h in hits})
    missing = [h['contig_id'] for h in hits if h['contig_id'] not in seqs]
    if missing:
        print(f"Warning: {len(missing)} contig(s) not found in {faa_file}, skipped.",
              file=sys.stderr)

    base_name = os.path.splitext(os.path.basename(tab_file))[0]
    cutoff_tag = f"_r{ratio_cutoff}" if ratio_cutoff is not None else ""
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(tab_file)),
        f"{base_name}_{pfam_id.replace('.', '_')}{cutoff_tag}_{mode}.fasta"
    )

    n_written = 0
    n_truncated = 0
    n_skipped_truncated = 0
    n_duplicates = 0
    seen_seqs = set()
    with open(out_path, 'w') as out:
        for h in hits:
            protein = seqs.get(h['contig_id'])
            if protein is None:
                continue

            orf_start, orf_end, flanked = find_orf_bounds(protein, h['start'], h['end'])
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

            header = (f">{h['contig_id']} ratio={h['ratio']} pfam={h['pfam_full']} "
                      f"domain={h['start']}-{h['end']} dom={h['dom_num']}/{h['dom_of']}")
            if truncated:
                header += " (truncated)"
                n_truncated += 1
            out.write(header + "\n")
            out.write(seq + "\n")
            n_written += 1

    print(f"Sequences      : {n_written} written"
          + (f" ({n_truncated} truncated)" if n_truncated else ""))
    if n_duplicates:
        print(f"Deduplicated   : {n_duplicates} identical sequence(s) dropped "
              f"(kept highest-ratio instance)")
    if n_skipped_truncated:
        print(f"Skipped        : {n_skipped_truncated} not flanked by stop codons "
              f"(use --include_truncated to keep)")
    print(f"Output         : {out_path}")


if __name__ == '__main__':
    main()
