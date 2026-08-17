#!/usr/bin/env python3
"""Concatenate per-sample phage FASTAs: rename contigs, length-filter, optional dedup.

Runs in the yggdrasil env (needs only biopython). Produces the master contig
FASTA plus the contig->sample and contig->length maps used downstream.
"""
import argparse
import csv
import hashlib
import sys

from Bio import SeqIO


def wrap(seq, w=80):
    return "\n".join(seq[i : i + w] for i in range(0, len(seq), w))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--min-length", type=int, default=5000)
    ap.add_argument("--dedup", action="store_true")
    ap.add_argument("--out-fna", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-lengths", required=True)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.samples, newline=""), delimiter="\t"))
    seen = set()
    n_in = n_out = n_short = n_dup = 0
    with open(a.out_fna, "w") as fa, open(a.out_map, "w") as mp, open(a.out_lengths, "w") as ln:
        mp.write("contig\tsample\n")
        ln.write("contig\tlength\n")
        for r in rows:
            s, g = r["sample"], r.get("genome", "")
            if not g:
                continue
            try:
                recs = list(SeqIO.parse(g, "fasta"))
            except Exception as e:
                print(f"[preprocess] WARN: cannot read {g}: {e}", file=sys.stderr)
                continue
            if not recs:
                print(f"[preprocess] WARN: empty {g}", file=sys.stderr)
                continue
            for rec in recs:
                n_in += 1
                seq = str(rec.seq).upper()
                if len(seq) < a.min_length:
                    n_short += 1
                    continue
                if a.dedup:
                    h = hashlib.md5(seq.encode()).hexdigest()
                    if h in seen:
                        n_dup += 1
                        continue
                    seen.add(h)
                cid = f"{s}__{rec.id}"
                fa.write(f">{cid}\n{wrap(seq)}\n")
                mp.write(f"{cid}\t{s}\n")
                ln.write(f"{cid}\t{len(seq)}\n")
                n_out += 1
    print(f"[preprocess] in={n_in} out={n_out} short={n_short} dup={n_dup}", file=sys.stderr)


if __name__ == "__main__":
    main()
