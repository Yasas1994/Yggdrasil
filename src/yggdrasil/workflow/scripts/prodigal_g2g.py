#!/usr/bin/env python3
"""gene2genome.tsv (protein_id, genome_id, keywords) for the prodigal branch.

prodigal keeps the original contig id and appends _<idx>; genome = contig with
the trailing _<idx> removed. keywords = best PHROG target from phrog.m8 (or None).
"""
import argparse
import csv
import sys


def best_hits(m8):
    hit = {}
    for r in csv.reader(open(m8), delimiter="\t"):
        if len(r) < 5:
            continue
        try:
            b = float(r[4])
        except ValueError:
            continue
        if r[0] not in hit or b > hit[r[0]][0]:
            hit[r[0]] = (b, r[1])
    return {k: v[1] for k, v in hit.items()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--m8", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    hits = best_hits(a.m8) if a.m8 else {}
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "genome_id", "keywords"])
        for line in open(a.faa):
            if not line.startswith(">"):
                continue
            gid = line[1:].split()[0]
            ctg = gid.rsplit("_", 1)[0]
            w.writerow([gid, ctg, hits.get(gid, "None") or "None"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
