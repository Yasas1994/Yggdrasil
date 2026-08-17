#!/usr/bin/env python3
"""Pick one representative genome per vclust cluster.

Representative = best CheckV quality tier, tie-broken by longest contig.
Tolerant of missing/empty CheckV output (falls back to length-only).

Outputs:
  votu_clusters.tsv       contig, votu, is_representative(0/1)
  votu_representatives.fna
"""
import argparse
import csv
import sys

from Bio import SeqIO

QUALITY_RANK = {
    "complete": 0,
    "high-quality": 1,
    "medium-quality": 2,
    "low-quality": 3,
    "not-determined": 4,
}


def read_clusters(path):
    """Return list of (member, votu). Accepts 2-col TSV (member,cluster) with
    optional header, auto-detected.
    """
    pairs = []
    with open(path, newline="") as fh:
        rows = [r for r in (line.rstrip("\n").split("\t") for line in fh) if r and r[0]]
    if not rows:
        return pairs
    start = 0
    header = [c.lower() for c in rows[0]]
    if any(h in ("cluster", "votu", "member", "sequence", "id", "contig") for h in header):
        start = 1
        # try to locate columns
        try:
            ci = next(i for i, h in enumerate(header) if h in ("cluster", "votu", "cluster_id"))
            mi = next(i for i, h in enumerate(header) if h in ("member", "sequence", "id", "contig"))
        except StopIteration:
            ci, mi = (1, 0) if len(header) >= 2 else (0, 0)
    else:
        mi, ci = 0, 1
    for r in rows[start:]:
        if len(r) <= max(mi, ci):
            continue
        pairs.append((r[mi], r[ci]))
    return pairs


def read_quality(path):
    rank = {}
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                cid = row.get("contig_id") or row.get("contig") or row.get("id")
                q = (row.get("checkv_quality") or "").strip().lower()
                if cid:
                    rank[cid] = QUALITY_RANK.get(q, 4)
    except FileNotFoundError:
        pass
    return rank


def read_lengths(path):
    d = {}
    with open(path, newline="") as fh:
        next(fh, None)
        for line in fh:
            c, _, l = line.rstrip("\n").partition("\t")
            if c:
                d[c] = int(l or 0)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--quality", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-fasta", required=True)
    a = ap.parse_args()

    pairs = read_clusters(a.clusters)
    if not pairs:
        print("[select_representatives] no clusters parsed", file=sys.stderr)
        open(a.out_clusters, "w").write("contig\tvotu\tis_representative\n")
        open(a.out_fasta, "w").close()
        return

    quality = read_quality(a.quality)
    lengths = read_lengths(a.lengths)

    by_votu: dict[str, list[str]] = {}
    for member, votu in pairs:
        by_votu.setdefault(votu, []).append(member)

    reps = set()
    with open(a.out_clusters, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig", "votu", "is_representative"])
        for votu, members in sorted(by_votu.items()):
            rep = min(members, key=lambda c: (quality.get(c, 4), -lengths.get(c, 0)))
            reps.add(rep)
            for c in members:
                w.writerow([c, votu, int(c == rep)])

    with open(a.out_fasta, "w") as out:
        for rec in SeqIO.parse(a.fasta, "fasta"):
            if rec.id in reps:
                out.write(f">{rec.id}\n{str(rec.seq)}\n")
    print(f"[select_representatives] votus={len(by_votu)} representatives={len(reps)}", file=sys.stderr)


if __name__ == "__main__":
    main()
