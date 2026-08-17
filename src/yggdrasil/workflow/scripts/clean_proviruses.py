#!/usr/bin/env python3
"""Remove host regions from CheckV proviruses.

CheckV `viruses.fna` holds non-provirus viral contigs; `proviruses.fna` holds
provirus sequences with host flanks already excised (headers carry the original
coordinates, e.g. `>contig_1 1-37388/47765`). This script concatenates both into
a cleaned contig set and normalizes provirus IDs back to the original contig IDs
so downstream sample mapping, clustering, and annotation stay consistent.
"""

import argparse
import re
import sys


def parse_fasta(path):
    """Yield (header, sequence) tuples."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)


def norm_provirus_id(header: str) -> str:
    """Strip CheckV provirus suffix and coordinate annotation.

    `Kenneth004__Kenneth004_37_1 1-37388/47765` -> `Kenneth004__Kenneth004_37`
    `Kenneth004__Kenneth004_42`               -> `Kenneth004__Kenneth004_42`
    """
    # drop the coordinate comment (everything after the first whitespace)
    tok = header.split()[0]
    # CheckV appends _1, _2, ... for multiple proviruses; strip the trailing _N
    return re.sub(r"_\d+$", "", tok)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--viruses", required=True, help="CheckV viruses.fna")
    ap.add_argument("--proviruses", required=True, help="CheckV proviruses.fna")
    ap.add_argument("--c2s", required=True, help="contig2sample.tsv from preprocess")
    ap.add_argument("--out-fna", required=True)
    ap.add_argument("--out-c2s", required=True)
    ap.add_argument("--out-lengths", required=True)
    args = ap.parse_args()

    records = []
    n_virus = n_provirus = 0

    for hdr, seq in parse_fasta(args.viruses):
        records.append((hdr, seq))
        n_virus += 1

    for hdr, seq in parse_fasta(args.proviruses):
        new_id = norm_provirus_id(hdr)
        records.append((new_id, seq))
        n_provirus += 1

    # write cleaned FASTA
    with open(args.out_fna, "w") as out:
        for hdr, seq in records:
            out.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")

    # lengths
    with open(args.out_lengths, "w") as out:
        out.write("contig\tlength\n")
        for hdr, seq in records:
            out.write(f"{hdr}\t{len(seq)}\n")

    # contig2sample: keep only contigs present in the cleaned set, preserve order
    cleaned_ids = {hdr for hdr, _ in records}
    with open(args.c2s) as inp, open(args.out_c2s, "w") as out:
        header = inp.readline()
        out.write(header)
        for line in inp:
            cid = line.split("\t", 1)[0]
            if cid in cleaned_ids:
                out.write(line)

    print(
        f"[clean] viruses={n_virus} proviruses={n_provirus} total={len(records)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
