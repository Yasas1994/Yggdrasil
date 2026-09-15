#!/usr/bin/env python3
"""Split a multi-FASTA into chunks of at most K sequences.

Stdlib-only (runs in the yggdrasil env). Used by checkpoint rules to fan out
monolithic steps (Pharokka, PhaStyle, BACPHLIP, iPHoP) into per-chunk jobs.
The number of chunks follows from the input size: ceil(n_seqs / seqs_per_chunk).

Writes <outdir>/<prefix>_chunkNNN.fna (NNN zero-padded) plus <outdir>/chunks.txt
listing the chunk ids (chunkNNN), one per line. An empty input yields zero
chunks and an empty chunks.txt; downstream merge rules must tolerate that.
"""
import argparse
import os
import sys


def iter_fasta(path):
    """Yield (header_line, seq_lines) records; header_line includes '>'."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header, seq = line, []
            elif header is not None:
                seq.append(line)
    if header is not None:
        yield header, seq


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--seqs-per-chunk", type=int, required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    a = ap.parse_args()

    if a.seqs_per_chunk < 1:
        print("[split_fasta] --seqs-per-chunk must be >= 1", file=sys.stderr)
        return 1
    os.makedirs(a.outdir, exist_ok=True)

    ids = []
    fh = None
    n = total = 0
    for header, seq in iter_fasta(a.fasta):
        if n % a.seqs_per_chunk == 0:
            if fh is not None:
                fh.close()
            cid = f"chunk{len(ids) + 1:03d}"
            ids.append(cid)
            fh = open(os.path.join(a.outdir, f"{a.prefix}_{cid}.fna"), "w")
        fh.write(header)
        fh.writelines(seq)
        n += 1
        total += 1
    if fh is not None:
        fh.close()
    with open(os.path.join(a.outdir, "chunks.txt"), "w") as fh:
        fh.write("".join(f"{c}\n" for c in ids))
    print(f"[split_fasta] seqs={total} chunks={len(ids)} outdir={a.outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
