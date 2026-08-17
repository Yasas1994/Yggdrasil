#!/usr/bin/env python3
"""Write genome_id, length (Kb) for each record in a FASTA (vContact3 -l input)."""
import argparse
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--kb", action="store_true", help="emit length in kilobases")
    args = ap.parse_args()

    name, n = None, 0
    with open(args.fasta) as fh, open(args.out, "w") as out:
        out.write("genome_id\tlength\n")
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    out.write(f"{name}\t{n / 1000 if args.kb else n:.3f}\n")
                name, n = line[1:].split()[0], 0
            else:
                n += len(line)
        if name is not None:
            out.write(f"{name}\t{n / 1000 if args.kb else n:.3f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
