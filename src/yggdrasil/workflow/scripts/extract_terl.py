#!/usr/bin/env python3
"""Extract terminase large-subunit (TerL) proteins from an annotate run.

Inputs (produced by rules/annotate.smk):
  gene_families.tsv  columns: contig, gene_id, family
                     `family` is an ID only (e.g. phrog_742), not free text.
  proteins.faa       header's first whitespace token == gene_id; the pharokka path
                     may also carry a free-text description after the id.

TerL detection ORs two signals (missing either one just yields fewer hits):
  1. regex match of --pattern in the protein FASTA defline (fires on pharokka output);
  2. family-id -> description match via optional --family-map (a PHROG/KO annotation
     TSV: family_id<TAB>description). Without a map, family IDs carry no text, so only
     signal (1) fires. If too few TerL proteins are found (<4) the phylo rule writes a
     placeholder tree instead of failing.
"""
import argparse
import re
import sys

import pandas as pd
from Bio import SeqIO


def load_family_map(path):
    m = {}
    if not path:
        return m
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pattern", default=r"(?i)terminase|terl|large terminase|terminase large")
    ap.add_argument("--family-map", default="")
    a = ap.parse_args()

    rx = re.compile(a.pattern)
    fam_desc = load_family_map(a.family_map)

    terl_ids = set()
    try:
        fam = pd.read_csv(a.families, sep="\t", dtype=str).fillna("")
    except Exception:
        fam = pd.DataFrame(columns=["gene_id", "family"])
    for _, row in fam.iterrows():
        fid = str(row.get("family", ""))
        desc = fam_desc.get(fid, "")  # bare IDs (phrog_742) match only with a map
        if rx.search(desc):
            terl_ids.add(str(row.get("gene_id", "")))

    n = 0
    with open(a.out, "w") as out:
        try:
            recs = list(SeqIO.parse(a.proteins, "fasta"))
        except Exception:
            recs = []
        for rec in recs:
            if rec.id in terl_ids or rx.search(rec.description or ""):
                out.write(f">{rec.id}\n{str(rec.seq)}\n")
                n += 1
    print(f"[extract_terl] TerL proteins={n}", file=sys.stderr)


if __name__ == "__main__":
    main()
