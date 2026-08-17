#!/usr/bin/env python3
"""Build annotation-derived tables from Pharokka output.

Pharokka renames genes to <HASH>_CDS_NNNN (one hash per *run*, shared across
all contigs); the authoritative gene -> contig and gene -> PHROG mapping lives
in pharokka_cds_final_merged_output.tsv.

  --out           gene_families.tsv  (contig, gene_id, family, product, function) PHROG hits only
  --out-g2g       gene2genome.tsv    (protein_id, genome_id) for vContact3
"""
import argparse
import csv
import os
import sys

NO_HIT = {"no_phrog", "no hit", "nan", "none", "-", ""}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="pharokka output directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-g2g", default=None)
    args = ap.parse_args()

    merged = os.path.join(args.indir, "pharokka_cds_final_merged_output.tsv")
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig", "gene_id", "family", "product", "function"])
        if not os.path.exists(merged):
            if args.out_g2g:
                with open(args.out_g2g, "w", newline="") as g2:
                    csv.writer(g2, delimiter="\t").writerow(["protein_id", "genome_id"])
            return 0
        with open(merged, newline="") as inf:
            r = csv.DictReader(inf, delimiter="\t")
            cols = list(r.fieldnames or [])
            gc = "gene" if "gene" in cols else cols[0]
            cc = "contig" if "contig" in cols else None
            pcols = [c for c in cols if "phrog" in c.lower()]
            pc = "phrog" if "phrog" in cols else (pcols[0] if pcols else None)
            ac = "annot" if "annot" in cols else None       # product (e.g. "terminase large subunit")
            fc = "category" if "category" in cols else None  # function category (e.g. "head and packaging")
            g2 = csv.writer(open(args.out_g2g, "w", newline=""), delimiter="\t") if args.out_g2g else None
            if g2:
                g2.writerow(["protein_id", "genome_id", "keywords"])
            for row in r:
                g = row.get(gc, "")
                c = row.get(cc, "") if cc else g.rsplit("_", 1)[0]
                v = (row.get(pc, "") or "").strip() if pc else ""
                hit = v.lower() not in NO_HIT
                fam = ("PHROG_" + v) if hit else "None"
                if g2 and g:
                    g2.writerow([g, c, fam])
                if g and hit:
                    prod = row.get(ac, "") if ac else ""
                    func = row.get(fc, "") if fc else ""
                    w.writerow([c, g, fam, prod, func])
    return 0


if __name__ == "__main__":
    sys.exit(main())
