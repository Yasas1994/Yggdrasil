#!/usr/bin/env python3
"""Map vContact3 output onto representative contig ids -> taxonomy.tsv (id, taxonomy).

vContact3 3.x writes exports/final_assignments.csv with per-rank *_prediction
columns (Reference == False rows are the user genomes). We assemble a lineage
string from those prediction columns in rank order. A tolerant fallback handles
older/other layouts (a single taxonomy column, or genus/family/... columns).
One row is emitted per representative; anything unplaced -> 'Unclassified'.
"""
import argparse
import glob
import os
import sys

import pandas as pd

ID_COLS = ("genome", "genome_id", "sequence", "seq_id", "id", "name", "contig")
TAX_COLS = ("taxonomy", "taxon", "genus", "family", "order", "vc subcluster", "vc",
            "subcluster", "cluster", "prediction", "assignment")
RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "subfamily",
         "genus", "species")
BAD = ("", "nan", "none", "unclassified", "unknown", "-")


def rep_ids(fasta):
    return [line[1:].split()[0] for line in open(fasta) if line.startswith(">")]


def candidates(vc_dir):
    preferred = [
        os.path.join(vc_dir, "exports", "final_assignments.csv"),
        os.path.join(vc_dir, "genome_by_genome_profile.csv"),
    ]
    rest = [p for p in glob.glob(os.path.join(vc_dir, "**", "*.csv"), recursive=True)
            if os.path.isfile(p) and p not in preferred]
    return [p for p in preferred if os.path.isfile(p)] + rest


def lower_map(cols):
    return {c.lower().strip(): c for c in cols}


def prediction_cols(df, lm):
    """vContact3 *_prediction columns, ordered by taxonomic rank."""
    by_rank = {r: lm.get(r + "_prediction") for r in RANKS}
    ordered = [by_rank[r] for r in RANKS if by_rank.get(r)]
    if ordered:
        return ordered
    # any other *_prediction column, stable order
    return [lm[k] for k in lm if k.endswith("_prediction")]


def extract(path):
    df = pd.read_csv(path)
    if df.empty:
        return None
    lm = lower_map(df.columns)
    idc = next((lm[k] for k in ID_COLS if k in lm), None) or df.columns[0]

    # user genomes only when a Reference flag exists (drops the DB reference rows)
    refc = lm.get("reference")
    if refc is not None:
        mask = df[refc].astype(str).str.strip().str.lower().isin(("false", "0", "no"))
        if mask.any():
            df = df[mask]

    tcs = prediction_cols(df, lm)
    if not tcs:
        if "taxonomy" in lm and lm["taxonomy"] != idc:
            tcs = [lm["taxonomy"]]
        else:
            tcs = [lm[k] for k in TAX_COLS if k in lm and lm[k] != idc]
    if not tcs:
        return None

    out = {}
    for _, row in df.iterrows():
        gid = str(row[idc])
        if not gid or gid.lower() == "nan":
            continue
        parts = [str(row[c]).strip() for c in tcs
                 if str(row[c]).strip().lower() not in BAD]
        out[gid] = ";".join(parts) if parts else "Unclassified"
    return out or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcontact3-dir", required=True)
    ap.add_argument("--representatives", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    reps = rep_ids(a.representatives)
    assignments = {}
    for p in candidates(a.vcontact3_dir):
        try:
            got = extract(p)
        except Exception as e:  # noqa: BLE001 - keep scanning other files
            print(f"[normalize_taxonomy] skip {p}: {e}", file=sys.stderr)
            continue
        if got:
            assignments = got
            print(f"[normalize_taxonomy] using {p} ({len(got)} genomes)", file=sys.stderr)
            break

    with open(a.out, "w", newline="") as fh:
        fh.write("id\ttaxonomy\n")
        for rid in reps:
            fh.write(f"{rid}\t{assignments.get(rid, 'Unclassified')}\n")
    print(f"[normalize_taxonomy] rows={len(reps)} placed={len(set(reps) & set(assignments))}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
