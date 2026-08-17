#!/usr/bin/env python3
"""Build vOTU x sample (presence/absence) and gene-family x sample (count) matrices (wide)."""
import argparse
import csv
import os
import sys

import pandas as pd


def read_map(path, key, val):
    d = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            d[row[key]] = row[val]
    return d


def sample_columns(samples_path):
    return [r["sample"] for r in csv.DictReader(open(samples_path, newline=""), delimiter="\t")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contig2sample", required=True)
    ap.add_argument("--votu-clusters", required=True)
    ap.add_argument("--gene-families", required=True)   # contig, gene_id, family (may be empty)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-votu", required=True)
    ap.add_argument("--out-gene", required=True)
    a = ap.parse_args()

    samples = sample_columns(a.samples)
    c2s = read_map(a.contig2sample, "contig", "sample")

    # ---- vOTU x sample ----
    cl = pd.read_csv(a.votu_clusters, sep="\t")
    cl["sample"] = cl["contig"].map(c2s)
    vm = (
        cl.dropna(subset=["sample"])
        .groupby(["votu", "sample"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=samples, fill_value=0)
    )
    # Counts of member contigs are not biologically meaningful here; collapse
    # to presence/absence (1 = vOTU detected in the sample).
    vm = (vm > 0).astype(int)
    vm.index.name = "votu"
    vm.reset_index().to_csv(a.out_votu, sep="\t", index=False)

    # ---- gene family x sample ----
    if os.path.exists(a.gene_families) and os.path.getsize(a.gene_families) > 0:
        gf = pd.read_csv(a.gene_families, sep="\t")
        if {"contig", "family"}.issubset(gf.columns):
            gf["sample"] = gf["contig"].map(c2s)
            gm = (
                gf.dropna(subset=["sample", "family"])
                .groupby(["family", "sample"])
                .size()
                .unstack(fill_value=0)
                .reindex(columns=samples, fill_value=0)
            )
            gm.index.name = "family"
            gm.reset_index().to_csv(a.out_gene, sep="\t", index=False)
        else:
            print("[build_count_tables] gene_families missing contig/family columns", file=sys.stderr)
            pd.DataFrame({"family": []}).to_csv(a.out_gene, sep="\t", index=False)
    else:
        pd.DataFrame({"family": []}).to_csv(a.out_gene, sep="\t", index=False)
    print(f"[build_count_tables] votu_matrix={vm.shape}", file=sys.stderr)


if __name__ == "__main__":
    main()
