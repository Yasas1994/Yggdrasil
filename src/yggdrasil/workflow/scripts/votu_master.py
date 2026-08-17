#!/usr/bin/env python3
"""Assemble the master vOTU table by joining all per-feature outputs.

All joins are on the representative contig id. Optional inputs (taxonomy,
lifestyle, AMG, host) are tolerated as missing -> filled with NA.
"""
import argparse
import csv
import os

import pandas as pd


def read_tsv(path):
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        return pd.read_csv(path, sep="\t")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--votu-clusters", required=True)
    ap.add_argument("--contig2sample", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--quality", default="")
    ap.add_argument("--taxonomy", default="")
    ap.add_argument("--lifestyle", default="")
    ap.add_argument("--amg", default="")
    ap.add_argument("--host", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cl = pd.read_csv(a.votu_clusters, sep="\t")
    reps = cl.loc[cl["is_representative"] == 1, ["votu", "contig"]].rename(
        columns={"contig": "representative"}
    )

    c2s = pd.read_csv(a.contig2sample, sep="\t")
    n_members = cl.groupby("votu").size().rename("n_members")
    n_samples = cl.merge(c2s, on="contig", how="left").groupby("votu")["sample"].nunique().rename("n_samples")

    lengths = pd.read_csv(a.lengths, sep="\t").rename(columns={"contig": "representative", "length": "rep_length"})

    master = reps.merge(n_members, on="votu").merge(n_samples, on="votu").merge(lengths, on="representative", how="left")

    def join(df, path, on, cols):
        sub = read_tsv(path)
        if sub is None or on not in sub.columns:
            for c in cols:
                df[c] = pd.NA
            return df
        keep = [on] + [c for c in cols if c in sub.columns]
        return df.merge(sub[keep].drop_duplicates(on), left_on="representative", right_on=on, how="left").drop(
            columns=[on], errors="ignore"
        )

    master = join(master, a.quality, "contig_id", ["checkv_quality", "completeness", "contamination"])
    master = join(master, a.taxonomy, "id", ["taxonomy"])
    master = join(master, a.lifestyle, "id", ["lifestyle", "lifestyle_score", "bacphlip_lifestyle", "bacphlip_confidence", "bacphlip_virulent", "bacphlip_temperate"])
    master = join(master, a.amg, "id", ["amg_flag", "n_amg_genes", "amg_classes"])
    master = join(master, a.host, "id", ["host", "host_score"])

    master.to_csv(a.out, sep="\t", index=False)
    print(f"[votu_master] rows={len(master)}", file=__import__("sys").stderr)


if __name__ == "__main__":
    main()
