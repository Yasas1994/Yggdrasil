#!/usr/bin/env python3
"""Coerce PhaStyle raw output + BACPHLIP output -> lifestyle.tsv.

PhaStyle writes: sequence_id, fasta_id, predicted_label, score_temperate,
score_virulent. We take fasta_id as the contig id, predicted_label as the call,
and the predicted class's probability as the score. BACPHLIP (id,
bacphlip_virulent, bacphlip_temperate, bacphlip_lifestyle, bacphlip_confidence)
is joined on id when provided. Matching is substring-based so minor header drift
across builds is tolerated.
"""
import argparse
import csv
import os
import sys

import pandas as pd


def find(cols, *keys):
    low = {c: c.lower() for c in cols}
    for k in keys:
        for c in cols:
            if k in low[c]:
                return c
    return None


def norm_label(val):
    v = (val or "").lower()
    if "temp" in v or "lysog" in v:
        return "temperate"
    if "vir" in v or "lytic" in v:
        return "virulent"
    return ""


def parse_phastyle(path):
    rows = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
    out = []
    if not rows:
        return pd.DataFrame(columns=["id", "lifestyle", "lifestyle_score"])
    cols = list(rows[0].keys())
    idc = find(cols, "fasta_id", "contig") or find(cols, "name", "seq", "id")
    pc = find(cols, "label", "prediction", "lifestyle", "class")
    tc = find(cols, "temperate", "lysog")
    vc = find(cols, "virulent", "lytic")
    for row in rows:
        rid = row.get(idc, "") if idc else ""
        life = norm_label(row.get(pc, "")) if pc else ""
        t = row.get(tc, "") if tc else ""
        v = row.get(vc, "") if vc else ""
        score = ""
        if t != "" or v != "":
            try:
                tf, vf = float(t or 0), float(v or 0)
            except ValueError:
                tf = vf = 0.0
            if not life:
                life = "temperate" if tf >= vf else "virulent"
            score = f"{max(tf, vf):.6g}"
        if life in ("temperate", "virulent") and rid:
            out.append((rid, life, score))
    return pd.DataFrame(out, columns=["id", "lifestyle", "lifestyle_score"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--bacphlip", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = parse_phastyle(a.raw)
    if a.bacphlip and os.path.exists(a.bacphlip) and os.path.getsize(a.bacphlip) > 0:
        bp = pd.read_csv(a.bacphlip, sep="\t")
        if "id" in bp.columns and not bp.empty:
            df = df.merge(bp.drop_duplicates("id"), on="id", how="outer")
    for c in ["bacphlip_virulent", "bacphlip_temperate", "bacphlip_lifestyle", "bacphlip_confidence"]:
        if c not in df.columns:
            df[c] = pd.NA
    df = df[["id", "lifestyle", "lifestyle_score", "bacphlip_lifestyle",
             "bacphlip_confidence", "bacphlip_virulent", "bacphlip_temperate"]]
    df.to_csv(a.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
