#!/usr/bin/env python3
"""Normalize an iPHoP host-prediction CSV -> iphop.tsv (id, host, host_score).

iPHoP 1.3 writes Host_prediction_to_genus_m<score>.csv (and a full-list variant).
Column names drift between releases, so matching is substring-based:
  id    = first column (Virus/query)
  host  = column containing 'genus', else 'host'
  score = column containing 'score'/'confidence', else the 3rd column
One row per query is emitted; queries absent from the CSV are simply not listed
(the master join fills those with NA).
"""
import argparse
import glob
import os
import sys

import pandas as pd


def pick(path):
    df = pd.read_csv(path)
    if df.empty:
        return None
    cols = list(df.columns)
    idc = cols[0]
    h = next((c for c in cols if "genus" in c.lower()), None)
    h = h or next((c for c in cols if "host" in c.lower()), cols[min(1, len(cols) - 1)])
    s = next((c for c in cols if "score" in c.lower() or "confidence" in c.lower()), None)
    s = s or cols[min(2, len(cols) - 1)]
    out = pd.DataFrame({"id": df[idc], "host": df[h], "host_score": df[s]})
    # iPHoP can emit several confident hosts per phage; keep the top-confidence call
    # so the master join stays one-row-per-id (the full multi-call CSV stays in out_dir).
    out["host_score"] = pd.to_numeric(out["host_score"], errors="coerce")
    out = out.sort_values("host_score", ascending=False).drop_duplicates("id", keep="first")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="iPHoP out_dir")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # prefer the genus-level filtered list; fall back to any host CSV
    cands = (
        sorted(glob.glob(os.path.join(a.indir, "Host_prediction_to_genus_*.csv")))
        + sorted(glob.glob(os.path.join(a.indir, "Host_prediction_*.csv")))
        + sorted(glob.glob(os.path.join(a.indir, "*host*.csv")))
    )
    cands = list(dict.fromkeys(cands))  # order-preserving dedup

    wrote = False
    with open(a.out, "w", newline="") as fh:
        fh.write("id\thost\thost_score\n")
        for p in cands:
            try:
                out = pick(p)
            except Exception as e:  # noqa: BLE001 - try next candidate
                print(f"[iphop_normalize] skip {p}: {e}", file=sys.stderr)
                continue
            if out is None or out.empty:
                continue
            out.to_csv(fh, sep="\t", index=False, header=False)
            wrote = True
            print(f"[iphop_normalize] {os.path.basename(p)} -> {len(out)} hosts", file=sys.stderr)
            break
    if not wrote:
        print("[iphop_normalize] no host predictions found; wrote header only", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
