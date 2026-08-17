#!/usr/bin/env python3
"""Run BACPHLIP on a multi-FASTA of phage genomes and normalize its output.

Stdlib-only (runs inside the BACPHLIP container, which has python + hmmer but
not necessarily pandas). BACPHLIP writes one ``.bacphlip`` table per record in
multi-fasta mode; we glob every ``*.bacphlip`` file, coerce the Virulent/
Temperate probabilities, and emit a single clean TSV:

    id  bacphlip_virulent  bacphlip_temperate  bacphlip_lifestyle  bacphlip_confidence

lifestyle = argmax(virulent, temperate); confidence = max probability.
Column discovery is substring-based so minor header drift across builds is OK.
"""
import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile


def find_col(cols, *keys):
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


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def parse_bacphlip(path):
    out = []
    with open(path, newline="") as fh:
        rows = list(csv.DictReader((r for r in fh if not r.startswith("#")), delimiter="\t"))
    if not rows:
        return out
    cols = list(rows[0].keys())
    idc = find_col(cols, "contig", "sequence", "genome", "name", "id") or cols[0]
    vc = find_col(cols, "virulent", "lytic")
    tc = find_col(cols, "temperate", "lysog")
    pc = find_col(cols, "lifestyle", "prediction", "label", "class")
    for row in rows:
        v = fnum(row.get(vc)) if vc else float("nan")
        t = fnum(row.get(tc)) if tc else float("nan")
        life = norm_label(row.get(pc)) if pc else ""
        if not life:
            life = "virulent" if (v or 0) >= (t or 0) else "temperate"
        conf = max(x for x in (v, t) if x == x) if any(x == x for x in (v, t)) else float("nan")
        out.append((row.get(idc, ""), v, t, life, conf))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", default="1")
    a = ap.parse_args()

    header = ["id", "bacphlip_virulent", "bacphlip_temperate", "bacphlip_lifestyle", "bacphlip_confidence"]
    wd = tempfile.mkdtemp(prefix="bacphlip_")
    seen, rows = {}, []
    try:
        src = os.path.join(wd, "votu.fna")
        shutil.copy(a.fasta, src)
        subprocess.run(["bacphlip", "-i", src, "--multi_fasta", "-f"], cwd=wd, check=True)
        for p in glob.glob(os.path.join(wd, "**", "*.bacphlip"), recursive=True):
            for rec in parse_bacphlip(p):
                rid = rec[0]
                if rid and rid not in seen:
                    seen[rid] = rec
        rows = list(seen.values())
        if not rows:
            print("[bacphlip] WARNING: no .bacphlip output produced", file=sys.stderr)
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            for rid, v, t, life, conf in rows:
                w.writerow([rid, f"{v:.6g}" if v == v else "", f"{t:.6g}" if t == t else "", life,
                            f"{conf:.6g}" if conf == conf else ""])
        print(f"[bacphlip] genomes={len(rows)} out={a.out}", file=sys.stderr)
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
