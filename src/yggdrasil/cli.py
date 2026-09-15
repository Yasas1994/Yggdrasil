"""Yggdrasil command-line interface."""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

from . import __version__, runner

GENOME_EXTS = (".fna", ".fa", ".fasta")
READ_EXTS = (".fastq.gz", ".fq.gz")


def _stem(p: Path) -> str:
    n = p.name
    for e in (".fastq.gz", ".fq.gz", ".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if n.endswith(e):
            return n[: -len(e)]
    return p.stem


def discover(genomes_dir: Path, reads_dir: Path | None) -> list[dict]:
    """Pair each genome with reads by sample prefix.

    Reads sample id = basename with trailing `_R[12]`/`.R[12]` removed.
    A genome is paired to the longest reads-key that is a prefix of its stem
    (handles stems like `Kenneth004_mobile` vs reads `Kenneth004`).
    """
    genomes = sorted(p for p in genomes_dir.iterdir() if p.name.endswith(GENOME_EXTS))
    reads: dict[str, dict[str, str]] = {}
    if reads_dir and reads_dir.exists():
        for p in sorted(reads_dir.iterdir()):
            if not p.name.endswith(READ_EXTS):
                continue
            s = _stem(p)
            m = re.search(r"[._]R?([12])$", s)
            if not m:
                continue
            sid = s[: m.start()]
            reads.setdefault(sid, {})[m.group(1)] = str(p.resolve())

    rows = []
    for g in genomes:
        gs = _stem(g)
        match = max((k for k in reads if gs == k or gs.startswith(k + "_")), key=len, default="")
        r = reads.get(match, {})
        rows.append({
            "sample": match or gs,
            "genome": str(g.resolve()),
            "r1": r.get("1", ""),
            "r2": r.get("2", ""),
        })
    return rows


def cmd_init(a) -> int:
    wd = Path(a.outdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    rows = discover(Path(a.genomes), Path(a.reads) if a.reads else None)
    if not rows:
        print(f"[yggdrasil] no genomes found in {a.genomes}", file=sys.stderr)
        return 1
    samples = wd / "samples.tsv"
    with samples.open("w") as fh:
        fh.write("sample\tgenome\tr1\tr2\n")
        for r in rows:
            fh.write(f"{r['sample']}\t{r['genome']}\t{r['r1']}\t{r['r2']}\n")
    cfg = wd / "config.yaml"
    if not cfg.exists():
        shutil.copy(runner.default_config(), cfg)
    paired = sum(1 for r in rows if r["r1"])
    print(f"[yggdrasil] wrote {samples} ({len(rows)} genomes, {paired} with reads)")
    print(f"[yggdrasil] wrote {cfg} — edit, then `yggdrasil run -w {wd}`")
    return 0


def cmd_run(a) -> int:
    extra = ["--use-conda"]
    if a.use_singularity:
        extra.append("--use-singularity")
    if a.dry_run:
        extra.append("-n")
    extra += a.extra
    return runner.run(Path(a.workdir), Path(a.config) if a.config else None, extra,
                      executor=a.executor, jobs=a.jobs, partition=a.partition,
                      cores=a.cores)


def cmd_setup(a) -> int:
    extra = ["--use-conda"] + a.extra
    return runner.setup_databases(Path(a.workdir), Path(a.config) if a.config else None, extra,
                                  cores=a.cores)


def cmd_config(a) -> int:
    print(Path(runner.default_config()).read_text())
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="yggdrasil", description="Downstream analysis of metagenomic phage genomes")
    p.add_argument("--version", action="version", version=f"yggdrasil {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("init", help="scaffold a run directory (samples.tsv + config.yaml)")
    pi.add_argument("-g", "--genomes", required=True, help="directory of per-sample phage FASTA files")
    pi.add_argument("-r", "--reads", help="directory of QC'd reads ({sample}_R{1,2}.fastq.gz)")
    pi.add_argument("-o", "--outdir", default=".", help="run/working directory (default: cwd)")
    pi.set_defaults(func=cmd_init)

    pr = sub.add_parser("run", help="run the workflow")
    pr.add_argument("-w", "--workdir", default=".", help="run directory containing config.yaml")
    pr.add_argument("-c", "--config", help="config.yaml (default: <workdir>/config.yaml)")
    pr.add_argument("--cores", type=int, default=1)
    pr.add_argument("--executor", choices=["local", "slurm"], default="local",
                    help="local: run on this machine with --cores; slurm: dispatch each rule as a SLURM job")
    pr.add_argument("--jobs", type=int, default=100, help="slurm mode: max concurrent jobs")
    pr.add_argument("--partition", default="batch", help="slurm mode: partition for submitted jobs")
    pr.add_argument("--use-singularity", action="store_true", help="also use Singularity/Apptainer containers (PhaStyle)")
    pr.add_argument("-n", "--dry-run", action="store_true")
    pr.add_argument("extra", nargs=argparse.REMAINDER, help="extra args forwarded to snakemake")
    pr.set_defaults(func=cmd_run)

    ps = sub.add_parser("setup-databases", help="download tool databases once")
    ps.add_argument("-w", "--workdir", default=".")
    ps.add_argument("-c", "--config")
    ps.add_argument("--cores", type=int, default=1)
    ps.add_argument("extra", nargs=argparse.REMAINDER)
    ps.set_defaults(func=cmd_setup)

    pc = sub.add_parser("config", help="print the bundled default config")
    pc.set_defaults(func=cmd_config)
    return p


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    fwd: list[str] = []
    if "--" in argv:
        i = argv.index("--")
        argv, fwd = argv[:i], argv[i + 1 :]
    a = build_parser().parse_args(argv)
    if hasattr(a, "extra"):
        a.extra = (a.extra or []) + fwd
    sys.exit(a.func(a))
