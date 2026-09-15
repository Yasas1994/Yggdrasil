"""Locate the bundled workflow and invoke Snakemake."""
from __future__ import annotations

import subprocess
import sys
from importlib.resources import files
from pathlib import Path

PKG = "yggdrasil"


def workflow_dir() -> Path:
    return Path(str(files(PKG) / "workflow"))


def snakefile() -> str:
    return str(workflow_dir() / "Snakefile")


def default_config() -> str:
    return str(workflow_dir() / "config" / "config.yaml")


def run(workdir: Path, user_config: Path | None, extra: list[str],
        executor: str = "local", jobs: int = 100, partition: str = "batch",
        cores: int | None = None) -> int:
    """Run the workflow. Returns Snakemake's exit code.

    Snakemake deep-merges --configfile in order, so the bundled defaults are
    passed first and the user's config second (user wins). If no config is
    given, <workdir>/config.yaml is used when present.

    executor=slurm dispatches each rule as its own SLURM job (needs
    snakemake-executor-plugin-slurm); executor=local runs locally with --cores.
    """
    workdir = workdir.resolve()
    if user_config is None:
        candidate = workdir / "config.yaml"
        user_config = candidate if candidate.exists() else None
    cmd = [
        "snakemake",
        "-s", snakefile(),
        "-d", str(workdir),
        "--configfile", default_config(),
    ]
    if user_config is not None:
        cmd += [str(user_config.resolve())]
    if executor == "slurm":
        cmd += ["--executor", "slurm", "--jobs", str(jobs),
                "--default-resources", f"slurm_partition={partition}"]
    elif cores is not None:
        cmd += ["--cores", str(cores)]
    cmd += extra
    print("[yggdrasil]", " ".join(cmd), file=sys.stderr)
    return subprocess.call(cmd)


def setup_databases(workdir: Path, user_config: Path | None, extra: list[str],
                    cores: int | None = None) -> int:
    return run(workdir, user_config, extra + ["setup_databases"], cores=cores)
