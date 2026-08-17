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


def run(workdir: Path, user_config: Path | None, extra: list[str]) -> int:
    """Run the workflow. Returns Snakemake's exit code.

    Snakemake deep-merges --configfile in order, so the bundled defaults are
    passed first and the user's config second (user wins). If no config is
    given, <workdir>/config.yaml is used when present.
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
    cmd += extra
    print("[yggdrasil]", " ".join(cmd), file=sys.stderr)
    return subprocess.call(cmd)


def setup_databases(workdir: Path, user_config: Path | None, extra: list[str]) -> int:
    return run(workdir, user_config, extra + ["setup_databases"])
