# Installation

Yggdrasil is a Python package that wraps a Snakemake workflow. The heavy
bioinformatics tools are **not** pip dependencies — each runs in its own conda
environment built by Snakemake on first use.

---

## Requirements

| Dependency | Why | Install |
|------------|-----|---------|
| Python ≥ 3.10 | the `yggdrasil` package | any |
| Snakemake ≥ 8 | workflow engine | `pip` (installed with the package) |
| conda / mamba | per-tool env isolation | miniforge / miniconda |
| Singularity or Apptainer | ProkBERT PhaStyle container | system package |
| R + vegan | ecology metrics | pulled via the `ecology` conda env |
| Quarto ≥ 1.4 | HTML report rendering | system package (required only if `report.enabled: true`) |

---

## pip (recommended)

```bash
git clone <repo-url> phage_pipe && cd phage_pipe
pip install .            # or: pip install -e .  (editable / development)
yggdrasil --version
```

This installs the `yggdrasil` CLI plus `snakemake`, `pandas`, `pyyaml` and
`biopython`. The analysis tools (CheckV, vclust, vContact3, Pharokka, CoverM,
…) are fetched into isolated conda envs the first time you `yggdrasil run`.

> **Tip:** pre-build the envs without running the pipeline:
> ```bash
> yggdrasil run -w run1 --cores 8 -- --conda-create-envs-only
> ```
> The `--` forwards everything after it straight to Snakemake.

---

## Databases

The core workflow needs three reference databases; the optional modules add a
few more. Download them once with:

```bash
yggdrasil setup-databases -w run1 --cores 8
```

| Database | Used by | Approx. size | Default location |
|----------|---------|--------------|------------------|
| CheckV   | QC      | ~6 GB        | `databases/checkv` |
| Pharokka / PHROG | annotation | ~2 GB | `databases/phrog` |
| vContact3 | taxonomy | ~7 GB | `databases/vcontact3` |
| iPHoP    | host prediction (optional) | ~300 GB (latest_rw full); ~9 GB test DB (iPHoP_db_rw_1.4_for-test) | `databases/iphop` |
| DRAM     | AMGs (optional) | ~10–20 GB | `databases/dram` |
| PhaStyle image | lifestyle | ~9 GB | `databases/phastyle.sif` |
| VirClust scripts | hierarchical clustering (optional) | ~1 MB (git clone) | `databases/virclust/repo` |
| VIRIDIC scripts | intergenomic similarity (optional) | ~1 MB (git clone) | `databases/viridic/repo` |

`setup-databases` runs each tool's own downloader inside its conda env
(`checkv download_database`, `install_databases.py`,
`vcontact3 prepare_databases`, `iphop download`) and clones the VirClust /
VIRIDIC standalone R scripts. Some steps print manual instructions when a
downloader is not available — follow them, then re-run.

The VirClust and VIRIDIC conda envs (`envs/virclust.yaml`, `envs/viridic.yaml`)
are heavy (R + BLAST + HMMER/HH-suite/MCL for VirClust) and build on first run
only when the module is enabled. VirClust runs branch A (PC-level clustering +
tree + core proteins); the HMM/PSC/PSSC and annotation branches are not wired
in (they need extra reference DBs).

Set `databases.dir` in `config.yaml` to share one copy across runs.

---

## Verify

```bash
yggdrasil --help
yggdrasil config                 # print the bundled default config
yggdrasil run -w run1 --dry-run  # validate the DAG without running tools
```
