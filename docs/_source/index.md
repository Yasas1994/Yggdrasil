
Yggdrasil : a reproducible Snakemake pipeline for downstream analysis of metagenomic phage genomes
===============

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License">
  </a>
  <a href="https://snakemake.github.io">
    <img src="https://img.shields.io/badge/snakemake-%E2%89%A88-149f8f" alt="Snakemake">
  </a>
  <a href="https://docs.python.org/3/">
    <img src="https://img.shields.io/badge/python-%E2%89%A53.10-3776ab" alt="Python">
  </a>
</p>

Yggdrasil takes per-sample phage genomes predicted from metagenomes (and
optionally the QC'd reads) and produces a quality-filtered, dereplicated,
annotated and taxonomically classified phage catalog with **per-sample
comparison outputs**: vOTU × sample and gene × sample count tables, ecological
metrics, and a single master vOTU table.

```{toctree}
:maxdepth 1
:caption contents

installation
usage
pipeline
outputs
configuration
```

## At a glance

```bash
pip install .
yggdrasil init -g phage_genomes/ -r QC/reads/ -o run1
yggdrasil setup-databases -w run1 --cores 8
yggdrasil run -w run1 --cores 32 --use-singularity
```

Every external tool runs in its own conda environment (env isolation); the
pipeline itself is a pip-installable package that wraps a bundled Snakemake
workflow.
