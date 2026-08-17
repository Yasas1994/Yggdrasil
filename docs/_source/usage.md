# Running Yggdrasil

This guide covers the end-to-end workflow, input layout, the CLI, and how to
turn on optional analyses.

For every tunable parameter see the [Configuration](configuration.md) page; for
what gets written see [Outputs](outputs.md).

---

## Table of contents

- [Quick start](#quick-start)
- [Input layout](#input-layout)
- [The run directory](#the-run-directory)
- [Executing the pipeline](#executing-the-pipeline)
- [Optional modules](#optional-modules)
- [Sample grouping and differential abundance](#sample-grouping-and-differential-abundance)
- [Command-line reference](#command-line-reference)

---

## Quick start

```bash
# 1) scaffold a run directory (pairs genomes with reads, writes samples.tsv + config.yaml)
yggdrasil init -g phage_genomes/ -r QC/reads/ -o run1

# 2) download the reference databases once (large; some steps may be manual)
yggdrasil setup-databases -w run1 --cores 8

# 3) run (builds per-tool conda envs on first use)
yggdrasil run -w run1 --cores 32 --use-singularity
```

Results land in `run1/results/` (or wherever `outdir` points).

---

## Input layout

Yggdrasil expects one multi-FASTA of phage contigs per sample and, optionally,
the matching QC'd paired reads:

```
phage_genomes/Kenneth004_mobile.fna
phage_genomes/Kenneth005_mobile.fna
…
QC/reads/Kenneth004_R1.fastq.gz
QC/reads/Kenneth004_R2.fastq.gz
…
```

`yggdrasil init` pairs a genome with its reads by **sample prefix**, so a
genome named `Kenneth004_mobile.fna` is correctly paired with
`Kenneth004_R{1,2}.fastq.gz`. The longest reads-key that is a prefix of the
genome stem wins.

- Genome FASTAs are discovered by extension: `.fna`, `.fa`, `.fasta` (and `.gz`).
- Reads are discovered by `_R1` / `_R2` (or `_1` / `_2`) in the filename.
- Samples without reads still run through every step except abundance.
- Empty genome files are skipped with a warning.

---

## The run directory

`yggdrasil init` writes two files into the run directory:

### `samples.tsv`

```text
sample	genome	r1	r2
Kenneth004	/abs/phage_genomes/Kenneth004_mobile.fna	/abs/QC/reads/Kenneth004_R1.fastq.gz	/abs/QC/reads/Kenneth004_R2.fastq.gz
Kenneth005	/abs/phage_genomes/Kenneth005_mobile.fna	/abs/QC/reads/Kenneth005_R1.fastq.gz	/abs/QC/reads/Kenneth005_R2.fastq.gz
```

Paths are absolute. Add extra columns (e.g. `treatment`) for
[grouping](#sample-grouping-and-differential-abundance).

### `config.yaml`

A copy of the bundled defaults. Edit it to change thresholds, threads, the
database location, and to toggle optional modules. The runner overlays your
`config.yaml` on the bundled defaults, so keys you omit keep their defaults.

---

## Executing the pipeline

```bash
yggdrasil run -w run1 --cores 32 --use-singularity
```

Common flags:

| Flag | Description | Default |
|------|-------------|---------|
| `-w, --workdir` | run directory containing `config.yaml` | `.` |
| `-c, --config` | path to `config.yaml` | `<workdir>/config.yaml` |
| `--cores` | max parallel jobs / threads | 1 |
| `--use-singularity` | run containerized steps (PhaStyle) | off |
| `-n, --dry-run` | build the DAG, run nothing | off |

Anything after a bare `--` is forwarded verbatim to Snakemake, e.g.:

```bash
# pre-build conda envs only
yggdrasil run -w run1 --cores 8 -- --conda-create-envs-only

# rerun a single rule
yggdrasil run -w run1 --cores 16 -- -R vclust --until select_representatives
```

> `--use-singularity` is required for the lifestyle step (ProkBERT PhaStyle
> ships as a container). Set `lifestyle.image` in `config.yaml` to a local
> `.sif` to avoid re-pulling:
> ```bash
> apptainer pull databases/phastyle.sif docker://obalasz/phastyle:latest
> ```

---

## Optional modules

Three heavier modules are **off by default**. Enable them in `config.yaml`:

```yaml
phylo: { enabled: true, threads: 16 }   # terminase (TerL) phylogeny (MAFFT + IQ-TREE2)
amg:   { enabled: true, threads: 16 }   # auxiliary metabolic genes (DRAM-v)
host:  { enabled: true, db: iphop }     # host prediction (iPHoP); db: iphop uses the latest downloaded DB; add -dbv iPHoP_db_rw_1.4_for-test for a tiny test DB
```

Each needs its database downloaded via `setup-databases`. AMG and host rules
fail with a clear message pointing at `setup-databases` if the DB is missing.

---

## Sample grouping and differential abundance

Add a metadata column to `samples.tsv` and point `ecology.group_col` at it:

```text
sample	genome	r1	r2	treatment
Kenneth004	…	…	…	control
Kenneth005	…	…	…	hfd
```

```yaml
ecology:
  group_col: treatment
```

When `group_col` is set, the ecology step also writes
`08_ecology/differential.tsv` (ANCOM-BC if installed, otherwise a base-R
Wilcoxon/Kruskal + BH fallback) and colours the ordination by group. Leave it
`null` to skip differential analysis.

---

## Command-line reference

```
yggdrasil [COMMAND] [OPTIONS]
```

| Command | Purpose |
|---------|---------|
| `yggdrasil init` | Scaffold a run directory (`samples.tsv` + `config.yaml`) |
| `yggdrasil run` | Run the workflow |
| `yggdrasil setup-databases` | Download tool databases once |
| `yggdrasil config` | Print the bundled default config |
| `yggdrasil --version` | Print the version |

### `yggdrasil init`

```
yggdrasil init -g GENOMES [-r READS] [-o WORKDIR]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-g, --genomes` | directory of per-sample phage FASTAs | required |
| `-r, --reads` | directory of QC'd reads (`{sample}_R{1,2}.fastq.gz`) | none |
| `-o, --outdir` | run / working directory | `.` |

### `yggdrasil run`

```
yggdrasil run [-w WORKDIR] [-c CONFIG] [--cores N] [--use-singularity] [-n] [-- SNAKEMAKE_ARGS...]
```

---

## Python integration

There is no public Python API. Drive the CLI via `subprocess`:

```python
import subprocess
subprocess.run(["yggdrasil", "run", "-w", "run1", "--cores", "32", "--use-singularity"], check=True)
```
