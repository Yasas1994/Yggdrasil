# AGENTS.md — Yggdrasil (phage_pipe)

Reproducible Snakemake pipeline for downstream analysis of metagenomic phage
genomes. Pip-installable Python package (`yggdrasil`) that bundles a Snakemake
workflow; every external bioinformatics tool runs in its own conda env (or a
Singularity/Apptainer image where noted).

Input: per-sample phage multi-FASTAs (`phage_genomes/*.fna`, ~24 "KennethNNN"
samples) + optional QC'd paired reads. Output: quality-filtered, dereplicated
(vOTU), annotated, taxonomically classified phage catalog with per-sample count
matrices, ecological metrics, and a master vOTU table.

## Commands

```bash
pip install -e ".[dev]"     # install package (editable) + pytest
pytest                      # run tests (pure-python smoke tests, no external tools)
yggdrasil init -g phage_genomes/ -r QC/reads/ -o run1   # scaffold run dir
yggdrasil setup-databases -w run1 --cores 8             # one-time DB download
yggdrasil run -w run1 --cores 32 --use-singularity      # full run (local)
yggdrasil run -w run1 --executor slurm --jobs 200 --partition batch --use-singularity  # SLURM cluster
yggdrasil run -w run1 --executor slurm-jobstep --jobs 20 --use-singularity  # inside an existing SLURM allocation
yggdrasil run -w run1 --cores 32 -n                     # dry-run: validate the DAG
yggdrasil config                                          # print default config
# extra snakemake args after `--`: yggdrasil run -w run1 -- --rerun-incomplete
```

Docs (Sphinx, MyST markdown) build with `make -C docs html`
(deps in `docs/requirements.txt`).

## Layout

- `src/yggdrasil/cli.py` — CLI (`init`, `run`, `setup-databases`, `config`).
  `discover()` pairs genomes with reads by sample prefix
  (`Kenneth004_mobile.fna` ↔ `Kenneth004_R{1,2}.fastq.gz`).
- `src/yggdrasil/runner.py` — locates the bundled workflow and invokes
  Snakemake; deep-merges bundled default config first, user config second
  (user wins). `--executor slurm` adds `--executor slurm --jobs N
  --default-resources slurm_partition=<p>`; local mode passes `--cores N`.
- `src/yggdrasil/workflow/` — the bundled Snakemake workflow (force-included
  into the wheel via pyproject.toml):
  - `Snakefile` — shared helpers (`env()`, `scr()`, `flag()`, sample accessors)
    in the global namespace; included `.smk` files rely on them. Targets are
    assembled conditionally from config flags.
  - `rules/*.smk` — core phases: preprocess → qc (CheckV) → cluster (vclust,
    95% ANI / 85% cov) → annotate (Pharokka/PHROG) → taxonomy (vContact3) →
    lifestyle (PhaStyle) → abundance (CoverM) → matrices → ecology → report.
  - `rules/optional/` — opt-in heavy modules toggled by `enabled:` in config:
    phylo (TerL tree), amg (Pharokka/phold), host (iPHoP), viridic, virclust.
  - `envs/*.yaml` — one conda env per tool; rules reference them via `env()`.
  - `scripts/` — pure-Python/R helpers invoked by rules (`preprocess.py`,
    `select_representatives.py`, `build_count_tables.py`, `votu_master.py`,
    `ecology.R`, `build_report.py`, …). These carry the load-bearing logic.
  - `config/config.yaml` — bundled defaults; `yggdrasil init` copies it into
    the run directory.
- `tests/test_smoke.py` — pytest smoke tests exercising `discover()` and the
  pure-python spine scripts end-to-end on synthetic fixtures.
- `plan.md` — original consolidated design plan (phases P1–P13, tool choices,
  rationale). Source of truth for design intent.
- `docs/_source/*.md` — Sphinx user docs; keep in sync with behavior changes.
- `run_all/`, `run_smoke/` — local run directories with configs, logs, results
  (gitignored patterns: `run*/`, `results/`, `*.log`, `databases/`).
- `notebooks/` — downstream R/Jupyter analysis (MaAsLin3, vOTU selection);
  not part of the pipeline.
- `phage_genomes/`, `databases/` — input FASTAs and downloaded tool DBs
  (large; `databases/` is gitignored).

## Conventions

- Python ≥ 3.10, src layout, hatchling build. Runtime deps: snakemake ≥ 8,
  pandas, pyyaml, biopython. Heavy tools come from conda envs, never pip.
- Contig identity: every contig ID is `{sample}__{original_id}` and carries
  its source sample through all steps; `contig → sample` maps feed all tables.
- Count tables are wide format (features × samples) for R/vegan/phyloseq.
- vOTU representative selection is quality-first (CheckV Complete >
  High-quality > longest), not longest-only.
- `preprocess.dedup` defaults to `false` — keep it off for ecological counts.
- Ecology input switches automatically: CoverM TPM when `abundance.enabled` and
  at least one sample has reads, otherwise the vOTU × sample presence/absence
  matrix (no CoverM jobs run).
- PhaStyle runs via Singularity image (`--use-singularity`); phold AMG option
  wants a GPU.
- `ecology.group_col: null` (default) skips differential abundance; set it to
  a `samples.tsv` column to enable ANCOM-BC.

### Parallelization patterns (cluster scale)

- **Per-sample pattern** (qc.smk): preprocess.py `--out-per-sample` writes one
  renamed, length-filtered FASTA per sample (empty file for empty/missing
  genomes), so `checkv_sample` has a static parse-time input per sample and
  `checkv_merge` concatenates into the same output paths the single rule uses.
  TSV merge = `awk 'NR==1 || FNR>1'` (first header kept), FASTA merge = `cat`.
- **Chunking pattern** (annotate/lifestyle/host): a `checkpoint split_reps_*`
  rule runs scripts/split_fasta.py (`--seqs-per-chunk K`, chunk count =
  ceil(n_seqs/K)) writing `{prefix}_chunkNNN.fna` + `chunks.txt` into a split
  dir; per-chunk rules consume `reps_{chunk}.fna`; an input function reads
  `chunks.txt` via `checkpoints.<name>.get().output[0]` so the DAG re-evaluates
  after the checkpoint. A merge rule then rebuilds the exact output paths the
  single-run rules produce, so downstream rules are layout-agnostic.
- Both patterns are config-gated at module level
  (`if bool(config["qc"].get("per_sample", True)): ... else: ...`,
  `if int(config[...].get("seqs_per_chunk", 2000)) > 0:`) so exactly one rule
  defines each output. `0` / `false` selects the original single-run rules.
- New heavy rules must declare `resources: mem_mb=..., runtime=...` (minutes)
  so `--executor slurm` submissions request sensible allocations; threads map
  to cpus via the SLURM executor plugin.

## Working rules

- Validate workflow changes with `pytest` plus
  `yggdrasil run -w run_smoke -n` (dry-run the DAG) before a real run.
  Checkpointed rules show as "Defined but x jobs will only be evaluated after
  checkpoint" in dry-runs — expected; check for WorkflowError and verify the
  non-checkpointed parts of the DAG.
- Helper scripts must stay runnable standalone (tests invoke them via
  `subprocess` with CLI args) and tool-free — external tools belong in rules,
  not scripts.
- New pipeline phase: add `rules/<phase>.smk` (use `env()`/`scr()`/`flag()`
  helpers), an `envs/<tool>.yaml`, an `enabled:`-style config toggle if heavy,
  and wire the target into `Snakefile` TARGETS.
- Update `README.md` (outputs table), `docs/_source/`, and this file when
  outputs, config keys, or commands change.
- Never commit `databases/`, `results/`, `run*/`, or `*.log` (see `.gitignore`).
