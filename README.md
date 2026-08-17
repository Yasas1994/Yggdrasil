# Yggdrasil

Reproducible Snakemake pipeline for downstream analysis of metagenomic phage
genomes. Pip-installable; every external tool runs in its own conda environment.

Takes per-sample phage FASTAs (+ optional QC'd reads) and produces a
quality-filtered, dereplicated (vOTU), annotated, taxonomically classified phage
catalog with **per-sample comparison outputs**: vOTU × sample and gene × sample
count tables, ecological metrics, and a master vOTU table.

## Install

```bash
git clone <repo-url> phage_pipe && cd phage_pipe
pip install .                 # or: pip install -e .  (editable/dev)
# heavy tools come from conda envs, not pip:
conda install -c conda-forge -c bioconda snakemake   # if not already present
```

## Quickstart

```bash
# 1) scaffold a run directory (pairs genomes with reads by sample prefix)
yggdrasil init -g phage_genomes/ -r QC/reads/ -o run1

# 2) one-time database download (large; some steps may be manual — see output)
yggdrasil setup-databases -w run1 --cores 8

# 3) run (builds per-tool conda envs on first use)
yggdrasil run -w run1 --cores 32 --use-singularity
```

`init` writes `run1/samples.tsv` (sample, genome, r1, r2) and `run1/config.yaml`.
Edit `config.yaml` to toggle modules/thresholds, then `run`. Add a metadata column
to `samples.tsv` and set `ecology.group_col` to enable differential abundance.

Reads are paired to genomes by sample prefix, so `Kenneth004_mobile.fna` pairs with
`Kenneth004_R{1,2}.fastq.gz`.

## Outputs (under `<outdir>`, default `results/`)

| Path | Description |
|---|---|
| `01_qc/checkv/quality_summary.tsv` | CheckV completeness/contamination |
| `02_cluster/votu_clusters.tsv` | contig → vOTU → representative |
| `02_cluster/votu_representatives.fna` | one rep per vOTU |
| `03_annotate/{proteins.faa,genes.gff,gene_families.tsv}` | genes + PHROG/KO/VOG |
| `04_taxonomy/taxonomy.tsv` | vContact3 taxonomy per representative |
| `05_lifestyle/lifestyle.tsv` | ProkBERT PhaStyle temperate/virulent |
| `06_abundance/coverm.tsv` | CoverM abundance (if reads provided) |
| `07_matrices/votu_sample_counts.tsv` | vOTU × sample counts |
| `07_matrices/gene_sample_counts.tsv` | gene family × sample counts |
| `08_ecology/{alpha,beta}.tsv`, `ordination.png` | diversity + ordination |
| `09_report/votu_master.tsv` | one row per vOTU (everything joined) |
| `09_report/phage_report.html` | browsable summary |

Optional (off by default): `opt_phylo/` (TerL tree), `opt_amg/` (DRAM-v AMGs),
`opt_host/` (iPHoP hosts). Enable in `config.yaml`.

## Pipeline

```
preprocess → CheckV → vclust (95% ANI / 85% cov) → representatives
   → Prodigal/Pharokka + PHROG → vContact3 → PhaStyle → CoverM
   → count matrices → ecology → master report
```

## Tools & citations

CheckV (Nayfach 2021), vclust (Zieleziński 2025), vContact3, PHROG (Terzian 2021),
Pharokka (Bouras 2023), ProkBERT PhaStyle (Juhász 2025), DRAM-v (Shaffer 2020),
iPHoP (Roux 2023), CoverM, vegan.

## Dev

```bash
pip install -e ".[dev]"
pytest
yggdrasil run -w run1 --dry-run     # validate the DAG without running tools
```

## License

MIT.
