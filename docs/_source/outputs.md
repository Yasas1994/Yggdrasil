# Outputs

All paths are relative to `<outdir>` (default `results/`).

---

## Per-sample comparison outputs

These are the tables most users want. All are wide (rows = features, columns =
samples) and import directly into R / `phyloseq` / `vegan`.

### `07_matrices/votu_sample_counts.tsv`

vOTU × sample contig counts.

```text
votu	Kenneth004	Kenneth005	Kenneth006
0	1	1	0
1	1	1	0
```

### `07_matrices/gene_sample_counts.tsv`

Gene-family (PHROG / KO / VOG) × sample counts.

```text
family	Kenneth004	Kenneth005
PHROG_0001	12	8
PHROG_0002	0	3
```

### `08_ecology/alpha.tsv`

Per-sample diversity.

| Column | Description |
|--------|-------------|
| `sample` | sample id |
| `richness` | number of observed vOTUs |
| `shannon` | Shannon index |
| `simpson` | Simpson index |
| `chao1` | Chao1 richness estimate |

### `08_ecology/beta.tsv`

Pairwise Bray–Curtis (and/or Jaccard) distance matrix, samples × samples.

### `08_ecology/ordination.png`

PCoA / NMDS ordination plot, coloured by `ecology.group_col` when set.

### `08_ecology/differential.tsv`

Only when `ecology.group_col` is set. Columns: `votu`, `method`, `pvalue`,
`padj`, `log2fc`, `mean_<group>`…. Uses ANCOM-BC if installed, else a base-R
fallback.

---

## Master tables

### `09_report/votu_master.tsv`

One row per vOTU with everything joined on the representative contig:

| Column | Source |
|--------|--------|
| `votu` | vclust cluster id |
| `representative` | representative contig id |
| `n_members` | contigs in the cluster |
| `n_samples` | samples the cluster appears in |
| `rep_length` | representative length (bp) |
| `checkv_quality`, `completeness`, `contamination` | CheckV |
| `taxonomy` | vContact3 |
| `lifestyle`, `lifestyle_score` | ProkBERT PhaStyle |
| `amg_flag` | DRAM-v (optional; NA if off) |
| `host`, `host_score` | iPHoP (optional; NA if off) |

### `09_report/phage_report.html`

A branded HTML report generated with
[bioinformatics-report-builder](https://github.com/Yasas1994/bioinformatics-report-builder)
and rendered by Quarto. It includes the vOTU master table, the vContact3
gene-sharing network (sample vOTUs highlighted and coloured by taxonomy), alpha
diversity, and the vOTU/gene count matrices. When the optional modules are
enabled, VIRIDIC intergenomic-similarity and VirClust hierarchical-clustering
sections (figure + cluster table) are appended. Linked figures and assets are
written to `09_report/assets/`. Open `phage_report.html` in a browser.

Requires **Quarto >=1.4** installed on the host.

---

## Per-step outputs

| Path | Description |
|------|-------------|
| `00_preprocess/all_contigs.fna` | concatenated, renamed, length-filtered contigs |
| `00_preprocess/contig2sample.tsv` | contig → sample map |
| `00_preprocess/lengths.tsv` | contig → length |
| `01_qc/checkv/quality_summary.tsv` | CheckV completeness / contamination / provirus |
| `02_cluster/votu_clusters.tsv` | contig → vOTU → `is_representative` |
| `02_cluster/votu_representatives.fna` | one representative per vOTU |
| `02_cluster/vclust/` | raw vclust ANI + cluster files |
| `03_annotate/proteins.faa` | called proteins (representatives) |
| `03_annotate/genes.gff` | gene calls |
| `03_annotate/gene_families.tsv` | contig, gene_id, family (PHROG/KO/VOG) |
| `04_taxonomy/taxonomy.tsv` | id, taxonomy (per representative) |
| `04_taxonomy/vcontact3/` | raw vContact3 network + profiles |
| `05_lifestyle/lifestyle.tsv` | id, lifestyle, lifestyle_score |
| `06_abundance/coverm.tsv` | representative × sample TPM (if reads) |
| `opt_phylo/terl.treefile` | TerL tree (optional) |
| `opt_amg/amg_summary.tsv` | id, amg_flag + DRAM-v columns (optional) |
| `opt_host/iphop.tsv` | id, host, host_score (optional) |
| `opt_viridic/clusters.csv` | VIRIDIC genome → species/genus cluster (optional) |
| `opt_viridic/similarity_matrix.tsv` | VIRIDIC pairwise intergenomic similarity matrix (optional) |
| `opt_viridic/heatmap.pdf` | VIRIDIC heatmap (optional) |
| `opt_virclust/tree.newick` | VirClust hierarchical tree (optional) |
| `opt_virclust/distance_matrix.tsv` | VirClust protein-content intergenomic distance matrix (optional) |
| `opt_virclust/genome_clusters.tsv` | VirClust genome → viral genome cluster (VGC) + sharing stats (optional) |
| `opt_virclust/core_proteins.faa` | VirClust core proteins across VGCs (optional) |
