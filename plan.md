# Phage Genome Downstream Analysis Pipeline — Consolidated Plan

> Snakemake workflow. Per-tool conda envs. Takes `phage_genomes/*.fna` (24 samples, ~5,040 contigs) and produces a quality-filtered, dereplicated, annotated, taxonomically classified phage catalog with **per-sample comparison outputs** (vOTU and gene count tables + ecological metrics).

## 1. Inputs

* `phage_genomes/KennethNNN_mobile.fna` — multi-FASTA per sample (24 files; Kenneth012 absent, Kenneth018 empty).
* Each contig carries its source sample through every step (ID = `{sample}__{original_id}`).
* **Reads (for P10):** `/media/yasas-wijesekara/Bender/ATLAS_mouse/mouse/QC/reads/{sample}_R{1,2}.fastq.gz`, already QC'd. Sample IDs match the genome files exactly (Kenneth004–028, minus 012; Kenneth018 reads are ~8 MB → empty, skipped).

## 2. Output contract (what users compare)

| Output | Purpose |
|---|---|
| `tables/votu_sample_counts.tsv` | vOTU × sample matrix (counts, and abundance if reads supplied) |
| `tables/gene_sample_counts.tsv` | gene/protein family × sample matrix (PHROG/KO/VOG) |
| `tables/votu_master.tsv` | one row per vOTU: rep, length, CheckV quality, taxonomy, lifestyle, host, AMG flag |
| `tables/sample_summary.tsv` | one row per sample: n contigs, n vOTUs, n genes, completeness stats |
| `ecology/alpha_diversity.tsv` | richness, Shannon, Simpson, Chao1 per sample |
| `ecology/beta_diversity.tsv` + PCoA/NMDS coords | Bray–Curtis / Jaccard distances, ordination |
| `report/phage_report.html` | integrated, browsable summary |

Count tables are wide format (rows = features, columns = samples) for direct import into R/phyloseq/vegan.

## 3. Phases

### P1 — Preprocess & sample tracking
* Concatenate FASTAs, prefix IDs with sample, length-filter (≥5 kb default, configurable), exact-sequence dedup (`seqkit rmdup`), emit a `contig → sample` map used by every downstream table.
* Warn + skip empty files (Kenneth018).

### P2 — QC
* **CheckV** `end_to_end`: completeness, contamination, proviral-region detection, MIUViG quality tiers.

### P3 — Dereplication (vOTUs)
* **vclust** `prefilter → align → cluster` at 95% ANI / 85% alignment fraction (ICTV/MIUViG species-level vOTU).
* Representative selection: CheckV `Complete` > `High-quality` > longest (not longest-only, to avoid picking a contaminated long contig).
* Outputs `votu_clusters.tsv` (contig → vOTU) and `votu_representatives.fna`.

### P4 — Gene calling & functional annotation
* **Prodigal** (or **PHANOTATE**, configurable) → proteins + GFF.
* **PHROG** search via **MMseqs2** (or **Pharokka** which bundles Prodigal + PHROG + tRNA/tmRNA + terminase). Choose one; Pharokka is the turn-key option.
* Gene-level table feeds the `gene_sample_counts.tsv` matrix (by PHROG/KO/VOG family).

### P5 — Taxonomy
* **vContact3** on representative proteins → viral clusters + ICTV-style taxonomy (Family/Genus where assignable; "Unclassified" otherwise).

### P6 — Hallmark-gene phylogeny (optional)
* Extract terminase large subunit (TerL) from annotations, align (MAFFT), tree (IQ-TREE2). Per-vOTU where present.

### P7 — AMG detection (optional)
* **DRAM-v** distill → AMG table + pathway. (Heavy; opt-in.)

### P8 — Lifestyle prediction
* **ProkBERT PhaStyle** (Juhász et al., 2025, *Bioinform. Adv.* vbaf188): nucleotide-language-model classifier, accurate on fragmented contigs (500 bp–10 kb), fast, no database, outperforms BACPHLIP/DeePhage/PhaTYP on out-of-sample phages.
* Supplied as Docker/Apptainer image (`obalasz/phage`); wrap via `--use-singularity` or a dedicated conda env if a pip package exists.
* Cross-check with integrase/repressor markers from P4 for confidence flagging.

### P9 — Host prediction (optional)
* **iPHoP** (`simroux/iphop`) — integrated k-mer + CRISPR + protein signals; ships a built-in `iphop download` command, so `setup_databases.smk` invokes it into `db_dir` rather than manual download.
* Note: full iPHoP DB is large (~tens of GB); a smaller split DB is available. Document the choice.
* Fallback **WIsH**/**RaFAH** if iPHoP is too heavy.

### P10 — Abundance (optional, requires reads)
* **CoverM** handles mapping + coverage end-to-end; no separate aligner rule. It wraps a configurable mapper (default `minimap2-sr`; also `bwa-mem`, `bwa-mem2`, `strobealign`) and emits TPM/RPKM/counts for vOTU reps → merged into `votu_sample_counts.tsv`.
* Rationale for minimap2 default: CoverM's native mapper, fast, accurate on short reads, one tool instead of two. Mapper is a one-line config change if a study prefers BWA-MEM2.
* If no reads: `votu_sample_counts.tsv` holds per-sample contig counts / presence-absence instead.

### P11 — Count matrices
* Build `votu_sample_counts.tsv` (from P3 membership + P10 abundance) and `gene_sample_counts.tsv` (from P4 family assignments). Wide format, samples as columns.

### P12 — Ecological metrics
* Alpha diversity (richness, Shannon, Simpson, Chao1) and beta diversity (Bray–Curtis, Jaccard) via `vegan`; ordination (PCoA/NMDS).
* Differential abundance (ANCOM-BC) and ordination coloring run **only if** `samples.tsv` provides a grouping column (set via `ecology.group_col`). Otherwise skipped — grouping is study-specific.

### P13 — Master report
* Single HTML (Snakemake `report` directive) linking the master table, count matrices, ecology plots, and per-genome annotations.

## 4. Workflow layout

```
workflow/
├── Snakefile
├── config/{config.yaml, samples.tsv}
├── envs/{checkv,vclust,vcontact3,pharokka,dramv,iphop,prodigal,mmseqs2,iqtree,seqkit,coverm,ecology,phastyle}.yaml
├── rules/{preprocess,qc,cluster,annotate,taxonomy,phylo,amg,lifestyle,host,abundance,matrices,ecology,report}.smk
└── scripts/{parse_checkv,select_representatives,build_count_tables,ecology.R,report}.py
```

* One conda env per tool (`--use-conda`); optional `--use-singularity` for full containerization (needed for PhaStyle).
* Heavy optional phases (P6/P7/P9/P10) toggled in `config.yaml`.
* Databases: one-time `setup_databases.smk` downloads CheckV/vclust/vContact3/PHROG/Pharokka/iPHoP DBs to a configurable `db_dir`.

## 5. config.yaml (sketch)

```yaml
samples: config/samples.tsv
outdir: results
db_dir: databases

preprocess: { min_length: 5000, dedup: true }
qc:         { threads: 16 }
cluster:    { ani: 95, align_frac: 0.85, threads: 32 }
annotate:   { gene_caller: prodigal, phrog: true, pharokka: true }
taxonomy:   { tool: vcontact3, threads: 16 }
lifestyle:  { tool: phastyle }          # ProkBERT PhaStyle (Docker/Singularity)

# opt-in heavy modules
phylo:      { enabled: false, threads: 16 }
amg:        { enabled: false, tool: dramv }
host:       { enabled: false, tool: iphop }
abundance:  { enabled: true, mapper: minimap2-sr, coverage: coverm }   # mapper: minimap2-sr | bwa-mem | bwa-mem2 | strobealign

ecology:
  group_col: null   # column in samples.tsv for differential abundance (optional)
```

## 6. Remaining decisions

1. **PhaStyle packaging:** Singularity/Apptainer image `obalasz/phage` via `--use-singularity` (recommended; no conda package exists). OK?
2. **iPHoP DB size:** default to the smaller split DB; full DB available on request (only if P9 enabled).
3. **KEGG/KOfam access** for AMG pathway completeness (P7), or skip (default: skip)?
4. **Representative selection:** quality-first (proposed) vs longest-only — OK as proposed?
5. **Grouping columns:** provide the `samples.tsv` column name(s) per study when needed; pipeline stays grouping-agnostic until then.

## 7. Next steps after approval

1. Finalize tool list + toggles; confirm P10 (reads) availability.
2. Write `config/samples.tsv` from `phage_genomes/`, `config/config.yaml`.
3. Scaffold Snakefile + per-phase rule files + `envs/*.yaml`.
4. Implement helper scripts (parse_checkv, select_representatives, build_count_tables, ecology.R).
5. Smoke-test on Kenneth004 + Kenneth005 (non-empty, small).
6. Add `setup_databases.smk` and a run README.

## References
* CheckV — Nayfach et al. 2021, *Nat Biotechnol* 39:578.
* vclust — Zieleziński et al. 2025, *Nat Methods* (refresh-bio/vclust).
* vContact3 — vcontact3.readthedocs.io.
* PHROG — Terzian et al. 2021, *NAR*.
* Pharokka — Bouras et al. 2023, *Bioinformatics*.
* ProkBERT PhaStyle — Juhász et al. 2025, *Bioinform. Adv.* 5:vbaf188.
* DRAM-v — Shaffer et al. 2020, *NAR Genom Bioinform*.
* iPHoP — Roux et al. 2023, *Nat Biotechnol*.
