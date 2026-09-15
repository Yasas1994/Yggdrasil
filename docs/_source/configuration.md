# Configuration

`yggdrasil init` copies the bundled defaults into `<workdir>/config.yaml`. The
runner deep-merges your file over the defaults, so you only need to specify the
keys you change.

```yaml
outdir: results
threads: 16
samples: samples.tsv

databases:
  dir: databases

preprocess:
  min_length: 5000
  dedup: false                # exact-sequence dedup; keep OFF for ecological counts

qc:
  threads: 16                 # single-run CheckV (per_sample: false)
  per_sample: true            # one CheckV job per sample (scales with sample count)
  sample_threads: 4           # threads per per-sample CheckV job

cluster:
  ani: 95                     # species-level vOTU ANI (percent)
  align_frac: 0.85            # minimum alignment fraction (0–1)
  threads: 32

annotate:
  phrog: true
  pharokka: true              # Pharokka (prodigal-gv + PHROG + tRNA); false -> prodigal-gv + MMseqs2/PHROG
  threads: 16
  seqs_per_chunk: 2000        # split vOTU reps into chunks of N seqs, one Pharokka job each; 0 = single run
  chunk_threads: 8            # threads per Pharokka chunk job

taxonomy:
  threads: 16
  db_version: 230             # vContact3 DB version (from prepare_databases)

lifestyle:
  tool: phastyle
  threads: 8
  seqs_per_chunk: 2000        # chunk reps for parallel PhaStyle/BACPHLIP jobs; 0 = single run
  image: docker://obalasz/phastyle:latest   # or a local .sif path

abundance:
  enabled: true
  mapper: minimap2-sr         # minimap2-sr | bwa-mem | bwa-mem2 | strobealign
  threads: 16

ecology:
  group_col: null             # samples.tsv column; null -> skip differential

report:
  enabled: true

# ---- opt-in heavy modules (off by default) ----
phylo:
  enabled: false
  threads: 16
amg:
  enabled: false
  threads: 16
host:
  enabled: false
  db: iphop                  # DB directory name under databases.dir (iPHoP unpacks into a versioned subdir)
  seqs_per_chunk: 2000       # chunk reps for parallel iPHoP jobs; 0 = single run
virclust:
  enabled: false             # protein-based hierarchical clustering of vOTUs
  threads: 16
  gene_code: 11
viridic:
  enabled: false             # ICTV-style intergenomic similarities + clustering
  threads: 16
  species_threshold: 95      # intergenomic similarity (%) for species clusters
  genus_threshold: 70        # intergenomic similarity (%) for genus clusters
```

---

## Reference

| Key | Type | Description |
|-----|------|-------------|
| `outdir` | str | output root, relative to the workdir |
| `samples` | str | path to `samples.tsv` |
| `databases.dir` | str | root for all tool databases |
| `preprocess.min_length` | int | drop contigs shorter than this (bp) |
| `preprocess.dedup` | bool | global exact-sequence dedup (off recommended) |
| `qc.per_sample` | bool | one CheckV job per sample + merge (true), or a single CheckV job over all contigs (false) |
| `qc.sample_threads` | int | threads per per-sample CheckV job |
| `cluster.ani` | int | vOTU ANI threshold, percent (→ 0–1 for vclust) |
| `cluster.align_frac` | float | minimum alignment fraction (0–1) |
| `annotate.pharokka` | bool | use Pharokka (prodigal-gv + PHROG); false uses prodigal-gv + MMseqs2/PHROG |
| `annotate.seqs_per_chunk` | int | vOTU reps per Pharokka chunk job (chunk count = ceil(n_reps / N)); 0 = single run |
| `annotate.chunk_threads` | int | threads per Pharokka chunk job |
| `taxonomy.db_version` | int | vContact3 prepared-DB version |
| `lifestyle.seqs_per_chunk` | int | reps per PhaStyle/BACPHLIP chunk job; 0 = single run |
| `lifestyle.image` | str | PhaStyle container: docker URI or local `.sif` |
| `abundance.enabled` | bool | map reads and quantify (needs reads) |
| `abundance.mapper` | str | CoverM mapper backend |
| `ecology.group_col` | str/null | samples.tsv column for differential + colouring |
| `phylo.enabled` | bool | TerL phylogeny |
| `amg.enabled` | bool | DRAM-v AMG detection |
| `host.enabled` | bool | iPHoP host prediction |
| `host.db` | str | DB directory name under `databases.dir` (iPHoP unpacks into a versioned subdir) |
| `host.seqs_per_chunk` | int | reps per iPHoP chunk job; 0 = single run |
| `virclust.enabled` | bool | VirClust protein-based hierarchical clustering + core proteins |
| `virclust.gene_code` | int | NCBI translation table for MetaGeneAnnotator (11 = bacteria/archaea) |
| `viridic.enabled` | bool | VIRIDIC intergenomic similarities + species/genus clustering |
| `viridic.species_threshold` | int | species-level intergenomic similarity cut-off (%) |
| `viridic.genus_threshold` | int | genus-level intergenomic similarity cut-off (%) |
