# Pipeline

Yggdrasil is a directed acyclic graph of Snakemake rules. Every contig carries
its source sample through the whole graph, which is what makes the per-sample
count tables fall out at the end.

```
preprocess → CheckV → vclust (95% ANI / 85% cov) → representatives
   → prodigal-gv/Pharokka + PHROG → vContact3 → PhaStyle → CoverM
   → count matrices → ecology → master report
```

---

## Stages

| # | Stage | Tool | Core? |
|---|-------|------|-------|
| P1 | Preprocess & sample tracking | Python + Biopython | yes |
| P2 | Quality control | CheckV | yes |
| P3 | Dereplication (vOTUs) | vclust | yes |
|| P4 | Gene calling & functional annotation | prodigal-gv / Pharokka + PHROG | yes |
| P5 | Taxonomy | vContact3 | yes |
| P6 | Lifestyle (temperate / virulent) | ProkBERT PhaStyle | yes |
| P7 | Abundance | CoverM (minimap2) | yes (if reads) |
| P8 | Count matrices | Python (pandas) | yes |
| P9 | Ecological metrics | R + vegan | yes |
|| P10 | Master report | bioinformatics-report-builder + Quarto | yes |
| — | TerL phylogeny | MAFFT + IQ-TREE2 | optional |
| — | AMG detection | DRAM-v | optional |
| — | Host prediction | iPHoP | optional |
| — | Intergenomic similarity + species/genus clustering | VIRIDIC | optional |
| — | Protein-based hierarchical clustering + core proteins | VirClust | optional |

---

## Key decisions

### Dereplication threshold

vOTUs are defined by vclust at **95% average nucleotide identity over 85%
alignment fraction** (the MIUViG / ICTV species-level convention):

```yaml
cluster:
  ani: 95          # percent — converted to a 0–1 fraction for vclust
  align_frac: 0.85
```

### Representative selection

Each vOTU's representative is the contig with the **best CheckV quality tier**,
tie-broken by length — not the longest contig, which can be a contaminated
outlier:

`Complete` > `High-quality` > `Medium-quality` > `Low-quality` > longest.

If CheckV output is unavailable for a contig, selection falls back to length.

### Annotation on representatives

Genes are called and annotated on the **vOTU representatives only** (not every
contig). This keeps vContact3, Pharokka and PhaStyle tractable on thousands of
contigs while still characterising every cluster.

### Lifestyle from nucleotides

Lifestyle uses **ProkBERT PhaStyle**, a genomic language model that classifies
directly from nucleotide sequence and stays accurate on fragmented contigs
(500 bp – 10 kb). It runs in the `obalasz/phastyle` container with the
`neuralbioinfo/PhaStyle-mini` model (fetched from Hugging Face on first run).

---

## Tools & citations

| Tool | Reference |
|------|-----------|
| CheckV | Nayfach *et al.* 2021, *Nat Biotechnol* 39:578 |
| vclust | Zieleziński *et al.* 2025, *Nat Methods* |
| vContact3 | vcontact3.readthedocs.io |
| PHROG | Terzian *et al.* 2021, *NAR* |
| Pharokka | Bouras *et al.* 2023, *Bioinformatics* |
| ProkBERT PhaStyle | Juhász *et al.* 2025, *Bioinform. Adv.* 5:vbaf188 |
| DRAM-v | Shaffer *et al.* 2020, *NAR Genom Bioinform* |
| iPHoP | Roux *et al.* 2023, *Nat Biotechnol* |
| VIRIDIC | Moraru *et al.* 2020, *Viruses* 12(11):1268 |
| VirClust | Moraru 2023, *Viruses* 15(4):1007 |
| CoverM | wwood/CoverM |
| vegan | Oksanen *et al.*, R package |
