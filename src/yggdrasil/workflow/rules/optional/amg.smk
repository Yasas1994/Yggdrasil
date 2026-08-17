# Optional AMG (auxiliary metabolic gene) detection - lightweight, DRAM-v-free.
#
# DRAM-v (~150-200 GB DB) is the de facto AMG standard but is too heavy here. Instead,
# PHROG product/function annotations from Pharokka are mapped to curated, literature-
# grounded AMG metabolic classes by scripts/amg_classify.py. Optionally, phold (Bouras et al.
# 2026, NAR gkaf1448) can be run on top of Pharokka to re-annotate hypothetical proteins with
# ProstT5 + Foldseek structural homology (~16 GB DB, GPU-accelerated). When phold is enabled,
# its per-CDS predictions are merged with the Pharokka table before classification.
# Because phold is sensitive to long hypothetical proteins, it is opt-in via amg.phold.
#
# AMG metabolic classes: photosynthesis, carbon/PPP, nucleotide, phosphate, sulfur, nitrogen,
# CAZy, peptidoglycan, cofactors; PHROG "moron/AMG/host-takeover" as a catch-all.
#
# Outputs:
#   opt_amg/amg_summary.tsv     id, amg_flag, n_amg_genes, amg_classes   (per vOTU rep; report/votu_master contract)
#   opt_amg/amg_genes.tsv       gene-level AMG calls (gene, contig, product, function, amg_class, confidence)
#   opt_amg/amg_per_sample.tsv  per-sample x AMG-class summary (n_genes, n_votus, tpm)

_AMG = config.get("amg", {})
_AMG_THREADS = _AMG.get("threads", 16)
_AMG_PHOLD = bool(_AMG.get("phold", True))
_dbv = str(_AMG.get("db", "") or "")
# db may be a bare name under databases.dir (default "phold") or an absolute/relative path.
_PHOLD_DB = _dbv if _dbv.startswith("/") else f"{config['databases']['dir']}/{_dbv or 'phold'}"
_GPU = bool(_AMG.get("gpu", True))


if _AMG_PHOLD:

    rule amg_phold_db:
        output:
            f"{_PHOLD_DB}/.done",
        params:
            db=_PHOLD_DB,
            gpu="--foldseek_gpu" if _GPU else "",
        threads: _AMG_THREADS
        conda:
            env("phold")
        shell:
            # phold install: downloads ProstT5 model + foldseek structure DB (~16 GB) to {params.db}.
            "phold install -d {params.db} {params.gpu} -t {threads} && touch {output}"


    rule amg_phold:
        input:
            gbk=f"{OUT}/03_annotate/pharokka/pharokka.gbk",
            db=f"{_PHOLD_DB}/.done",
        output:
            preds=f"{OUT}/opt_amg/phold/phold_per_cds_predictions.tsv",
        params:
            outdir=f"{OUT}/opt_amg/phold",
            db=_PHOLD_DB,
            gpu="--foldseek_gpu" if _GPU else "--cpu",
        threads: _AMG_THREADS
        conda:
            env("phold")
        shell:
            # phold run = predict (ProstT5 3Di) + compare (foldseek). --hyps only re-annotates
            # hypothetical proteins from the Pharokka GenBank; this is the realistic use case for
            # metagenomic vOTUs (the long structural proteins are already PHROG-annotated and are the
            # bottleneck of a full run). --autotune picks the optimal ProstT5 batch size. -f overwrites.
            "phold run -i {input.gbk} -o {params.outdir} -d {params.db} -t {threads} {params.gpu} --hyps --autotune -f"


def _amg_inputs(wildcards):
    d = {
        "g2c": f"{OUT}/03_annotate/gene2genome.tsv",
        "clusters": f"{OUT}/02_cluster/votu_clusters.tsv",
        "c2s": f"{OUT}/01_qc/cleaned_contig2sample.tsv",
    }
    if _AMG_PHOLD:
        # phold predictions take precedence; Pharokka gene_families fill the rest
        d["ann"] = f"{OUT}/opt_amg/phold/phold_per_cds_predictions.tsv"
        d["ann2"] = f"{OUT}/03_annotate/gene_families.tsv"
    else:
        # Pharokka-only AMG classification (fast, no GPU/structure DB needed)
        d["ann"] = f"{OUT}/03_annotate/gene_families.tsv"
    # coverm becomes a real input only when abundance ran, so the rule re-fires when it
    # changes; the path is always passed as a param and tolerated as missing by the script.
    if flag("abundance"):
        d["coverm"] = f"{OUT}/06_abundance/coverm.tsv"
    return d


rule amg_classify:
    input:
        unpack(_amg_inputs),
    output:
        summary=f"{OUT}/opt_amg/amg_summary.tsv",
        genes=f"{OUT}/opt_amg/amg_genes.tsv",
        per_sample=f"{OUT}/opt_amg/amg_per_sample.tsv",
    params:
        script=scr("amg_classify.py"),
        phold_arg=lambda wc: f"--ann2 {OUT}/03_annotate/gene_families.tsv" if _AMG_PHOLD else "",
        coverm=f"{OUT}/06_abundance/coverm.tsv",
    conda:
        env("report")
    shell:
        "python {params.script} --ann {input.ann} {params.phold_arg} --g2c {input.g2c} "
        "--clusters {input.clusters} --c2s {input.c2s} --coverm {params.coverm} "
        "--out-summary {output.summary} --out-genes {output.genes} --out-sample {output.per_sample}"
