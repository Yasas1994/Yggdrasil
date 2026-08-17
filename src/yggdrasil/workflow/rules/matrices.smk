rule count_tables:
    input:
        c2s=f"{OUT}/01_qc/cleaned_contig2sample.tsv",
        clusters=f"{OUT}/02_cluster/votu_clusters.tsv",
        gf=f"{OUT}/03_annotate/gene_families.tsv",
        samples=config["samples"],
    output:
        votu=f"{OUT}/07_matrices/votu_sample_presence_absence.tsv",
        gene=f"{OUT}/07_matrices/gene_sample_counts.tsv",
    params:
        script=scr("build_count_tables.py"),
    shell:
        "python {params.script} --contig2sample {input.c2s} --votu-clusters {input.clusters} "
        "--gene-families {input.gf} --samples {input.samples} "
        "--out-votu {output.votu} --out-gene {output.gene}"
