# vclust: prefilter -> align -> cluster, then quality-aware representative selection.
# CLI verified against refresh-bio/vclust (align writes ani.tsv + ani.ids.tsv sidecar;
# cluster reads --ids that sidecar). Thresholds are 0-1 fractions in vclust.

rule vclust:
    input:
        f"{OUT}/01_qc/cleaned_contigs.fna",
    output:
        clusters=f"{OUT}/02_cluster/vclust/clusters.tsv",
    params:
        outdir=f"{OUT}/02_cluster/vclust",
        ani=lambda wildcards: config["cluster"]["ani"] / 100.0,
        af=config["cluster"]["align_frac"],
    threads: config["cluster"]["threads"]
    conda:
        env("vclust")
    shell:
        "vclust prefilter -i {input} -o {params.outdir}/fltr.txt -t {threads} "
        "&& vclust align -i {input} -o {params.outdir}/ani.tsv --filter {params.outdir}/fltr.txt -t {threads} "
        "&& vclust cluster -i {params.outdir}/ani.tsv -o {output.clusters} --ids {params.outdir}/ani.ids.tsv "
        "--algorithm single --metric ani --ani {params.ani} --qcov {params.af}"


rule select_representatives:
    input:
        clusters=f"{OUT}/02_cluster/vclust/clusters.tsv",
        quality=f"{OUT}/01_qc/checkv/quality_summary.tsv",
        lengths=f"{OUT}/01_qc/cleaned_lengths.tsv",
        fasta=f"{OUT}/01_qc/cleaned_contigs.fna",
    output:
        clusters=f"{OUT}/02_cluster/votu_clusters.tsv",
        fasta=f"{OUT}/02_cluster/votu_representatives.fna",
    params:
        script=scr("select_representatives.py"),
    shell:
        "python {params.script} --clusters {input.clusters} --quality {input.quality} "
        "--lengths {input.lengths} --fasta {input.fasta} "
        "--out-clusters {output.clusters} --out-fasta {output.fasta}"
