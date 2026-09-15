def _votu_master_inputs(wildcards):
    d = {
        "clusters": f"{OUT}/02_cluster/votu_clusters.tsv",
        "c2s": f"{OUT}/01_qc/cleaned_contig2sample.tsv",
        "lengths": f"{OUT}/01_qc/cleaned_lengths.tsv",
        "quality": f"{OUT}/01_qc/checkv/quality_summary.tsv",
        "taxonomy": f"{OUT}/04_taxonomy/taxonomy.tsv",
        "lifestyle": f"{OUT}/05_lifestyle/lifestyle.tsv",
    }
    # Optional modules become real inputs when enabled, so the master rebuilds
    # exactly when their outputs change (host/amg paths are still passed as params
    # and tolerated as missing by votu_master.py when the module is off).
    if flag("host"):
        d["host"] = f"{OUT}/opt_host/iphop.tsv"
    if flag("amg"):
        d["amg"] = f"{OUT}/opt_amg/amg_summary.tsv"
    return d


rule votu_master:
    input:
        unpack(_votu_master_inputs),
    output:
        f"{OUT}/09_report/votu_master.tsv",
    params:
        script=scr("votu_master.py"),
        amg=f"{OUT}/opt_amg/amg_summary.tsv",
        host=f"{OUT}/opt_host/iphop.tsv",
    shell:
        "python {params.script} --votu-clusters {input.clusters} --contig2sample {input.c2s} "
        "--lengths {input.lengths} --quality {input.quality} --taxonomy {input.taxonomy} "
        "--lifestyle {input.lifestyle} --amg {params.amg} --host {params.host} --out {output}"


def _report_inputs(wildcards):
    d = {
        "master": f"{OUT}/09_report/votu_master.tsv",
        "votu": f"{OUT}/07_matrices/votu_sample_presence_absence.tsv",
        "gene": f"{OUT}/07_matrices/gene_sample_counts.tsv",
        "alpha": f"{OUT}/08_ecology/alpha.tsv",
        # vcontact3 rule outputs a directory(); depend on the directory and let
        # build_report.py read the inner files (graphml / final_assignments) at runtime.
        "vc_dir": f"{OUT}/04_taxonomy/vcontact3",
    }
    if flag("viridic"):
        d["viridic_sim"] = f"{OUT}/opt_viridic/similarity_matrix.tsv"
        d["viridic_clusters"] = f"{OUT}/opt_viridic/clusters.csv"
    if flag("virclust"):
        d["virclust_tree"] = f"{OUT}/opt_virclust/tree.newick"
        d["virclust_clusters"] = f"{OUT}/opt_virclust/genome_clusters.tsv"
    if flag("amg"):
        d["amg_summary"] = f"{OUT}/opt_amg/amg_summary.tsv"
        d["amg_genes"] = f"{OUT}/opt_amg/amg_genes.tsv"
        d["amg_sample"] = f"{OUT}/opt_amg/amg_per_sample.tsv"
    return d


rule report:
    input:
        unpack(_report_inputs),
    output:
        html=f"{OUT}/09_report/phage_report.html",
        assets=directory(f"{OUT}/09_report/assets"),
    params:
        script=scr("build_report.py"),
        vc_graphml=f"{OUT}/04_taxonomy/vcontact3/exports/networks/part1.graphml",
        vc_refs=f"{OUT}/04_taxonomy/vcontact3/exports/final_assignments.csv",
        viridic_sim=f"{OUT}/opt_viridic/similarity_matrix.tsv",
        viridic_clusters=f"{OUT}/opt_viridic/clusters.csv",
        virclust_tree=f"{OUT}/opt_virclust/tree.newick",
        virclust_clusters=f"{OUT}/opt_virclust/genome_clusters.tsv",
        amg_votu=f"{OUT}/opt_amg/amg_summary.tsv",
        amg_genes=f"{OUT}/opt_amg/amg_genes.tsv",
        amg_sample=f"{OUT}/opt_amg/amg_per_sample.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    conda:
        env("report")
    shell:
        "python {params.script} --master {input.master} --votu {input.votu} "
        "--gene {input.gene} --alpha {input.alpha} --vc-graphml {params.vc_graphml} "
        "--vc-refs {params.vc_refs} --viridic-sim {params.viridic_sim} "
        "--viridic-clusters {params.viridic_clusters} --virclust-tree {params.virclust_tree} "
        "--virclust-clusters {params.virclust_clusters} --amg-votu {params.amg_votu} "
        "--amg-genes {params.amg_genes} --amg-sample {params.amg_sample} --out {output.html}"
