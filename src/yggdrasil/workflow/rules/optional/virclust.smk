# VirClust (Moraru 2023): protein-based hierarchical clustering of the dereplicated
# vOTU representatives -> hierarchical tree + viral genome clusters (VGCs) + core proteins.
# Gated by config virclust.enabled. Standalone R scripts cloned by setup_databases.sh
# into databases.dir/virclust/repo. Branch A (PC-level) only: skips the HMM/PSC/PSSC
# and annotation branches (those need extra DBs). VirClust_MASTER.R locates its siblings
# via this.path::here(), so it can be invoked from any cwd.


rule virclust:
    input:
        reps=f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        tree=f"{OUT}/opt_virclust/tree.newick",
        dist=f"{OUT}/opt_virclust/distance_matrix.tsv",
        clusters=f"{OUT}/opt_virclust/genome_clusters.tsv",
        core=f"{OUT}/opt_virclust/core_proteins.faa",
    params:
        master=lambda wildcards: f"{config['databases']['dir']}/virclust/repo/VirClust/vir_clust_standalone/VirClust_MASTER.R",
        work=lambda wildcards: f"{OUT}/opt_virclust/run",
        gene_code=config["virclust"]["gene_code"],
    threads: config["virclust"]["threads"]
    conda:
        env("virclust")
    shell:
        "test -f {params.master} || "
        "{{ echo 'VirClust scripts missing; run: yggdrasil setup-databases' >&2; exit 1; }}; "
        "WORK=$(realpath -m {params.work}); REPS=$(realpath {input.reps}); "
        "TREE=$(realpath -m {output.tree}); DIST=$(realpath -m {output.dist}); CLUST=$(realpath -m {output.clusters}); "
        "CORE=$(realpath -m {output.core}); "
        "Rscript {params.master} projdir=$WORK infile=$REPS "
        "sing=conda condaenvpath= "
        "step1A=T step2A=T step3A=T step4A=T step5A=T step3A_Plot=F step4A_Plot=F "
        "cpu={threads} gene_code={params.gene_code} show_tree=no show_heat=no "
        "&& cp $WORK/04a-06a_genome_clustering_PC/04/hc_tree.newick $TREE "
        "&& cp $WORK/04a-06a_genome_clustering_PC/04/MyDistPCs_MA.tsv $DIST "
        "&& cp $WORK/04a-06a_genome_clustering_PC/05/virDF.tsv $CLUST "
        "&& cp $WORK/07/core_a/core_prots_for_annots_all.faa $CORE"
