# VIRIDIC (Moraru et al. 2020): ICTV-style intergenomic similarities + species/genus
# clustering of the dereplicated vOTU representatives. Gated by config viridic.enabled.
# Standalone R scripts are cloned by setup_databases.sh into databases.dir/viridic/repo.
# The master script uses getwd(), so we cd into the scripts dir before invoking.


rule viridic:
    input:
        reps=f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        clusters=f"{OUT}/opt_viridic/clusters.csv",
        sim=f"{OUT}/opt_viridic/similarity_matrix.tsv",
        heatmap=f"{OUT}/opt_viridic/heatmap.pdf",
    params:
        scripts=lambda wildcards: f"{config['databases']['dir']}/viridic/repo/VIRIDIC/stand_alone/viridic_scripts",
        work=lambda wildcards: f"{OUT}/opt_viridic/run",
        thsp=config["viridic"]["species_threshold"],
        thgen=config["viridic"]["genus_threshold"],
        heatmap_script=scr("viridic_heatmap.R"),
    threads: config["viridic"]["threads"]
    conda:
        env("viridic")
    shell:
        "test -f {params.scripts}/00_viridic_master.R || "
        "{{ echo 'VIRIDIC scripts missing; run: yggdrasil setup-databases' >&2; exit 1; }}; "
        "WORK=$(realpath -m {params.work}); REPS=$(realpath {input.reps}); "
        "CLUST=$(realpath -m {output.clusters}); SIM=$(realpath -m {output.sim}); HM=$(realpath -m {output.heatmap}); "
        # steps=sim_clust: run similarity + clustering only. The native heatmap
        # step (03_viridic_heatmap.R) renders a per-cell ComplexHeatmap vector
        # PDF that exhausts memory on large vOTU sets; we render our own
        # rasterized PDF below instead.
        "cd {params.scripts} && "
        "Rscript 00_viridic_master.R projdir=$WORK in=$REPS "
        "ncor={threads} thsp={params.thsp} thgen={params.thgen} steps=sim_clust "
        "&& cp $WORK/04_VIRIDIC_out/clusters.csv $CLUST "
        "&& cp $WORK/04_VIRIDIC_out/sim_MA_genCol.csv $SIM "
        "&& Rscript {params.heatmap_script} $SIM $CLUST $HM"
