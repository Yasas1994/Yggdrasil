# vContact3 gene-sharing taxonomy for vOTU representatives.
# Proteins + a gene-to-genome map feed `vcontact3 run`; a pure-python normalizer
# maps whatever assignment vContact3 emits back onto the representative contig ids.

rule vcontact3:
    input:
        faa=f"{OUT}/03_annotate/proteins.faa",
        g2g=f"{OUT}/03_annotate/gene2genome.tsv",
        reps=f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        directory(f"{OUT}/04_taxonomy/vcontact3"),
    params:
        outdir=f"{OUT}/04_taxonomy/vcontact3",
        taxdir=f"{OUT}/04_taxonomy",
        db=lambda wc: f"{config['databases']['dir']}/vcontact3",
        db_version=lambda wc: config["taxonomy"].get("db_version", 230),
        len_script=scr("fasta_lengths.py"),
        t2p=scr("tsv2parquet.py"),
    threads: config["taxonomy"]["threads"]
    resources:
        mem_mb=32000,
        runtime=1440,
    conda:
        env("vcontact3")
    shell:
        # vcontact3 3.1.x: --gene2genome (protein_id,genome_id,keywords) and --len-nucleotide
        # (genome_id,length in Kb) are required with --proteins, and BOTH must be parquet
        # (vcontact3 reads them with pd.read_parquet). Intermediates live in taxdir (NOT the
        # output dir) because vcontact3 copies the inputs into --output; same-path = error.
        # gene2genome is emitted by annotate so gene -> original-contig survives pharokka's
        # per-run gene renaming. pandas pinned <3 (vcontact3 3.1.x str->float assignment).
        "mkdir -p {params.outdir}"
        " && python {params.len_script} --fasta {input.reps} --out {params.taxdir}/rep_lengths.tsv --kb"
        " && python {params.t2p} {input.g2g} {params.taxdir}/gene2genome.parquet"
        " && python {params.t2p} {params.taxdir}/rep_lengths.tsv {params.taxdir}/rep_lengths.parquet"
        " && vcontact3 run --proteins {input.faa} --gene2genome {params.taxdir}/gene2genome.parquet"
        " --len-nucleotide {params.taxdir}/rep_lengths.parquet"
        " --db-path {params.db} --db-version {params.db_version} --db-domain prokaryotes"
        " --output {params.outdir} --exports profiles graphml --threads {threads} -f"


rule normalize_taxonomy:
    input:
        vc=f"{OUT}/04_taxonomy/vcontact3",
        reps=f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        f"{OUT}/04_taxonomy/taxonomy.tsv",
    params:
        script=scr("normalize_taxonomy.py"),
    shell:
        "python {params.script} --vcontact3-dir {input.vc} --representatives {input.reps} --out {output}"
