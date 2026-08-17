rule checkv:
    input:
        f"{OUT}/00_preprocess/all_contigs.fna",
    output:
        summary=f"{OUT}/01_qc/checkv/quality_summary.tsv",
        viruses=f"{OUT}/01_qc/checkv/viruses.fna",
        proviruses=f"{OUT}/01_qc/checkv/proviruses.fna",
        contamination=f"{OUT}/01_qc/checkv/contamination.tsv",
    params:
        outdir=f"{OUT}/01_qc/checkv",
        db=lambda wildcards: f"{config['databases']['dir']}/checkv",
    threads: config["qc"]["threads"]
    conda:
        env("checkv")
    shell:
        # checkv download_database unpacks to a versioned subdir (checkv-db-vX.Y);
        # -d must point at the dir that contains genome_db/. Resolve it at runtime.
        "CKDB=$(ls -d {params.db}/checkv-db-v* 2>/dev/null | head -n1); CKDB=${{CKDB:-{params.db}}}; "
        "checkv end_to_end {input} {params.outdir} -t {threads} -d $CKDB "
        "&& test -s {output.summary} && test -s {output.viruses} && test -s {output.proviruses}"


rule clean_proviruses:
    input:
        viruses=f"{OUT}/01_qc/checkv/viruses.fna",
        proviruses=f"{OUT}/01_qc/checkv/proviruses.fna",
        c2s=f"{OUT}/00_preprocess/contig2sample.tsv",
    output:
        fna=f"{OUT}/01_qc/cleaned_contigs.fna",
        c2s=f"{OUT}/01_qc/cleaned_contig2sample.tsv",
        lengths=f"{OUT}/01_qc/cleaned_lengths.tsv",
    params:
        script=scr("clean_proviruses.py"),
    shell:
        # Remove host regions from CheckV proviruses. IDs are normalized back to the
        # original contig IDs so downstream sample mapping stays intact.
        "python {params.script} --viruses {input.viruses} --proviruses {input.proviruses} "
        "--c2s {input.c2s} --out-fna {output.fna} --out-c2s {output.c2s} --out-lengths {output.lengths}"
