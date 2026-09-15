# CheckV quality control. Two mutually exclusive layouts, chosen at parse time:
#   qc.per_sample: true  -> one checkv_sample job per sample (scales with sample
#                           count) + checkv_merge concatenating into the same four
#                           01_qc/checkv/ paths the single rule writes.
#   qc.per_sample: false -> the original single checkv job over all contigs.
# Only one branch defines its rules, so downstream rules see exactly one
# producer for each 01_qc/checkv/ file either way.

if bool(config["qc"].get("per_sample", True)):

    rule checkv_sample:
        input:
            # directory() output of preprocess; per-sample FASTA is resolved in
            # the shell (files inside a directory() output can't be rule inputs).
            per_sample=f"{OUT}/00_preprocess/per_sample",
        output:
            summary=f"{OUT}/01_qc/per_sample/{{sample}}/quality_summary.tsv",
            contamination=f"{OUT}/01_qc/per_sample/{{sample}}/contamination.tsv",
            viruses=f"{OUT}/01_qc/per_sample/{{sample}}/viruses.fna",
            proviruses=f"{OUT}/01_qc/per_sample/{{sample}}/proviruses.fna",
        params:
            outdir=f"{OUT}/01_qc/per_sample/{{sample}}",
            db=lambda wildcards: f"{config['databases']['dir']}/checkv",
        threads: config["qc"].get("sample_threads", 4)
        resources:
            mem_mb=8000,
            runtime=120,
        conda:
            env("checkv")
        shell:
            # Empty per-sample FASTA (sample had no usable contigs): emit header-only
            # TSVs (CheckV 1.0.x end_to_end column layout) and empty FASTAs without
            # calling checkv, so checkv_merge always has one file set per sample.
            # The DB dir resolves like the single rule: checkv download_database
            # unpacks to a versioned checkv-db-vX.Y subdir holding genome_db/.
            "FNA={input.per_sample}/{wildcards.sample}.fna; "
            "if [ ! -s $FNA ] || ! grep -q '^>' $FNA; then "
            "mkdir -p {params.outdir} && "
            "printf 'contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\n' > {output.summary} && "
            "printf 'contig_id\tcontig_length\ttotal_genes\tviral_genes\thost_genes\tprovirus\tproviral_length\thost_length\tregion_types\tregion_lengths\tregion_coords_bp\tregion_coords_genes\tregion_viral_genes\tregion_host_genes\n' > {output.contamination} && "
            "touch {output.viruses} {output.proviruses}; "
            "else "
            "CKDB=$(ls -d {params.db}/checkv-db-v* 2>/dev/null | head -n1); CKDB=${{CKDB:-{params.db}}}; "
            "checkv end_to_end $FNA {params.outdir} -t {threads} -d $CKDB "
            # summary must be non-empty (header + one row per contig); the FASTAs
            # only need to EXIST -- a sample can legitimately have zero proviruses
            # (empty proviruses.fna) and test -s would wrongly fail the job.
            "&& test -s {output.summary} && test -f {output.viruses} && test -f {output.proviruses}; "
            "fi"


    rule checkv_merge:
        input:
            summary=expand(f"{OUT}/01_qc/per_sample/{{sample}}/quality_summary.tsv", sample=sample_list()),
            contamination=expand(f"{OUT}/01_qc/per_sample/{{sample}}/contamination.tsv", sample=sample_list()),
            viruses=expand(f"{OUT}/01_qc/per_sample/{{sample}}/viruses.fna", sample=sample_list()),
            proviruses=expand(f"{OUT}/01_qc/per_sample/{{sample}}/proviruses.fna", sample=sample_list()),
        output:
            summary=f"{OUT}/01_qc/checkv/quality_summary.tsv",
            viruses=f"{OUT}/01_qc/checkv/viruses.fna",
            proviruses=f"{OUT}/01_qc/checkv/proviruses.fna",
            contamination=f"{OUT}/01_qc/checkv/contamination.tsv",
        shell:
            # awk 'NR==1 || FNR>1': keep the first file's header, drop the rest.
            "awk 'NR==1 || FNR>1' {input.summary} > {output.summary} "
            "&& awk 'NR==1 || FNR>1' {input.contamination} > {output.contamination} "
            "&& cat {input.viruses} > {output.viruses} "
            "&& cat {input.proviruses} > {output.proviruses}"


else:

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
        resources:
            mem_mb=32000,
            runtime=1440,
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
