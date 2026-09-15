# ProkBERT PhaStyle lifestyle (temperate/virulent) for vOTU representatives.
# Runs inside the obalasz/phastyle image; a host-side normalizer coerces the raw
# prediction table to id,lifestyle,lifestyle_score. The model
# (neuralbioinfo/PhaStyle-mini) is fetched from Hugging Face on first run and cached.
#
# lifestyle.seqs_per_chunk > 0 (default) splits the reps FASTA into chunks
# (checkpoint split_reps_lifestyle) and runs PhaStyle + BACPHLIP per chunk;
# lifestyle_merge concatenates (header-dedup) into the same paths the single
# rules write. 0 = the original single jobs over all representatives.

_L = config["lifestyle"]
_L_THREADS = _L.get("threads", 8)

if int(_L.get("seqs_per_chunk", 2000)) > 0:

    checkpoint split_reps_lifestyle:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            directory(f"{OUT}/05_lifestyle/split"),
        params:
            script=scr("split_fasta.py"),
            k=int(_L.get("seqs_per_chunk", 2000)),
        shell:
            "python {params.script} --fasta {input} --seqs-per-chunk {params.k} "
            "--outdir {output} --prefix reps"


    def _lifestyle_chunk_ids(wildcards):
        ck = checkpoints.split_reps_lifestyle.get(**wildcards).output[0]
        return Path(ck, "chunks.txt").read_text().split()


    rule phastyle_chunk:
        input:
            f"{OUT}/05_lifestyle/split/reps_{{chunk}}.fna",
        output:
            f"{OUT}/05_lifestyle/chunks/{{chunk}}.phastyle.tsv",
        threads: _L_THREADS
        resources:
            mem_mb=16000,
            runtime=720,
        container:
            _L.get("image", "docker://obalasz/phastyle:latest")
        shell:
            # Same PhaStyle CLI as the single rule, one chunk at a time.
            "python /home/prokbert/PhaStyle/bin/PhaStyle.py "
            "--fastain {input} --out {output} "
            "--ftmodel neuralbioinfo/PhaStyle-mini "
            "--batch-size 128 -c {threads}"


    rule bacphlip_chunk:
        input:
            f"{OUT}/05_lifestyle/split/reps_{{chunk}}.fna",
        output:
            f"{OUT}/05_lifestyle/chunks/{{chunk}}.bacphlip.tsv",
        threads: _L_THREADS
        resources:
            mem_mb=8000,
            runtime=360,
        conda:
            env("bacphlip")
        params:
            script=scr("bacphlip_predict.py"),
        shell:
            "python {params.script} --fasta {input} --out {output} --threads {threads}"


    rule lifestyle_merge:
        input:
            phastyle=lambda wc: expand(
                f"{OUT}/05_lifestyle/chunks/{{chunk}}.phastyle.tsv",
                chunk=_lifestyle_chunk_ids(wc),
            ),
            bacphlip=lambda wc: expand(
                f"{OUT}/05_lifestyle/chunks/{{chunk}}.bacphlip.tsv",
                chunk=_lifestyle_chunk_ids(wc),
            ),
        output:
            phastyle=f"{OUT}/05_lifestyle/phastyle_raw.tsv",
            bacphlip=f"{OUT}/05_lifestyle/bacphlip.tsv",
        shell:
            # awk 'NR==1 || FNR>1': keep the first file's header, drop the rest.
            # Zero chunks (empty reps) -> empty files; lifestyle_normalize tolerates both.
            "if [ -n '{input.phastyle}' ]; then awk 'NR==1 || FNR>1' {input.phastyle} > {output.phastyle}; else touch {output.phastyle}; fi "
            "&& if [ -n '{input.bacphlip}' ]; then awk 'NR==1 || FNR>1' {input.bacphlip} > {output.bacphlip}; else touch {output.bacphlip}; fi"


else:

    rule phastyle:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            f"{OUT}/05_lifestyle/phastyle_raw.tsv",
        threads: _L_THREADS
        resources:
            mem_mb=16000,
            runtime=720,
        container:
            _L.get("image", "docker://obalasz/phastyle:latest")
        shell:
            # PhaStyle.py CLI (container build): --fastain/--out/--ftmodel/--batch-size/-c.
            # (The repo README shows --per_device_eval_batch_size; the image uses --batch-size.)
            # CPU is fine; add --nv via singularity args if a GPU is available.
            "python /home/prokbert/PhaStyle/bin/PhaStyle.py "
            "--fastain {input} --out {output} "
            "--ftmodel neuralbioinfo/PhaStyle-mini "
            "--batch-size 128 -c {threads}"


    rule bacphlip:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            f"{OUT}/05_lifestyle/bacphlip.tsv",
        threads: _L_THREADS
        resources:
            mem_mb=8000,
            runtime=360,
        conda:
            env("bacphlip")
        params:
            script=scr("bacphlip_predict.py"),
        shell:
            "python {params.script} --fasta {input} --out {output} --threads {threads}"


rule lifestyle_normalize:
    input:
        phastyle=f"{OUT}/05_lifestyle/phastyle_raw.tsv",
        bacphlip=f"{OUT}/05_lifestyle/bacphlip.tsv",
    output:
        f"{OUT}/05_lifestyle/lifestyle.tsv",
    params:
        script=scr("lifestyle_normalize.py"),
    shell:
        "python {params.script} --raw {input.phastyle} --bacphlip {input.bacphlip} --out {output}"
