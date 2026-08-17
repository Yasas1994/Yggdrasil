# ProkBERT PhaStyle lifestyle (temperate/virulent) for vOTU representatives.
# Runs inside the obalasz/phastyle image; a host-side normalizer coerces the raw
# prediction table to id,lifestyle,lifestyle_score. The model
# (neuralbioinfo/PhaStyle-mini) is fetched from Hugging Face on first run and cached.

rule phastyle:
    input:
        f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        f"{OUT}/05_lifestyle/phastyle_raw.tsv",
    threads: config["lifestyle"].get("threads", 8)
    container:
        config["lifestyle"].get("image", "docker://obalasz/phastyle:latest")
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
    threads: config["lifestyle"].get("threads", 8)
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
