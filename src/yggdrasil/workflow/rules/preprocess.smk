rule preprocess:
    input:
        samples=config["samples"],
    output:
        fna=f"{OUT}/00_preprocess/all_contigs.fna",
        cmap=f"{OUT}/00_preprocess/contig2sample.tsv",
        lengths=f"{OUT}/00_preprocess/lengths.tsv",
        per_sample=directory(f"{OUT}/00_preprocess/per_sample"),
    params:
        script=scr("preprocess.py"),
        min_length=config["preprocess"]["min_length"],
    shell:
        "python {params.script} --samples {input.samples} --min-length {params.min_length} "
        + ("--dedup " if config["preprocess"].get("dedup") else "")
        + "--out-fna {output.fna} --out-map {output.cmap} --out-lengths {output.lengths} "
        "--out-per-sample {output.per_sample}"
