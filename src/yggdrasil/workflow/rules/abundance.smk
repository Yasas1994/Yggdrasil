# Quantify QC'd reads against vOTU representatives with CoverM (CoverM wraps
# the mapper). One CoverM run per sample, then merged into a contig x sample
# matrix. Flag names follow wwood/CoverM 0.7.x; adjust if your version differs.

rule coverm_sample:
    input:
        r1=lambda wc: reads_of(wc.sample, 1),
        r2=lambda wc: reads_of(wc.sample, 2),
        ref=f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        temp(f"{OUT}/06_abundance/per_sample/{{sample}}.tsv"),
    params:
        mapper=config["abundance"].get("mapper", "minimap2-sr"),
    threads: config["abundance"].get("threads", 8)
    conda:
        env("coverm")
    shell:
        "coverm contig --coupled {input.r1} {input.r2} --reference {input.ref} "
        "--mapper {params.mapper} --methods tpm -t {threads} -o {output}"


rule coverm_merge:
    input:
        lambda wc: expand(
            f"{OUT}/06_abundance/per_sample/{{sample}}.tsv",
            sample=samples_with_reads(),
        ),
    output:
        f"{OUT}/06_abundance/coverm.tsv",
    run:
        # Runs inside the rule job (compute node env has pandas); the driver
        # process itself stays pandas-free so it can live on the login node.
        import pandas as pd
        merged = None
        for s in samples_with_reads():
            df = pd.read_csv(f"{OUT}/06_abundance/per_sample/{s}.tsv", sep="\t")
            df = df.rename(columns={df.columns[0]: "contig"})        # CoverM: 'Contig'
            val = [c for c in df.columns if c != "contig"][0]         # tpm column
            df = df[["contig", val]].rename(columns={val: s})
            merged = df if merged is None else merged.merge(df, on="contig", how="outer")
        if merged is None:
            merged = pd.DataFrame({"contig": []})
        merged.to_csv(output[0], sep="\t", index=False)
