# Ecological summaries from the vOTU x sample ABUNDANCE matrix (CoverM mean
# coverage). Counts are not meaningful here; diversity/ordination are computed
# on abundance. differential.tsv is produced by a separate rule that the master
# Snakefile only requests when ecology.group_col is set, so a null group_col
# never triggers it.

rule ecology:
    input:
        abundance=f"{OUT}/06_abundance/coverm.tsv",
        samples=config["samples"],
    output:
        alpha=f"{OUT}/08_ecology/alpha.tsv",
        beta=f"{OUT}/08_ecology/beta.tsv",
        ordination=f"{OUT}/08_ecology/ordination.png",
    params:
        script=scr("ecology.R"),
        gc=lambda wc: (
            f"--group-col {config['ecology']['group_col']}"
            if config.get("ecology", {}).get("group_col") else ""
        ),
    conda:
        env("ecology")
    shell:
        "Rscript {params.script} --abundance {input.abundance} --samples {input.samples} "
        "{params.gc} --out-alpha {output.alpha} --out-beta {output.beta} "
        "--out-ordination {output.ordination}"


rule ecology_differential:
    input:
        abundance=f"{OUT}/06_abundance/coverm.tsv",
        samples=config["samples"],
    output:
        f"{OUT}/08_ecology/differential.tsv",
    params:
        script=scr("ecology.R"),
        # Only demanded when group_col is set (gated in the Snakefile), so it is non-null here.
        group_col=lambda wc: config.get("ecology", {}).get("group_col"),
    conda:
        env("ecology")
    shell:
        "Rscript {params.script} --abundance {input.abundance} --samples {input.samples} "
        "--group-col {params.group_col} --out-differential {output}"
