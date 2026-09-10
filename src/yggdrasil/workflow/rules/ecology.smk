# Ecological summaries from a vOTU x sample matrix: CoverM TPM when the
# abundance module is enabled AND at least one sample has reads; otherwise the
# vOTU x sample presence/absence matrix (diversity on 0/1 is standard; no
# CoverM jobs run). Counts are not meaningful here; diversity/ordination are
# computed on abundance (or presence/absence). differential.tsv is produced by
# a separate rule that the master Snakefile only requests when
# ecology.group_col is set, so a null group_col never triggers it.

_ECO_MATRIX = (
    f"{OUT}/06_abundance/coverm.tsv"
    if flag("abundance") and samples_with_reads()
    else f"{OUT}/07_matrices/votu_sample_presence_absence.tsv"
)

rule ecology:
    input:
        abundance=_ECO_MATRIX,
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
        abundance=_ECO_MATRIX,
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
