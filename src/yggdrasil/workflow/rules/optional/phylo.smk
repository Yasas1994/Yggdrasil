# Optional terminase large-subunit (TerL) phylogeny. Gated by config phylo.enabled in
# the master Snakefile. mafft>=7.5 + IQ-TREE>=2.2 (binary `iqtree2`).
# Inputs come from rules/annotate.smk (gene_families.tsv / proteins.faa).

rule phylo_extract_terl:
    input:
        families=f"{OUT}/03_annotate/gene_families.tsv",
        proteins=f"{OUT}/03_annotate/proteins.faa",
    output:
        f"{OUT}/opt_phylo/terl.faa",
    params:
        script=scr("extract_terl.py"),
    shell:
        "python {params.script} --families {input.families} --proteins {input.proteins} --out {output}"


rule phylo_align:
    input:
        f"{OUT}/opt_phylo/terl.faa",
    output:
        aln=f"{OUT}/opt_phylo/terl.aln.faa",
    params:
        log=f"{OUT}/opt_phylo/mafft.log",
    threads: config["phylo"]["threads"]
    conda:
        env("mafft")
    shell:
        "n=$(grep -c '^>' {input} 2>/dev/null || echo 0); "
        "if [ \"${{n:-0}}\" -lt 4 ]; then : > {output.aln}; "
        "else mafft --auto --thread {threads} {input} > {output.aln} 2> {params.log}; fi"


rule phylo_tree:
    input:
        aln=f"{OUT}/opt_phylo/terl.aln.faa",
        terl=f"{OUT}/opt_phylo/terl.faa",
    output:
        treefile=f"{OUT}/opt_phylo/terl.treefile",
        note=f"{OUT}/opt_phylo/terl.note.txt",
    params:
        prefix=f"{OUT}/opt_phylo/terl",
    threads: config["phylo"]["threads"]
    conda:
        env("iqtree")
    shell:
        "if [ ! -s {input.aln} ]; then "
        ": > {output.treefile}; "
        "n=$(grep -c '^>' {input.terl} 2>/dev/null || echo 0); "
        "echo \"TerL phylogeny skipped: only ${{n:-0}} TerL sequence(s); need >=4.\" > {output.note}; "
        "else "
        "iqtree2 -s {input.aln} -m LG+G -fast -nt {threads} -pre {params.prefix} -redo -quiet; "
        "echo \"TerL phylogeny complete (iqtree2 -m LG+G -fast).\" > {output.note}; fi"
