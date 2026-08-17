# Gene calling + functional annotation of vOTU representatives (cost control).
# Branch on annotate.pharokka:
#   ON  -> pharokka (prodigal-gv genes + PHROG via pyhmmer/mmseqs) in one env.
#   OFF -> prodigal-gv calls genes; mmseqs2 searches the PHROG profile DB.
# The pharokka post-process is a script (pharokka_gene_families.py) because
# pharokka renames genes and the merged CDS table holds the authoritative
# gene -> contig + PHROG mapping; inline shell quoting is too fragile here.

_A = config["annotate"]
_THREADS = _A.get("threads", 16)
_DB = config["databases"]["dir"]


if bool(_A.get("pharokka")):

    rule annotate_pharokka:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            faa=f"{OUT}/03_annotate/proteins.faa",
            gff=f"{OUT}/03_annotate/genes.gff",
            gf=f"{OUT}/03_annotate/gene_families.tsv",
            g2g=f"{OUT}/03_annotate/gene2genome.tsv",
        params:
            outdir=f"{OUT}/03_annotate/pharokka",
            db=lambda wc: f"{_DB}/phrog",
            script=scr("pharokka_gene_families.py"),
        threads: _THREADS
        conda:
            env("pharokka")
        shell:
            # pharokka>=1.7: -g prodigal-gv for viral gene calling; proteins land in prodigal-gv.faa.
            # -f overwrites the outdir so rule re-runs (code/input change) don't error.
            "pharokka.py -i {input} -o {params.outdir} -d {params.db} -t {threads} -g prodigal-gv -f"
            " && cp $(ls {params.outdir}/*.faa | grep -v terL | head -n1) {output.faa}"
            " && cp {params.outdir}/pharokka.gff {output.gff}"
            " && python {params.script} --indir {params.outdir} --out {output.gf} --out-g2g {output.g2g}"


else:

    rule call_genes:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            faa=f"{OUT}/03_annotate/proteins.faa",
            gff=f"{OUT}/03_annotate/genes.gff",
        threads: _THREADS
        conda:
            env("prodigal-gv")
        shell:
            # prodigal-gv: viral gene calling; -p meta for mixed contigs; -a proteins, -o GFF features.
            "prodigal-gv -i {input} -a {output.faa} -o {output.gff} -f gff -p meta -q"


    if bool(_A.get("phrog")):

        rule search_phrog:
            input:
                f"{OUT}/03_annotate/proteins.faa",
            output:
                f"{OUT}/03_annotate/phrog.m8",
            params:
                db=lambda wc: f"{_DB}/phrog",
                tmp=f"{OUT}/03_annotate/mmseqs_tmp",
            threads: _THREADS
            conda:
                env("mmseqs2")
            shell:
                # mmseqs2>=15: phrog DB is an MMseqs2 profile DB prefix at {db}/phrog.
                "mmseqs easy-search {input} {params.db} {output} {params.tmp}"
                " --threads {threads} --format-output \"query,target,pident,evalue,bits\""


    else:

        rule search_phrog:
            output:
                f"{OUT}/03_annotate/phrog.m8",
            shell:
                "touch {output}"


    rule build_families:
        input:
            faa=f"{OUT}/03_annotate/proteins.faa",
            m8=f"{OUT}/03_annotate/phrog.m8",
        output:
            f"{OUT}/03_annotate/gene_families.tsv",
        shell:
            # gene_id == prodigal protein header token (seqid_idx); contig = drop trailing _idx.
            "python -c 'import csv,sys\n"
            "faa,m8,out=sys.argv[1:4]\n"
            "hits=dict()\n"
            "for r in csv.reader(open(m8),delimiter=\"\t\"):\n"
            "    if len(r)<5: continue\n"
            "    b=float(r[4])\n"
            "    if r[0] not in hits or b>hits[r[0]][0]: hits[r[0]]=(b,r[1])\n"
            "w=csv.writer(open(out,\"w\",newline=\"\"),delimiter=\"\t\")\n"
            "w.writerow([\"contig\",\"gene_id\",\"family\"])\n"
            "for line in open(faa):\n"
            "    if line.startswith(\">\"):\n"
            "        gid=line[1:].split()[0]\n"
            "        w.writerow([gid.rsplit(\"_\",1)[0],gid,hits.get(gid,(0,\"\"))[1]])\n"
            "' {input.faa} {input.m8} {output}"


    rule build_g2g_prodigal:
        input:
            faa=f"{OUT}/03_annotate/proteins.faa",
            m8=f"{OUT}/03_annotate/phrog.m8",
        output:
            f"{OUT}/03_annotate/gene2genome.tsv",
        params:
            script=scr("prodigal_g2g.py"),
        shell:
            "python {params.script} --faa {input.faa} --m8 {input.m8} --out {output}"
