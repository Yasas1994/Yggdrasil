# Gene calling + functional annotation of vOTU representatives (cost control).
# Branch on annotate.pharokka:
#   ON  -> pharokka (prodigal-gv genes + PHROG via pyhmmer/mmseqs) in one env.
#          annotate.seqs_per_chunk > 0 (default) splits the reps FASTA into
#          chunks (checkpoint split_reps_annotate) and runs pharokka per chunk,
#          then annotate_merge rebuilds a merged pharokka dir; 0 = one big run.
#   OFF -> prodigal-gv calls genes; mmseqs2 searches the PHROG profile DB.
# The pharokka post-process is a script (pharokka_gene_families.py) because
# pharokka renames genes and the merged CDS table holds the authoritative
# gene -> contig + PHROG mapping; inline shell quoting is too fragile here.

_A = config["annotate"]
_THREADS = _A.get("threads", 16)
_DB = config["databases"]["dir"]


if bool(_A.get("pharokka")):

    if int(_A.get("seqs_per_chunk", 2000)) > 0:

        checkpoint split_reps_annotate:
            input:
                f"{OUT}/02_cluster/votu_representatives.fna",
            output:
                directory(f"{OUT}/03_annotate/split"),
            params:
                script=scr("split_fasta.py"),
                k=int(_A.get("seqs_per_chunk", 2000)),
            shell:
                # Chunk count is derived from the input size (ceil(n/K)); chunks.txt
                # lists the chunk ids and is read by _pharokka_chunks() below once
                # the checkpoint has run (standard data-dependent fan-out idiom).
                "python {params.script} --fasta {input} --seqs-per-chunk {params.k} "
                "--outdir {output} --prefix reps"


        def _pharokka_chunks(wildcards):
            ck = checkpoints.split_reps_annotate.get(**wildcards).output[0]
            ids = Path(ck, "chunks.txt").read_text().split()
            return expand(f"{OUT}/03_annotate/chunks/{{chunk}}", chunk=ids)


        rule annotate_pharokka_chunk:
            input:
                f"{OUT}/03_annotate/split/reps_{{chunk}}.fna",
            output:
                directory(f"{OUT}/03_annotate/chunks/{{chunk}}"),
            params:
                db=lambda wc: f"{_DB}/phrog",
            threads: _A.get("chunk_threads", 8)
            resources:
                mem_mb=32000,
                runtime=1440,
            conda:
                env("pharokka")
            shell:
                # Full pharokka run per chunk; -f overwrites so re-runs don't error.
                "pharokka.py -i {input} -o {output} -d {params.db} -t {threads} -g prodigal-gv -f"


        rule annotate_merge:
            input:
                _pharokka_chunks,
            output:
                faa=f"{OUT}/03_annotate/proteins.faa",
                gff=f"{OUT}/03_annotate/genes.gff",
                gf=f"{OUT}/03_annotate/gene_families.tsv",
                g2g=f"{OUT}/03_annotate/gene2genome.tsv",
            params:
                outdir=f"{OUT}/03_annotate/pharokka",
                chunks_dir=f"{OUT}/03_annotate/chunks",
                script=scr("pharokka_gene_families.py"),
            shell:
                # Fabricate a merged pharokka outdir from the per-chunk runs: cds TSV
                # with a single header (awk 'NR==1 || FNR>1'), plus concatenated
                # faa (excluding terL, as in the single rule) / gff / gbk. Gene ids
                # are per-run-hashed <HASH>_CDS_NNNN, so chunk hashes don't collide.
                # The same post-process script as the single rule runs on the merged dir.
                # Zero chunks (empty reps) -> empty merged files, header-only tables.
                "mkdir -p {params.outdir}; "
                "TSVS=''; FAAS=''; GFFS=''; GBKS=''; "
                "for d in {params.chunks_dir}/*/; do "
                "[ -d \"$d\" ] || continue; "
                "t=\"$d/pharokka_cds_final_merged_output.tsv\"; [ -s \"$t\" ] && TSVS=\"$TSVS $t\"; "
                "f=$(ls \"$d\"/*.faa 2>/dev/null | grep -v terL | head -n1); [ -n \"$f\" ] && FAAS=\"$FAAS $f\"; "
                "[ -s \"$d/pharokka.gff\" ] && GFFS=\"$GFFS $d/pharokka.gff\"; "
                "[ -s \"$d/pharokka.gbk\" ] && GBKS=\"$GBKS $d/pharokka.gbk\"; "
                "done; "
                "if [ -n \"$TSVS\" ]; then awk 'NR==1 || FNR>1' $TSVS > {params.outdir}/pharokka_cds_final_merged_output.tsv; fi; "
                "if [ -n \"$FAAS\" ]; then cat $FAAS > {params.outdir}/pharokka.faa; else touch {params.outdir}/pharokka.faa; fi; "
                "if [ -n \"$GFFS\" ]; then cat $GFFS > {params.outdir}/pharokka.gff; else touch {params.outdir}/pharokka.gff; fi; "
                "if [ -n \"$GBKS\" ]; then cat $GBKS > {params.outdir}/pharokka.gbk; else touch {params.outdir}/pharokka.gbk; fi; "
                "cp {params.outdir}/pharokka.faa {output.faa} && cp {params.outdir}/pharokka.gff {output.gff} "
                "&& python {params.script} --indir {params.outdir} --out {output.gf} --out-g2g {output.g2g}"


    else:

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
            resources:
                mem_mb=32000,
                runtime=1440,
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
