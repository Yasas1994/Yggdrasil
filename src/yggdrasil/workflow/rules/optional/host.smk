# Optional iPHoP (simroux/iphop) host prediction. Gated by config host.enabled.
# iPHoP 1.4 CLI: `iphop predict --fa_file <fna> --db_dir <db> --out_dir <dir> --num_threads <n>`.
# iPHoP unpacks into a VERSIONED subdir (<db_dir>/iPHoP_db_<ver>_rw/ or Test_db_rw_v1.4/) that holds
# the db/ + db_infos/ folders. The rule auto-resolves that subdir under <databases.dir>/<host.db>, so
# the config path stays stable across DB versions.
# The genus-level CSV (Host_prediction_to_genus_m<score>.csv) is normalized to id,host,host_score
# by scripts/iphop_normalize.py so report.smk's votu_master can join on representative id.
# iPHoP's final classifier uses TensorFlow and grabs the GPU by default; we force CPU (small-VRAM
# cards OOM otherwise) via CUDA_VISIBLE_DEVICES="".
#
# host.seqs_per_chunk > 0 (default) splits the reps FASTA into chunks (checkpoint
# split_reps_host) and runs iPHoP per chunk; host_merge concatenates the per-chunk
# genus CSVs into opt_host/run/ and normalizes exactly like the single rule.
# 0 = the original single job over all representatives.

_H = config["host"]

if int(_H.get("seqs_per_chunk", 2000)) > 0:

    checkpoint split_reps_host:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            directory(f"{OUT}/opt_host/split"),
        params:
            script=scr("split_fasta.py"),
            k=int(_H.get("seqs_per_chunk", 2000)),
        shell:
            "python {params.script} --fasta {input} --seqs-per-chunk {params.k} "
            "--outdir {output} --prefix reps"


    def _host_chunk_ids(wildcards):
        ck = checkpoints.split_reps_host.get(**wildcards).output[0]
        return Path(ck, "chunks.txt").read_text().split()


    rule iphop_chunk:
        input:
            f"{OUT}/opt_host/split/reps_{{chunk}}.fna",
        output:
            directory(f"{OUT}/opt_host/chunks/{{chunk}}"),
        params:
            db=lambda wildcards: f"{config['databases']['dir']}/{_H.get('db', 'iphop')}",
        threads: _H["threads"]
        resources:
            mem_mb=64000,
            runtime=2880,
        conda:
            env("iphop")
        shell:
            "mkdir -p {output}; "
            "DBROOT=$(d=$(find {params.db} -maxdepth 2 -type d -name db_infos 2>/dev/null | head -1); printf '%s' \"${{d%/db_infos}}\"); "
            "if [ -z \"$DBROOT\" ]; then echo \"[host] iPHoP DB not found under {params.db} (need a versioned subdir with db/ + db_infos/). Run: iphop download --db_dir {params.db}\" >&2; exit 1; fi; "
            "CUDA_VISIBLE_DEVICES=\"\" TF_FORCE_GPU_ALLOW_GROWTH=\"true\" iphop predict --fa_file {input} --db_dir \"$DBROOT\" --out_dir {output} --num_threads {threads}"


    rule host_merge:
        input:
            lambda wc: expand(f"{OUT}/opt_host/chunks/{{chunk}}", chunk=_host_chunk_ids(wc)),
        output:
            f"{OUT}/opt_host/iphop.tsv",
        params:
            outdir=f"{OUT}/opt_host/run",
            chunks_dir=f"{OUT}/opt_host/chunks",
            script=scr("iphop_normalize.py"),
        shell:
            # Concat per-chunk genus CSVs (single header via awk 'NR==1 || FNR>1')
            # into the merged run dir, then normalize like the single-run rule.
            # Zero chunks / no predictions -> iphop_normalize writes a header-only table.
            "mkdir -p {params.outdir}; "
            "CSVS=''; for d in {params.chunks_dir}/*/; do [ -d \"$d\" ] || continue; "
            "for f in \"$d\"/Host_prediction_to_genus_m*.csv; do [ -s \"$f\" ] && CSVS=\"$CSVS $f\"; done; done; "
            "if [ -n \"$CSVS\" ]; then awk 'NR==1 || FNR>1' $CSVS > {params.outdir}/Host_prediction_to_genus_m90.csv; fi; "
            "python {params.script} --indir {params.outdir} --out {output}"


else:

    rule host_iphop:
        input:
            f"{OUT}/02_cluster/votu_representatives.fna",
        output:
            f"{OUT}/opt_host/iphop.tsv",
        params:
            outdir=f"{OUT}/opt_host/run",
            db=lambda wildcards: f"{config['databases']['dir']}/{_H.get('db', 'iphop')}",
            script=scr("iphop_normalize.py"),
        threads: _H["threads"]
        resources:
            mem_mb=64000,
            runtime=2880,
        conda:
            env("iphop")
        shell:
            "mkdir -p {params.outdir}; "
            "DBROOT=$(d=$(find {params.db} -maxdepth 2 -type d -name db_infos 2>/dev/null | head -1); printf '%s' \"${{d%/db_infos}}\"); "
            "if [ -z \"$DBROOT\" ]; then echo \"[host] iPHoP DB not found under {params.db} (need a versioned subdir with db/ + db_infos/). Run: iphop download --db_dir {params.db}\" >&2; exit 1; fi; "
            "if [ ! -s {input} ] || ! grep -q '^>' {input}; then printf 'id\thost\thost_score\n' > {output}; exit 0; fi; "
            "CUDA_VISIBLE_DEVICES=\"\" TF_FORCE_GPU_ALLOW_GROWTH=\"true\" iphop predict --fa_file {input} --db_dir \"$DBROOT\" --out_dir {params.outdir} --num_threads {threads}; "
            "python {params.script} --indir {params.outdir} --out {output}"
