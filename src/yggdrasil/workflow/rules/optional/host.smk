# Optional iPHoP (simroux/iphop) host prediction. Gated by config host.enabled.
# iPHoP 1.4 CLI: `iphop predict --fa_file <fna> --db_dir <db> --out_dir <dir> --num_threads <n>`.
# iPHoP unpacks into a VERSIONED subdir (<db_dir>/iPHoP_db_<ver>_rw/ or Test_db_rw_v1.4/) that holds
# the db/ + db_infos/ folders. The rule auto-resolves that subdir under <databases.dir>/<host.db>, so
# the config path stays stable across DB versions.
# The genus-level CSV (Host_prediction_to_genus_m<score>.csv) is normalized to id,host,host_score
# by scripts/iphop_normalize.py so report.smk's votu_master can join on representative id.
# iPHoP's final classifier uses TensorFlow and grabs the GPU by default; we force CPU (small-VRAM
# cards OOM otherwise) via CUDA_VISIBLE_DEVICES="".

rule host_iphop:
    input:
        f"{OUT}/02_cluster/votu_representatives.fna",
    output:
        f"{OUT}/opt_host/iphop.tsv",
    params:
        outdir=f"{OUT}/opt_host/run",
        db=lambda wildcards: f"{config['databases']['dir']}/{config['host'].get('db', 'iphop')}",
        script=scr("iphop_normalize.py"),
    threads: config["host"]["threads"]
    conda:
        env("iphop")
    shell:
        "mkdir -p {params.outdir}; "
        "DBROOT=$(d=$(find {params.db} -maxdepth 2 -type d -name db_infos 2>/dev/null | head -1); printf '%s' \"${{d%/db_infos}}\"); "
        "if [ -z \"$DBROOT\" ]; then echo \"[host] iPHoP DB not found under {params.db} (need a versioned subdir with db/ + db_infos/). Run: iphop download --db_dir {params.db}\" >&2; exit 1; fi; "
        "if [ ! -s {input} ] || ! grep -q '^>' {input}; then printf 'id\thost\thost_score\n' > {output}; exit 0; fi; "
        "CUDA_VISIBLE_DEVICES=\"\" TF_FORCE_GPU_ALLOW_GROWTH=\"true\" iphop predict --fa_file {input} --db_dir \"$DBROOT\" --out_dir {params.outdir} --num_threads {threads}; "
        "python {params.script} --indir {params.outdir} --out {output}"
