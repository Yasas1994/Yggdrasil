"""Smoke tests for the spine: sample discovery + pure-python helper scripts.

These run without any external bioinformatics tool or database, exercising the
load-bearing logic that the delegated rules feed into.
"""
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import pytest

from yggdrasil.cli import discover

SCR = Path(str(files("yggdrasil") / "workflow" / "scripts"))


def run(name, *args):
    return subprocess.run(
        [sys.executable, str(SCR / name), *args],
        check=True, capture_output=True, text=True,
    )


def write(p, s):
    p.write_text(s)
    return p


def test_discover_pairs_mobile_suffix(tmp_path):
    g, r = tmp_path / "g", tmp_path / "r"
    g.mkdir(); r.mkdir()
    write(g / "Kenneth004_mobile.fna", ">c1\nAAAA\n")
    write(r / "Kenneth004_R1.fastq.gz", "@r1\nAAAA\n+\n!!!!\n")
    write(r / "Kenneth004_R2.fastq.gz", "@r2\nAAAA\n+\n!!!!\n")
    rows = discover(g, r)
    assert len(rows) == 1
    assert rows[0]["sample"] == "Kenneth004"
    assert rows[0]["r1"].endswith("Kenneth004_R1.fastq.gz")
    assert rows[0]["r2"].endswith("Kenneth004_R2.fastq.gz")


def test_discover_handles_no_reads(tmp_path):
    g = tmp_path / "g"; g.mkdir()
    write(g / "S1.fna", ">c1\nAAAA\n")
    rows = discover(g, None)
    assert rows[0]["sample"] == "S1" and rows[0]["r1"] == ""


def test_spine_scripts_end_to_end(tmp_path):
    # --- fixtures: 2 samples, contigs of mixed length + one exact duplicate ---
    g = tmp_path / "g"; g.mkdir()
    longA, longB = "A" * 6000, "C" * 6000
    write(g / "S1.fna", f">contig1\n{longA}\n>short\n{'T'*100}\n")
    write(g / "S2.fna", f">contig9\n{longB}\n>dup\n{longA}\n")  # dup == S1 contig1
    samples = tmp_path / "samples.tsv"
    write(samples, "sample\tgenome\tr1\tr2\n"
                   f"S1\t{g/'S1.fna'}\t\t\nS2\t{g/'S2.fna'}\t\t\n")

    # preprocess (dedup on to drop the cross-sample duplicate)
    fna, cmap, lengths = (tmp_path / n for n in ("all.fna", "c2s.tsv", "len.tsv"))
    run("preprocess.py", "--samples", str(samples), "--min-length", "5000",
        "--dedup", "--out-fna", str(fna), "--out-map", str(cmap), "--out-lengths", str(lengths))
    ids = [l[1:] for l in fna.read_text().splitlines() if l.startswith(">")]
    assert "S1__contig1" in ids and "S2__contig9" in ids
    assert not any("short" in i for i in ids)          # length-filtered
    assert len(ids) == 2                               # duplicate dropped

    # fake vclust clusters (member, cluster) + CheckV quality
    cl = tmp_path / "clusters.tsv"
    write(cl, "member\tcluster\nS1__contig1\tvOTU_1\nS2__contig9\tvOTU_1\n")
    qual = tmp_path / "quality_summary.tsv"
    write(qual, "contig_id\tcheckv_quality\tcompleteness\tcontamination\n"
                "S1__contig1\tHigh-quality\t95\t0\nS2__contig9\tComplete\t100\t0\n")

    # select_representatives: S2__contig9 is Complete -> representative
    vc, reps = tmp_path / "votu_clusters.tsv", tmp_path / "reps.fna"
    run("select_representatives.py", "--clusters", str(cl), "--quality", str(qual),
        "--lengths", str(lengths), "--fasta", str(fna),
        "--out-clusters", str(vc), "--out-fasta", str(reps))
    assert "S2__contig9" in reps.read_text()
    vcl = vc.read_text()
    assert "S2__contig9\tvOTU_1\t1" in vcl
    assert "S1__contig1\tvOTU_1\t0" in vcl

    # build_count_tables: gene_families with one family hit on S1 contig
    gf = tmp_path / "gene_families.tsv"
    write(gf, "contig\tgene_id\tfamily\nS1__contig1\tg1\tPHROG_0001\n"
              "S2__contig9\tg2\tPHROG_0001\nS2__contig9\tg3\tPHROG_0002\n")
    votu_mat, gene_mat = tmp_path / "votu.tsv", tmp_path / "gene.tsv"
    run("build_count_tables.py", "--contig2sample", str(cmap), "--votu-clusters", str(vc),
        "--gene-families", str(gf), "--samples", str(samples),
        "--out-votu", str(votu_mat), "--out-gene", str(gene_mat))
    vh = votu_mat.read_text().splitlines()
    assert vh[0].split("\t") == ["votu", "S1", "S2"]
    assert vh[1].split("\t") == ["vOTU_1", "1", "1"]   # one contig in each sample
    gh = gene_mat.read_text().splitlines()
    rows = {r.split("\t")[0]: r.split("\t")[1:] for r in gh[1:]}
    assert rows["PHROG_0001"] == ["1", "1"]
    assert rows["PHROG_0002"] == ["0", "1"]

    # votu_master: join with taxonomy + lifestyle (host/amg absent -> NA)
    tax = tmp_path / "taxonomy.tsv"; write(tax, "id\ttaxonomy\nS2__contig9\tCaudoviricetes\n")
    life = tmp_path / "lifestyle.tsv"; write(life, "id\tlifestyle\tlifestyle_score\nS2__contig9\tvirulent\t0.97\n")
    master = tmp_path / "master.tsv"
    run("votu_master.py", "--votu-clusters", str(vc), "--contig2sample", str(cmap),
        "--lengths", str(lengths), "--quality", str(qual), "--taxonomy", str(tax),
        "--lifestyle", str(life), "--out", str(master))
    header = master.read_text().splitlines()[0].split("\t")
    row = master.read_text().splitlines()[1].split("\t")
    d = dict(zip(header, row))
    assert d["representative"] == "S2__contig9"
    assert d["n_members"] == "2" and d["n_samples"] == "2"
    assert d["checkv_quality"] == "Complete"
    assert d["taxonomy"] == "Caudoviricetes"
    assert d["lifestyle"] == "virulent"
