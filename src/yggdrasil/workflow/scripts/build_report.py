#!/usr/bin/env python3
"""Build the Yggdrasil phage report using bioinformatics-report-builder.

Inputs:
  --master        vOTU master TSV
  --votu          vOTU x sample count TSV
  --gene          gene family x sample count TSV
  --alpha         alpha diversity TSV
  --vc-graphml    vContact3 network graphml
  --vc-refs       vContact3 final_assignments.csv
  --vclust        vclust cluster TSV (optional, for sample vOTU list)
  --out           output HTML path
  --output-dir    directory to save QMD/assets (defaults to parent of --out)

The script writes a .qmd under output_dir, renders it with Quarto, and
creates an assets/ folder containing figures and linked resources.
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from bioinformatics_report import Report


# --------------------------------------------------------------------------- #
# Network helpers
# --------------------------------------------------------------------------- #

def _load_taxonomy(ref_path: Path) -> dict[str, dict[str, str]]:
    """Return Genome -> taxonomic assignments from vContact3 final_assignments."""
    taxa: dict[str, dict[str, str]] = {}
    if not ref_path.exists():
        return taxa
    df = pd.read_csv(ref_path)
    for _, row in df.iterrows():
        genome = str(row.get("Genome", ""))
        if not genome:
            continue
        taxa[genome] = {
            "GenomeName": str(row.get("GenomeName", "")),
            "Reference": bool(row.get("Reference", True)),
            "family": str(row.get("family_prediction", "")),
            "genus": str(row.get("genus_prediction", "")),
            "order": str(row.get("order_prediction", "")),
        }
    return taxa


def _sample_nodes(taxa: dict[str, dict[str, str]]) -> set[str]:
    return {g for g, d in taxa.items() if not d.get("Reference", True)}


def _taxonomy_label(d: dict[str, str]) -> str:
    """Use the lowest non-novel, non-empty prediction as the node label."""
    for key in ("family", "genus", "order"):
        val = d.get(key, "")
        if val and not pd.isna(val) and val.lower() not in {"nan", "none", ""}:
            return val
    return "Unclassified"


def _build_network_figure(
    graphml_path: Path,
    ref_path: Path,
    out_png: Path,
    min_shared_genes: int = 10,
    top_edges_per_sample: int = 8,
    max_refs: int = 150,
    figsize: tuple[int, int] = (7, 7),
    dpi: int = 150,
) -> None:
    """Create a filtered vContact3 gene-sharing network figure.

    - Keep all sample vOTU nodes.
    - For each sample vOTU, keep the strongest edges (shared_genes) up to
      top_edges_per_sample and above min_shared_genes.
    - Keep the top max_refs reference neighbours by cumulative weight.
    - Color every node by its lowest resolved taxonomy (family/genus/order).
    - Mark sample vOTUs with a black edge and larger size.
    """
    if not graphml_path.exists() or not ref_path.exists():
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.text(0.5, 0.5, "vContact3 network not available", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return

    taxa = _load_taxonomy(ref_path)
    samples = _sample_nodes(taxa)

    G = nx.read_graphml(graphml_path)

    # Gather incident edges for each sample node, sorted by shared_genes.
    sample_edges: dict[str, list[tuple[str, int]]] = {s: [] for s in samples}
    for u, v, d in G.edges(data=True):
        w = d.get("shared_genes", 0)
        if w < min_shared_genes:
            continue
        if u in samples:
            sample_edges[u].append((v, w))
        if v in samples:
            sample_edges[v].append((u, w))

    selected_edges: set[tuple[str, str]] = set()
    ref_weight: Counter[str] = Counter()
    for s, edges in sample_edges.items():
        edges.sort(key=lambda x: x[1], reverse=True)
        for ref, w in edges[:top_edges_per_sample]:
            selected_edges.add(tuple(sorted((s, ref))))
            ref_weight[ref] += w

    # Also include edges between sample vOTUs when they share genes.
    for s1, s2 in [(u, v) for u, v, d in G.edges(data=True) if u in samples and v in samples]:
        w = G[s1][s2].get("shared_genes", 0)
        if w >= min_shared_genes:
            selected_edges.add(tuple(sorted((s1, s2))))

    # Keep top reference neighbours and all sample nodes.
    top_refs = {n for n, _ in ref_weight.most_common(max_refs)}
    kept_nodes = samples | top_refs

    H = nx.Graph()
    H.add_nodes_from(kept_nodes)
    for u, v in selected_edges:
        if u in kept_nodes and v in kept_nodes:
            H.add_edge(u, v)

    # Drop isolated references; always keep sample nodes.
    H.remove_nodes_from([n for n in list(H.nodes()) if n not in samples and H.degree(n) == 0])

    n_refs = H.number_of_nodes() - len(samples)
    n_refs = max(n_refs, 0)

    if H.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.text(0.5, 0.5, "No network edges above threshold", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return

    # Taxonomic color scheme.
    labels = {_taxonomy_label(taxa.get(n, {})): n for n in H.nodes()}
    labels = {}
    for n in H.nodes():
        labels[n] = _taxonomy_label(taxa.get(n, {}))

    taxa_counts = Counter(labels.values())
    sorted_taxa = sorted(taxa_counts, key=lambda x: taxa_counts[x], reverse=True)
    palette = (
        plt.cm.tab20.colors
        if len(sorted_taxa) <= 20
        else plt.cm.tab20.colors + plt.cm.tab20b.colors
    )
    color_map = {tax: palette[i % len(palette)] for i, tax in enumerate(sorted_taxa)}

    node_color = [color_map[labels[n]] for n in H.nodes()]
    node_size = [450 if n in samples else 160 for n in H.nodes()]
    edgecolors = ["black" if n in samples else "#888888" for n in H.nodes()]
    linewidths = [2.5 if n in samples else 0.5 for n in H.nodes()]

    pos = nx.spring_layout(H, k=0.45, iterations=80, seed=42)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.3, edge_color="#999999", width=0.5)
    nx.draw_networkx_nodes(
        H,
        pos,
        ax=ax,
        node_color=node_color,
        node_size=node_size,
        edgecolors=edgecolors,
        linewidths=linewidths,
        alpha=0.95,
    )

    # Legend: top taxa + sample marker (truncate long labels).
    max_label = 45
    handles = []
    for tax in sorted_taxa[:15]:
        label = tax if len(tax) <= max_label else tax[:max_label - 3] + "..."
        handles.append(mpatches.Patch(color=color_map[tax], label=label))
    handles.append(
        mpatches.Patch(
            facecolor="#ffffff", edgecolor="black", linewidth=1.5, label="sample vOTU"
        )
    )
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize="small",
        title="Taxonomy",
    )
    ax.set_title(
        f"vContact3 gene-sharing network (edges with shared genes ≥ {min_shared_genes}; "
        f"{len(samples)} sample vOTUs, {n_refs} reference neighbours)"
    )
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Optional-module figures (VIRIDIC / VirClust)
# --------------------------------------------------------------------------- #

def _build_viridic_heatmap(
    sim_path: Path,
    out_png: Path,
    figsize: tuple[int, int] = (10, 9),
    dpi: int = 150,
) -> None:
    """Render the VIRIDIC intergenomic-similarity matrix as a heatmap."""
    if not sim_path.exists() or os.path.getsize(sim_path) == 0:
        return
    df = pd.read_csv(sim_path, sep="\t", index_col=0)
    if df.empty:
        return
    n = len(df)
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(df.values, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    if n <= 60:
        ax.set_xticks(range(n))
        ax.set_xticklabels(df.columns, rotation=90, fontsize=5)
        ax.set_yticks(range(n))
        ax.set_yticklabels(df.index, fontsize=5)
        ax.tick_params(length=0)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="intergenomic similarity (%)")
    ax.set_title(f"VIRIDIC intergenomic similarity ({n} vOTUs)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _parse_newick(text: str) -> dict:
    """Minimal Newick parser -> nested {name, length, children} dicts."""
    text = text.strip().rstrip(";")
    pos = [0]

    def parse() -> dict:
        node: dict = {"name": "", "length": 0.0, "children": []}
        if pos[0] < len(text) and text[pos[0]] == "(":
            pos[0] += 1
            while pos[0] < len(text) and text[pos[0]] != ")":
                node["children"].append(parse())
                if pos[0] < len(text) and text[pos[0]] == ",":
                    pos[0] += 1
            pos[0] += 1  # ')'
        start = pos[0]
        while pos[0] < len(text) and text[pos[0]] not in ",():":
            pos[0] += 1
        node["name"] = text[start:pos[0]].strip().strip("'\"")
        if pos[0] < len(text) and text[pos[0]] == ":":
            pos[0] += 1
            start = pos[0]
            while pos[0] < len(text) and text[pos[0]] not in ",()":
                pos[0] += 1
            try:
                node["length"] = float(text[start:pos[0]])
            except ValueError:
                node["length"] = 0.0
        return node

    return parse()


def _build_virclust_tree(
    tree_path: Path,
    out_png: Path,
    figsize: tuple[int, int] = (9, 12),
    dpi: int = 150,
) -> None:
    """Render the VirClust hierarchical tree (Newick) as a rectangular phylogram."""
    if not tree_path.exists() or os.path.getsize(tree_path) == 0:
        return
    try:
        root = _parse_newick(tree_path.read_text())
    except Exception:
        return
    leaves: list[dict] = []

    def assign(node: dict, depth: float) -> float:
        node["x"] = depth
        if not node["children"]:
            node["y"] = float(len(leaves))
            leaves.append(node)
            return node["y"]
        ys = [assign(c, depth + c["length"]) for c in node["children"]]
        node["y"] = sum(ys) / len(ys)
        return node["y"]

    assign(root, 0.0)
    if not leaves:
        return

    fig, ax = plt.subplots(figsize=figsize)

    def draw(node: dict) -> None:
        if node["children"]:
            ys = [c["y"] for c in node["children"]]
            ax.plot([node["x"], node["x"]], [min(ys), max(ys)], color="black", lw=0.8)
            for c in node["children"]:
                ax.plot([node["x"], c["x"]], [c["y"], c["y"]], color="black", lw=0.8)
                draw(c)

    draw(root)
    xpad = max((l["x"] for l in leaves), default=1.0) * 0.02 + 1e-6
    for leaf in leaves:
        ax.text(leaf["x"] + xpad, leaf["y"], leaf["name"], fontsize=6, va="center")
    ax.set_ylim(len(leaves) - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xlabel("intergenomic distance (shared protein content)")
    ax.set_title(f"VirClust hierarchical clustering ({len(leaves)} vOTUs)")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Report sections
# --------------------------------------------------------------------------- #

def _build_amg_class_figure(genes_path: Path, out_path: Path) -> None:
    """Horizontal bar chart of AMG genes per metabolic class, split by confidence tier."""
    df = pd.read_csv(genes_path, sep="\t")
    if df.empty or "amg_class" not in df.columns:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No AMGs detected", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return
    order = df["amg_class"].value_counts().index.tolist()
    score = df["amg_score"] if "amg_score" in df.columns else pd.Series(2, index=df.index)
    hi = df[score <= 2].groupby("amg_class").size().reindex(order).fillna(0)
    lo = df[score > 2].groupby("amg_class").size().reindex(order).fillna(0)
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.42 * len(order))))
    ax.barh(order, hi, color="#2a7f4f", label="high/mod confidence (score 1-2)")
    ax.barh(order, lo, left=hi, color="#c7a008", label="moron/host-takeover (score 3)")
    ax.invert_yaxis()
    ax.set_xlabel("AMG genes")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _table_rows(path: Path, max_rows: int = 2000) -> tuple[list[str], list[list[str]]]:
    """Return headers and rows for the report builder add_table()."""
    if not path.exists() or os.path.getsize(path) == 0:
        return ([], [])
    df = pd.read_csv(path, sep="\t")
    df = df.head(max_rows)
    return list(df.columns), [[str(v) for v in row] for row in df.values]


def _metrics_from_master(master: Path) -> list[dict]:
    df = pd.read_csv(master, sep="\t")
    metrics = [
        {"label": "vOTUs", "value": str(len(df)), "sub": "representative genomes", "highlight": True},
    ]
    # Samples from the counts matrix if available.
    return metrics


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master", required=True, type=Path)
    ap.add_argument("--votu", required=True, type=Path)
    ap.add_argument("--gene", required=True, type=Path)
    ap.add_argument("--alpha", required=True, type=Path)
    ap.add_argument("--vc-graphml", required=True, type=Path)
    ap.add_argument("--vc-refs", required=True, type=Path)
    ap.add_argument("--viridic-sim", type=Path, default=None)
    ap.add_argument("--viridic-clusters", type=Path, default=None)
    ap.add_argument("--virclust-tree", type=Path, default=None)
    ap.add_argument("--virclust-clusters", type=Path, default=None)
    ap.add_argument("--amg-votu", type=Path, default=None)
    ap.add_argument("--amg-genes", type=Path, default=None)
    ap.add_argument("--amg-sample", type=Path, default=None)
    ap.add_argument("--output-dir", type=Path, default=None)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    output_dir = args.output_dir or args.out.parent
    assets_dir = output_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)

    # Build network figure.
    network_png = assets_dir / "vcontact_network.png"
    _build_network_figure(
        args.vc_graphml,
        args.vc_refs,
        network_png,
        min_shared_genes=10,
        top_edges_per_sample=8,
        max_refs=150,
    )

    # Master table metrics.
    master_df = pd.read_csv(args.master, sep="\t")
    votu_df = pd.read_csv(args.votu, sep="\t")
    samples = list(votu_df.columns[1:]) if len(votu_df.columns) > 1 else []
    classified = (master_df.get("family", pd.Series(dtype=str)).fillna("") != "").sum()

    report = Report(
        title_line1="Yggdrasil",
        title_line2="Phage Downstream Analysis Report",
        pipeline="vclust · CheckV · vContact3 · PhaStyle · BACPHLIP",
        operator="Yggdrasil pipeline",
        reference="INPHARED / PHROG",
        cluster="",
        run_label="phage_pipe",
        run_status="Pipeline complete",
        sidebar_footer="Yggdrasil phage pipeline",
        footer_left="Yggdrasil phage pipeline",
        footer_right="Generated by bioinformatics-report-builder",
        page_width="1000px",
    )
    report.set_metadata(
        [
            ("Samples", str(len(samples))),
            ("vOTUs", str(len(master_df))),
            ("Classified", f"{classified} ({100*classified/len(master_df):.1f}%)"),
        ]
    )

    summary = report.add_section("01", "Executive Summary", count=f"{len(samples)} samples · {len(master_df)} vOTUs")
    # Inject CSS so long tables scroll vertically and keep a sticky header.
    summary.add_raw(
        "<style>"
        ".data-table{display:block;max-height:55vh;overflow-y:auto;overflow-x:auto;width:100%;border-collapse:collapse}"
        ".data-table thead th{position:sticky;top:0;background:var(--panel,#f3f3f3);z-index:2}"
        "</style>"
    )
    summary.set_overview(
        "Quality-controlled phage genomes were dereplicated into vOTUs, taxonomically "
        "classified with vContact3, and annotated for lifestyle and abundance."
    )
    summary.add_metrics(
        [
            {"label": "vOTUs", "value": str(len(master_df)), "sub": "representative genomes", "highlight": True},
            {"label": "Samples", "value": str(len(samples)), "sub": "quantified"},
            {"label": "Classified", "value": f"{100*classified/len(master_df):.1f}%", "sub": f"{classified} with family-level taxonomy"},
        ]
    )

    # vOTU master table.
    sec_master = report.add_section("02", "vOTU Master Table", count=f"{len(master_df)} rows")
    headers, rows = _table_rows(args.master, max_rows=1000)
    sec_master.add_table(
        headers=headers,
        rows=rows,
        caption="vOTU summary with quality, taxonomy, lifestyle and optional host/AMG calls.",
        sortable=True,
    )

    # Network section.
    sec_net = report.add_section("03", "vContact3 Network", count="gene-sharing network")
    sec_net.set_overview(
        "Gene-sharing network from vContact3. Reference genomes are shown as grey-filled "
        "nodes; sample-derived vOTUs are outlined in black and coloured by their lowest "
        "resolved taxonomy (family / genus / order)."
    )
    sec_net.add_figure(network_png, caption="vContact3 gene-sharing network coloured by taxonomy. Sample vOTUs are highlighted with black borders.", width="45%")

    # Alpha diversity.
    sec_alpha = report.add_section("04", "Alpha Diversity", count="Shannon / Simpson / Richness")
    headers, rows = _table_rows(args.alpha, max_rows=1000)
    sec_alpha.add_table(
        headers=headers,
        rows=rows,
        caption="Alpha diversity indices per sample.",
        sortable=True,
    )

    # vOTU presence/absence.
    sec_votu = report.add_section("05", "vOTU × Sample Presence/Absence", count="presence/absence matrix")
    headers, rows = _table_rows(args.votu, max_rows=2000)
    sec_votu.add_table(
        headers=headers,
        rows=rows,
        caption="vOTU detection (1 = present) across samples. Ecological metrics use CoverM abundance, not this matrix.",
        sortable=True,
    )

    sec_gene = report.add_section("06", "Gene Family × Sample Counts", count="functional matrix")
    headers, rows = _table_rows(args.gene, max_rows=2000)
    sec_gene.add_table(
        headers=headers,
        rows=rows,
        caption="PHROG gene family counts across samples.",
        sortable=True,
    )

    # Lifestyle: PhaStyle vs BACPHLIP.
    _life = master_df.get("lifestyle", pd.Series(dtype=str)).fillna("")
    _bpl = master_df.get("bacphlip_lifestyle", pd.Series(dtype=str)).fillna("")
    _pt = int((_life == "temperate").sum()); _pv = int((_life == "virulent").sum())
    _bt = int((_bpl == "temperate").sum()); _bv = int((_bpl == "virulent").sum())
    _both = _life.isin(["temperate", "virulent"]) & _bpl.isin(["temperate", "virulent"])
    _agree = int((_life[_both] == _bpl[_both]).sum()); _nboth = int(_both.sum())
    _agree_pct = (100.0 * _agree / _nboth) if _nboth else 0.0
    sec_life = report.add_section("07", "Lifestyle (PhaStyle vs BACPHLIP)", count="temperate / virulent")
    sec_life.set_overview(
        "Two independent lifestyle predictors: PhaStyle (ProkBERT) and BACPHLIP "
        "(conserved-domain random forest). BACPHLIP assumes complete phage genomes "
        "and defaults to virulent when no lysogeny domains are found; agreement is "
        "reported over vOTUs called by both methods."
    )
    sec_life.add_metrics(
        [
            {"label": "PhaStyle temperate", "value": str(_pt), "sub": f"{_pv} virulent"},
            {"label": "BACPHLIP temperate", "value": str(_bt), "sub": f"{_bv} virulent"},
            {"label": "Agreement", "value": f"{_agree_pct:.1f}%", "sub": f"{_agree}/{_nboth} both-called", "highlight": True},
        ]
    )
    sec_life.add_table(
        headers=["method", "temperate", "virulent", "total_called"],
        rows=[
            ["PhaStyle", _pt, _pv, _pt + _pv],
            ["BACPHLIP", _bt, _bv, _bt + _bv],
            ["both (agree)", _agree, "", _nboth],
        ],
        caption="Lifestyle calls per method over vOTU representatives.",
        sortable=False,
    )

    # AMG (auxiliary metabolic genes) - optional; phold + curated AMG-class mapping.
    if args.amg_genes and args.amg_genes.exists() and os.path.getsize(args.amg_genes) > 0:
        amg_df = pd.read_csv(args.amg_genes, sep="\t")
        n_amg = len(amg_df)
        n_votu_amg = amg_df["contig"].nunique() if "contig" in amg_df.columns else 0
        n_hi = int((amg_df["amg_score"] <= 2).sum()) if "amg_score" in amg_df.columns else n_amg
        n_cls = amg_df["amg_class"].nunique() if "amg_class" in amg_df.columns else 0
        amg_png = assets_dir / "amg_classes.png"
        _build_amg_class_figure(args.amg_genes, amg_png)
        sec_amg = report.add_section("08", "Auxiliary Metabolic Genes", count=f"{n_amg} AMG genes · {n_votu_amg} vOTUs")
        sec_amg.set_overview(
            "Auxiliary metabolic genes (AMGs) identified without DRAM-v (its ~150-200 GB DB is "
            "impractical here). AMG calls are derived from PHROG product/function annotations; "
            "when the optional phold module is enabled, Pharokka hypothetical proteins are further "
            "re-annotated by structural homology (ProstT5 + Foldseek, ~16 GB DB) and merged with "
            "the PHROG calls. Each gene's product/function is mapped to a curated, literature-grounded "
            "AMG class scheme. amg_score mirrors DRAM-v auxiliary_score (1=high, 2=moderate, "
            "3=PHROG moron/host-takeover); transposon/MGE genes are excluded. Calls are candidates "
            "- validate by genomic context before functional claims."
        )
        sec_amg.add_metrics(
            [
                {"label": "AMG genes", "value": str(n_amg), "sub": f"{n_hi} high/mod confidence", "highlight": True},
                {"label": "vOTUs with AMG", "value": str(n_votu_amg), "sub": f"of {len(master_df)} total"},
                {"label": "AMG classes", "value": str(n_cls), "sub": "metabolic categories"},
            ]
        )
        sec_amg.add_figure(amg_png, caption="AMG genes per metabolic class, coloured by confidence tier (DRAM-v-style amg_score).", width="60%")
        if args.amg_sample and args.amg_sample.exists() and os.path.getsize(args.amg_sample) > 0:
            headers, rows = _table_rows(args.amg_sample, max_rows=2000)
            sec_amg.add_table(headers=headers, rows=rows, caption="Per-sample × AMG-class summary (n_genes, n_votus, TPM abundance).", sortable=True)
        if args.amg_votu and args.amg_votu.exists() and os.path.getsize(args.amg_votu) > 0:
            headers, rows = _table_rows(args.amg_votu, max_rows=1000)
            sec_amg.add_table(headers=headers, rows=rows, caption="Per-vOTU AMG summary (n_amg_genes, n_amg_highconf, amg_classes).", sortable=True)

    # VIRIDIC intergenomic similarity (optional).
    if args.viridic_sim and args.viridic_sim.exists() and os.path.getsize(args.viridic_sim) > 0:
        viridic_png = assets_dir / "viridic_heatmap.png"
        _build_viridic_heatmap(args.viridic_sim, viridic_png)
        sec_vd = report.add_section("09", "VIRIDIC Intergenomic Similarity", count="ICTV-style distances")
        sec_vd.set_overview(
            "Pairwise intergenomic similarities between vOTU representatives computed "
            "with VIRIDIC (BLASTN-based, ICTV algorithm). The cluster table reports "
            "the species- and genus-level assignments at the configured similarity "
            "thresholds."
        )
        sec_vd.add_figure(viridic_png, caption="VIRIDIC intergenomic similarity heatmap.", width="55%")
        if args.viridic_clusters and args.viridic_clusters.exists():
            headers, rows = _table_rows(args.viridic_clusters, max_rows=2000)
            sec_vd.add_table(
                headers=headers,
                rows=rows,
                caption="VIRIDIC species/genus cluster assignments per vOTU.",
                sortable=True,
            )

    # VirClust hierarchical clustering (optional).
    if args.virclust_tree and args.virclust_tree.exists() and os.path.getsize(args.virclust_tree) > 0:
        virclust_png = assets_dir / "virclust_tree.png"
        _build_virclust_tree(args.virclust_tree, virclust_png)
        sec_vc = report.add_section("10", "VirClust Hierarchical Clustering", count="protein-based tree")
        sec_vc.set_overview(
            "Hierarchical clustering of vOTU representatives from shared protein "
            "content (VirClust, BLASTP protein clusters). Branch lengths are "
            "intergenomic distances; the table maps each vOTU to its viral genome "
            "cluster (VGC)."
        )
        sec_vc.add_figure(virclust_png, caption="VirClust protein-based hierarchical tree of vOTUs.", width="55%")
        if args.virclust_clusters and args.virclust_clusters.exists():
            headers, rows = _table_rows(args.virclust_clusters, max_rows=2000)
            sec_vc.add_table(
                headers=headers,
                rows=rows,
                caption="VirClust viral genome cluster (VGC) assignment and per-genome protein-sharing statistics.",
                sortable=True,
            )

    # Save QMD and render.
    qmd_path = report.save(output_dir, filename="phage_report.qmd")
    subprocess.run(["quarto", "render", str(qmd_path), "--to", "html"], check=True)

    # Quarto renders to <qmd stem>.html in the same directory.
    rendered_html = qmd_path.with_suffix(".html")
    if rendered_html != args.out:
        rendered_html.rename(args.out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
