#!/usr/bin/env python3
"""Classify phage genes into auxiliary metabolic gene (AMG) classes - lightweight, DRAM-v-free.

Mirrors the DRAM-v annotate->distill->filter->categorize workflow without the ~150-200 GB DB:
  annotate : phold (structure-informed) re-annotates Pharokka hypothetical proteins; the
             existing Pharokka gene_families.tsv provides product/function for the rest.
             Both tables are merged before classification.
  distill  : this script maps each product to a curated AMG metabolic class and assigns a
             DRAM-v-style confidence tier `amg_score` (1=high, 2=moderate, 3=marginal).
  filter   : transposon/MGE machinery is excluded (DRAM-v is_transposon/--remove_transposons);
             only amg_score 1-3 candidates are emitted.
  categorize : per-gene calls are summarized per vOTU rep and per sample.

Reads the best available per-gene functional annotation (phold per_cds_predictions.tsv and/or
Pharokka gene_families.tsv), maps each gene's product/PHROG-category to a curated,
literature-grounded AMG metabolic-class scheme, and emits:

  --out-summary  amg_summary.tsv      per vOTU rep: id, amg_flag, n_amg_genes, n_amg_highconf, amg_classes
  --out-genes    amg_genes.tsv        per AMG gene: gene, contig, product, function, amg_class, amg_score, confidence
  --out-sample   amg_per_sample.tsv   per sample x AMG class: n_genes, n_votus, tpm

AMG class scheme (regex on product, case-insensitive; first match wins) is grounded in the
phage-AMG literature: Thompson et al. 2011 (Cell), Brum & Sullivan 2015, Hurwitz & U'Ren 2016,
Kieft et al. 2021 (Nat Commun, sulfur), Luo et al. 2022 / Warwick-Dugdale 2019 (reviews),
Enav/Mann P-acquisition work. Class I = core metabolism (photosynthesis, carbon, nucleotide);
Class II = peripheral (transport, cofactors). Genes in the PHROG "moron, auxiliary metabolic
gene and host takeover" category matching no specific class -> "Host takeover / other (PHROG moron)"
(amg_score 3). AMG calls are CANDIDATES: best practice is to flag here and validate by gene
neighbourhood / host-context before functional claims.

Column detection is fuzzy so both phold and Pharokka tables parse. contig is taken from the
annotation table when present, else mapped from --g2c (gene2genome.tsv: protein_id -> genome_id).
"""
import argparse
import os
import re
import sys

import pandas as pd

# Curated AMG metabolic classes. Ordered: specific/core classes first so a multi-match gene is
# assigned the most informative class. Patterns are case-insensitive regexes on the product.
AMG_CLASSES = {
    "Photosynthesis": [
        r"\bpsbA\b", r"\bpsbD\b", r"\bpsbE\b", r"photosystem", r"\bhli\b",
        r"high[- ]light[- ]inducible", r"plastocyanin", r"\bpetE\b", r"\bpetF\b",
        r"phycocyanin", r"phycoerythrin", r"phycobilisome", r"\bcpcA\b", r"\bcpcB\b",
        r"\bcpeT\b", r"protochlorophyllide", r"\bpcyA\b", r"chlorophyll",
    ],
    "Carbon metabolism / PPP": [
        r"transaldolase", r"\btalC\b", r"transketolase", r"\btktA?\b",
        r"glucose-6-phosphate dehydrogenase", r"\bzwf\b", r"6-phosphogluconate dehydrogenase",
        r"\bgnd\b", r"\bcp12\b", r"fructose-bisphosphate aldolase", r"\bfbaA?\b",
        r"glyceraldehyde-3-phosphate dehydrogenase", r"\bgapA\b", r"phosphoglucomutase",
        r"phosphofructokinase", r"\bpfkA\b", r"triosephosphate isomerase", r"\btpiA\b",
        r"ribulose-phosphate 3-epimerase", r"\brpe\b", r"mannose-6-phosphate isomerase",
        r"\bmanA\b", r"glucokinase", r"\bglk\b",
    ],
    "Nucleotide metabolism": [
        r"ribonucleotide reductase", r"ribonucleoside-diphosphate reductase", r"\bnrdA\b",
        r"\bnrdB\b", r"\bnrdE\b", r"\bnrdF\b", r"\bnrdH\b", r"\bnrdI\b", r"thymidylate synthase",
        r"\bthyX\b", r"\bthyA\b", r"\btd\b", r"dUTPase", r"deoxyuridine triphosphatase",
        r"\bdut\b", r"deoxycytidylate deaminase", r"\bdcd\b", r"dihydrofolate reductase",
        r"\bfolA\b", r"\bdfr\b", r"thioredoxin", r"glutaredoxin", r"guanylate kinase",
        r"\bgmk\b", r"cytidine deaminase", r"\bcdd\b", r"nucleoside diphosphate kinase",
        r"\bndk\b", r"adenylate kinase", r"\badk\b", r"uridylate kinase", r"\bpyrH\b",
    ],
    "Phosphate metabolism": [
        r"\bpstS\b", r"phosphate[- ]binding", r"\bphoH\b", r"alkaline phosphatase", r"\bphoA\b",
        r"\bphoU\b", r"\bphoB\b", r"\bphoR\b", r"polyphosphate kinase", r"\bppk\b",
        r"exopolyphosphatase", r"\bppx\b", r"phytase", r"glycerophosphoryl diester",
        r"\bugpQ\b", r"inorganic pyrophosphatase", r"\bppa\b",
    ],
    "Sulfur metabolism": [
        r"\bcysH\b", r"phosphoadenosine phosphosulfate reductase", r"PAPS reductase",
        r"\bcysC\b", r"\bcysD\b", r"\bcysN\b", r"\bcysI\b", r"\bcysJ\b", r"sulfite reductase",
        r"\bdsrA\b", r"\bdsrB\b", r"dissimilatory", r"\bsoxB\b", r"\bsoxY\b", r"\bsoxZ\b",
        r"\baprA\b", r"\baprB\b", r"adenylylsulfate", r"adenosine-5'-phosphosulfate", r"\bsat\b",
        r"taurine dioxygenase", r"\btauD\b", r"alkanesulfonate", r"\bssuD\b", r"\bsqr\b",
    ],
    "Nitrogen metabolism": [
        r"ammonium transporter", r"\bamtB?\b", r"nitrate reductase", r"\bnarG\b", r"\bnarH\b",
        r"\bnarI\b", r"\bnasA\b", r"nitrite reductase", r"\bnirA\b", r"\bnirK\b", r"nitrogenase",
        r"\bnifH\b", r"urease", r"\bureC\b", r"glutamine synthetase", r"\bglnA\b",
        r"glycine cleavage", r"\bgcvT\b", r"\bgcvH\b", r"\bgcvP\b", r"cyanate hydratase",
        r"\bcynS\b",
    ],
    "Carbohydrate (CAZy)": [
        r"glycoside hydrolase", r"glycosyl hydrolase", r"glycosylhydrolase", r"chitinase",
        r"chitosanase", r"cellulase", r"endoglucanase", r"\bamylase\b", r"pullulanase",
        r"alginate lyase", r"pectate lyase", r"pectin lyase", r"xylanase", r"beta-glucosidase",
        r"\bglucanase\b", r"mannosidase", r"galactosidase", r"fucosidase", r"sialidase",
        r"neuraminidase", r"hyaluronidase", r"hyaluronate lyase", r"agarase", r"carrageenase",
        r"chondroitin", r"keratan", r"heparin",
    ],
    "Peptidoglycan / cell wall": [
        r"\bmurA\b", r"\bmurB\b", r"\bmurC\b", r"\bmurD\b", r"\bmurE\b", r"\bmurF\b",
        r"\bmurG\b", r"peptidoglycan biosynth", r"UDP-N-acetylmuramoyl", r"\bamiA\b",
        r"\bdacA\b", r"D-alanyl-D-alanine", r"lipopolysaccharide biosynth", r"\blptD\b",
        r"O-antigen", r"teichoic acid",
    ],
    "Cofactors / vitamins": [
        r"cobalamin", r"\bcobS\b", r"\bcobT\b", r"\bcobQ\b", r"vitamin B12", r"queuosine",
        r"\bqueC\b", r"\bqueD\b", r"\bqueE\b", r"\bqueF\b", r"GTP cyclohydrolase", r"\bfolE\b",
        r"biotin synth", r"\bbioB\b", r"thiamine", r"\bthiC\b", r"\bthiD\b", r"\bthiE\b",
        r"molybdopterin", r"\bmoaA\b", r"\bmoeA\b", r"nicotinate", r"\bnadD\b", r"\bnadE\b",
        r"pyridoxal", r"\bpdxA\b", r"riboflavin", r"\bribA\b",
    ],
}
_CATCHALL = "Host takeover / other (PHROG moron)"
_MORON_RE = re.compile(r"moron|auxiliary metabolic|host takeover", re.I)
# Mobile-genetic-element / transposon machinery falsely matches metabolic DBs (DRAM-v
# is_transposon / --remove_transposons). Integrase/transposase are lysogeny/MGE markers, not
# AMGs, so exclude them from AMG calls.
_TRANSPOSON_RE = re.compile(
    r"transposas|transposon|insertion sequence|IS element|resolvase|recombinase|excisionase|"
    r"\bxis\b|conjugativ|mobilization|\bmobA\b|\bmobB\b|relaxase|invertase|"
    r"site-specific integrase|\bintegrase\b",
    re.I,
)
# Precompile: class -> [regex]
_CLASS_RES = {cls: [re.compile(p, re.I) for p in pats] for cls, pats in AMG_CLASSES.items()}


def classify(product: str, function: str, confidence: str = "") -> tuple[str, int]:
    """Return (amg_class, amg_score); ('', 0) if not an AMG candidate.

    amg_score mirrors DRAM-v auxiliary_score tiers: 1 = specific metabolic class + high phold
    confidence; 2 = specific class, lower confidence; 3 = PHROG moron/AMG/host-takeover only.
    """
    prod = (product or "").strip()
    if not prod or _TRANSPOSON_RE.search(prod):
        return "", 0
    for cls, res in _CLASS_RES.items():
        if any(r.search(prod) for r in res):
            return cls, (1 if (confidence or "").lower() == "high" else 2)
    if _MORON_RE.search(function or ""):
        return _CATCHALL, 3
    return "", 0


def _pick(df: pd.DataFrame, candidates) -> str:
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols:
            return cols[cand.lower()]
    # fuzzy contains
    for c in df.columns:
        cl = c.lower()
        if any(k in cl for k in candidates):
            return c
    return ""


def normalize_annotations(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """Normalize a phold or Pharokka table to a common gene/contig/product/function/confidence schema."""
    gene_c = _pick(df, ["gene", "cds_id", "gene_id", "cds", "locus_tag", "protein_id", "query", "id"])
    contig_c = _pick(df, ["contig", "seqid", "locus", "scaffold", "genome_id", "genome", "sequence"])
    prod_c = _pick(df, ["product", "annot", "annotation", "description"])
    func_c = _pick(df, ["function", "category", "phrog_category", "phrog_function"])
    conf_c = _pick(df, ["annotation_confidence", "confidence"])
    print(f"[amg] {label} cols detected: gene={gene_c!r} contig={contig_c!r} product={prod_c!r} "
          f"function={func_c!r} confidence={conf_c!r}", file=sys.stderr)
    out = pd.DataFrame()
    if not gene_c:
        return out
    out["gene"] = df[gene_c].astype(str)
    out["contig"] = df[contig_c].astype(str) if contig_c else ""
    out["product"] = df[prod_c].astype(str) if prod_c else ""
    out["function"] = df[func_c].astype(str) if func_c else ""
    out["confidence"] = df[conf_c].astype(str) if conf_c else ""
    out = out[out["gene"].str.strip() != ""]
    return out.reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ann", required=True, help="phold per_cds_predictions.tsv (hypothetical proteins)")
    ap.add_argument("--ann2", default="", help="Pharokka gene_families.tsv (PHROG-annotated genes)")
    ap.add_argument("--g2c", default="", help="gene2genome.tsv (protein_id -> genome_id) fallback contig map")
    ap.add_argument("--clusters", default="", help="votu_clusters.tsv (contig, votu, is_representative)")
    ap.add_argument("--c2s", default="", help="contig2sample.tsv")
    ap.add_argument("--coverm", default="", help="coverm.tsv (rep x sample TPM); optional")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-genes", required=True)
    ap.add_argument("--out-sample", required=True)
    a = ap.parse_args()

    ann = normalize_annotations(pd.read_csv(a.ann, sep="\t", dtype=str).fillna(""), "ann")
    if a.ann2 and os.path.exists(a.ann2):
        ann2 = normalize_annotations(pd.read_csv(a.ann2, sep="\t", dtype=str).fillna(""), "ann2")
        combined = pd.concat([ann, ann2], ignore_index=True).drop_duplicates(subset=["gene"], keep="first")
    else:
        combined = ann
    if combined.empty:
        for pth, hdr in [(a.out_summary, "id\tamg_flag\tn_amg_genes\tn_amg_highconf\tamg_classes\n"),
                         (a.out_genes, "gene\tcontig\tproduct\tfunction\tamg_class\tamg_score\tconfidence\n"),
                         (a.out_sample, "sample\tamg_class\tn_genes\tn_votus\ttpm\n")]:
            with open(pth, "w") as fh:
                fh.write(hdr)
        print("[amg] no annotation rows; wrote empty AMG tables", file=sys.stderr)
        return 0

    # gene -> contig map (fallback when the annotation lacks a contig column)
    g2c = {}
    if combined["contig"].eq("").all() and a.g2c and os.path.exists(a.g2c):
        g = pd.read_csv(a.g2c, sep="\t", dtype=str).fillna("")
        pc = _pick(g, ["protein_id", "gene", "gene_id", "id"])
        gc = _pick(g, ["genome_id", "contig", "genome"])
        if pc and gc:
            g2c = dict(zip(g[pc], g[gc]))

    rows = []
    for _, r in combined.iterrows():
        gene = r["gene"].strip()
        if not gene:
            continue
        contig = r["contig"] if r["contig"] else g2c.get(gene, "")
        product = r["product"]
        function = r["function"]
        conf = r["confidence"]
        cls, score = classify(product, function, conf)
        if cls:
            rows.append((gene, contig, product, function, cls, score, conf))

    genes = pd.DataFrame(
        rows, columns=["gene", "contig", "product", "function", "amg_class", "amg_score", "confidence"]
    )
    genes.to_csv(a.out_genes, sep="\t", index=False)

    # ---- per-vOTU summary (id = representative contig) ----
    if genes.empty:
        summary = pd.DataFrame(columns=["id", "amg_flag", "n_amg_genes", "n_amg_highconf", "amg_classes"])
    else:
        agg = genes.groupby("contig").agg(
            n_amg_genes=("gene", "size"),
            n_amg_highconf=("amg_score", lambda s: int((s <= 2).sum())),
            amg_classes=("amg_class", lambda s: "; ".join(sorted(set(s)))),
        ).reset_index().rename(columns={"contig": "id"})
        agg.insert(1, "amg_flag", 1)
        summary = agg[["id", "amg_flag", "n_amg_genes", "n_amg_highconf", "amg_classes"]]
    summary.to_csv(a.out_summary, sep="\t", index=False)

    # ---- per-sample x AMG-class ----
    per_sample = _per_sample(genes, a)
    per_sample.to_csv(a.out_sample, sep="\t", index=False)

    print(f"[amg] AMG genes={len(genes)}  vOTUs_with_AMG={genes['contig'].nunique() if len(genes) else 0}  "
          f"classes={genes['amg_class'].nunique() if len(genes) else 0}", file=sys.stderr)
    return 0


def _per_sample(genes: pd.DataFrame, a) -> pd.DataFrame:
    cols = ["sample", "amg_class", "n_genes", "n_votus", "tpm"]
    if genes.empty:
        return pd.DataFrame(columns=cols)

    # Prefer coverm (rep x sample TPM): consistent with ecology; presence = TPM>0.
    if a.coverm and os.path.exists(a.coverm) and os.path.getsize(a.coverm) > 0:
        cm = pd.read_csv(a.coverm, sep="\t", index_col=0)
        cm.index = cm.index.astype(str)
        g = genes[genes["contig"].isin(cm.index)].copy()
        if g.empty:
            return pd.DataFrame(columns=cols)
        recs = []
        for cls, sub in g.groupby("amg_class"):
            per_rep = sub.groupby("contig").size()  # genes of this class per rep
            tpm = cm.loc[per_rep.index]  # rep x sample TPM
            for sample in tpm.columns:
                vals = tpm[sample]
                present = vals[vals > 0]
                if present.empty:
                    continue
                recs.append((sample, cls, int(per_rep[present.index].sum()),
                             int(present.size), float(vals.sum())))
        if not recs:
            return pd.DataFrame(columns=cols)
        return pd.DataFrame(recs, columns=cols)

    # Fallback: membership presence via clusters + contig2sample (no abundance).
    if a.clusters and a.c2s and os.path.exists(a.clusters) and os.path.exists(a.c2s):
        cl = pd.read_csv(a.clusters, sep="\t", dtype=str).fillna("")
        c2s = pd.read_csv(a.c2s, sep="\t", dtype=str).fillna("")
        rep2votu = dict(zip(cl.loc[cl.get("is_representative", "0") == "1", "contig"],
                            cl.loc[cl.get("is_representative", "0") == "1", "votu"]))
        mem = cl.merge(c2s, on="contig", how="left")
        votu2samples = mem.groupby("votu")["sample"].apply(lambda s: set(s.dropna())).to_dict()
        g = genes.copy()
        g["votu"] = g["contig"].map(rep2votu)
        g = g.dropna(subset=["votu"])
        recs = []
        for cls, sub in g.groupby("amg_class"):
            cnt = sub.groupby("votu").size()
            for votu, n in cnt.items():
                for sample in votu2samples.get(votu, ()):
                    recs.append((sample, cls, int(n), 1, float("nan")))
        if not recs:
            return pd.DataFrame(columns=cols)
        df = pd.DataFrame(recs, columns=cols)
        return df.groupby(["sample", "amg_class"], as_index=False).agg(
            n_genes=("n_genes", "sum"), n_votus=("n_votus", "sum"), tpm=("tpm", "sum"))
    return pd.DataFrame(columns=cols)


if __name__ == "__main__":
    sys.exit(main())
