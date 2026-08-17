#!/usr/bin/env Rscript
# Ecological summaries for Yggdrasil: alpha diversity, Bray-Curtis beta,
# PCoA ordination, and (optionally) group-wise differential abundance.
# Input is the vOTU x sample ABUNDANCE matrix (CoverM TPM); counts are
# not meaningful here. Base R + vegan only. ANCOM-BC (Bioconductor) is used if
# installed, else a transparent Wilcoxon/Kruskal fallback with BH adjustment.
suppressPackageStartupMessages(library(vegan))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NA) {
  i <- match(key, args)
  if (is.na(i) || i == length(args)) return(default)
  v <- args[i + 1]
  if (startsWith(v, "--")) return(default)  # value omitted; next token is a flag
  v
}

counts_path <- get_arg("--abundance")
if (is.na(counts_path)) counts_path <- get_arg("--counts")
samples_path <- get_arg("--samples")
group_col   <- get_arg("--group-col", NA)
out_alpha   <- get_arg("--out-alpha", NA)
out_beta    <- get_arg("--out-beta", NA)
out_ord     <- get_arg("--out-ordination", NA)
out_diff    <- get_arg("--out-differential", NA)

if (is.na(group_col) || group_col %in% c("", "NA", "null", "None")) group_col <- NA

# ---- abundance: votu (rows) x sample (cols), first column = id ----
cm <- as.matrix(read.table(counts_path, header = TRUE, sep = "\t",
                           row.names = 1, check.names = FALSE))
cm[is.na(cm)] <- 0
storage.mode(cm) <- "double"
samples <- colnames(cm)                 # sample ids = matrix columns
X <- t(cm)                              # vegan expects samples (rows) x taxa (cols)

# ---- optional grouping from samples.tsv ----
meta <- tryCatch(read.table(samples_path, header = TRUE, sep = "\t",
                            check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
grp <- NULL
if (!is.na(group_col) && !is.null(meta) && group_col %in% colnames(meta)) {
  g <- factor(meta[match(samples, meta[["sample"]]), group_col])
  if (length(g) == length(samples) && any(!is.na(g)) && nlevels(g) >= 2) grp <- g
}

# ---- alpha diversity: sample, richness, shannon, simpson ----
# chao1 is intentionally omitted: it is a count-based estimator (singleton/
# doubleton frequencies) and the input is fractional CoverM TPM, never counts.
if (!is.na(out_alpha)) {
  a <- data.frame(sample   = samples,
                  richness = as.numeric(specnumber(X)),
                  shannon  = as.numeric(diversity(X, "shannon")),
                  simpson  = as.numeric(diversity(X, "simpson")),
                  check.names = FALSE)
  write.table(a, out_alpha, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ---- beta diversity: pairwise Bray-Curtis, samples x samples ----
bray <- NULL
if (!is.na(out_beta) || !is.na(out_ord)) {
  bray <- vegdist(X, method = "bray")
}
if (!is.na(out_beta)) {
  M <- as.matrix(bray)
  write.table(data.frame(sample = rownames(M), M, check.names = FALSE),
              out_beta, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ---- ordination: PCoA (cmdscale) on Bray-Curtis ----
if (!is.na(out_ord)) {
  png(out_ord, width = 900, height = 700)
  k <- min(2, nrow(X) - 1L)
  ord <- if (k >= 1) tryCatch(cmdscale(bray, k = k, eig = TRUE), error = function(e) NULL) else NULL
  if (is.null(ord) || nrow(ord$points) < 1) {
    plot.new(); title(main = "Ordination unavailable (too few samples)")
  } else {
    pts <- ord$points
    if (ncol(pts) == 1) pts <- cbind(pts, 0)
    cols <- if (!is.null(grp)) as.numeric(grp) else rep(1, nrow(pts))
    plot(pts[, 1], pts[, 2], col = cols, pch = 19,
         xlab = "PCoA1", ylab = "PCoA2", main = "PCoA (Bray-Curtis)")
    text(pts[, 1], pts[, 2], labels = samples, pos = 3, cex = 0.7)
    if (!is.null(grp))
      legend("topright", legend = levels(grp), col = seq_along(levels(grp)),
             pch = 19, title = group_col)
  }
  dev.off()
}

# ---- differential abundance (only when an output path is given) ----
if (!is.na(out_diff)) {
  if (requireNamespace("ANCOMBC", quietly = TRUE))
    message("[ecology] ANCOMBC detected but base-R tests are used for a stable schema; see script to switch.")
  else
    message("[ecology] ANCOM-BC not installed; using base-R Wilcoxon/Kruskal + BH (install Bioconductor ANCOMBC to enable).")

  res <- NULL
  if (!is.null(grp)) {
    lv <- levels(grp); n <- nrow(cm)
    pv <- numeric(n); lfc <- numeric(n)
    for (i in seq_len(n)) {
      y <- cm[i, ]
      if (all(y == 0)) { pv[i] <- NA_real_; lfc[i] <- NA_real_; next }
      if (length(lv) == 2L) {
        tt <- tryCatch(wilcox.test(y[grp == lv[1]], y[grp == lv[2]]), error = function(e) NULL)
        pv[i] <- if (is.null(tt)) NA_real_ else tt$p.value
        m1 <- mean(y[grp == lv[1]]); m2 <- mean(y[grp == lv[2]])
        lfc[i] <- log2((m2 + 1) / (m1 + 1))
      } else {
        kt <- tryCatch(kruskal.test(y ~ grp), error = function(e) NULL)
        pv[i] <- if (is.null(kt)) NA_real_ else kt$p.value
        lfc[i] <- NA_real_
      }
    }
    res <- data.frame(votu = rownames(cm),
                      method = if (length(lv) == 2L) "wilcoxon" else "kruskal",
                      pvalue = pv, padj = p.adjust(pv, method = "BH"), log2fc = lfc,
                      check.names = FALSE)
    for (g in lv) res[[paste0("mean_", g)]] <- rowMeans(cm[, grp == g, drop = FALSE])
    res <- res[order(res$padj), ]
  } else {
    message("[ecology] differential skipped: no valid group column or <2 groups.")
    res <- data.frame(votu = character(), method = character(), pvalue = numeric(),
                      padj = numeric(), log2fc = numeric())
  }
  write.table(res, out_diff, sep = "\t", row.names = FALSE, quote = FALSE)
}
