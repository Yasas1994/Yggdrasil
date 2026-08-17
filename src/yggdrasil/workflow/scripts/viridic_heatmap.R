#!/usr/bin/env Rscript
# Lightweight VIRIDIC heatmap: render the intergenomic-similarity matrix as a
# *rasterized* PDF. Replaces VIRIDIC's native 03_viridic_heatmap.R step, whose
# per-cell ComplexHeatmap vector output blows up memory/file size for large
# vOTU sets (>= ~1000 genomes). Rows/cols are ordered by species cluster so the
# block structure of the similarity matrix is visible.

args <- commandArgs(trailingOnly = TRUE)
sim_path   <- args[1]  # similarity matrix TSV (sim_MA_genCol.csv)
clust_path <- args[2]  # clusters CSV (genome, species_cluster, genus_cluster)
out_pdf    <- args[3]

df <- read.delim(sim_path, row.names = 1, check.names = FALSE)
m <- as.matrix(df)
n <- nrow(m)

# Order genomes by species cluster (then genus) for block structure.
ord <- seq_len(n)
if (file.exists(clust_path) && file.size(clust_path) > 0) {
  cl <- read.delim(clust_path, stringsAsFactors = FALSE)
  cl <- cl[match(rownames(m), cl$genome), ]
  if (all(c("species_cluster", "genus_cluster") %in% colnames(cl))) {
    ord <- order(cl$species_cluster, cl$genus_cluster, na.last = TRUE)
  }
}
m <- m[ord, ord, drop = FALSE]

# Clip extreme size so the raster stays a sane resolution.
pdf(out_pdf, width = 10, height = 9, useDingbats = FALSE)
op <- par(mar = c(1, 1, 3, 9))
image(
  x = seq_len(n), y = seq_len(n), z = t(m[n:1, , drop = FALSE]),
  col = hcl.colors(256, "viridis"), zlim = c(0, 100),
  useRaster = TRUE, axes = FALSE,
  xlab = "", ylab = "",
  main = sprintf("VIRIDIC intergenomic similarity (%d vOTUs)", n)
)
# Colour legend.
par(xpd = TRUE)
legend_vals <- seq(0, 100, length.out = 256)
legend_x <- n * 1.04
legend_y <- seq(n * 0.15, n * 0.85, length.out = 256)
rect(legend_x, legend_y[-256], legend_x * 1.03, legend_y[-1],
     col = hcl.colors(256, "viridis"), border = NA)
text(legend_x * 1.06, seq(n * 0.15, n * 0.85, length.out = 5),
     labels = seq(0, 100, length.out = 5), cex = 0.8, adj = 0)
text(legend_x * 1.06, n * 0.92, "similarity (%)", cex = 0.8, adj = 0)
par(op)
dev.off()
