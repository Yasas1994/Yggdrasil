library(maaslin3)
library(dplyr)

feature_table <- read.delim("./feature_table.tsv", row.names = "sample")
metadata <- read.delim("./metadata.tsv", row.names = "sample")

fit <- maaslin3(
  input_data = feature_table,
  input_metadata = metadata,
  output = "./maaslin3_output",
  formula = "~ tumor_status + site + (1 | mouse_id)",
  normalization = "NONE",
  transform = "LOG",
  median_comparison_abundance = FALSE,
  augment = TRUE,
  standardize = TRUE,
  reference = c("tumor_status,Mock", "site,Ileum"),
  small_random_effects = TRUE
)

metadata_sens <- metadata[rownames(metadata) != "Kenneth011", , drop = FALSE]
feature_table_sens <- feature_table[rownames(feature_table) != "Kenneth011", , drop = FALSE]

# sanity check before running - confirm Kenneth011 is actually gone from both
stopifnot(!"Kenneth011" %in% rownames(metadata_sens))
stopifnot(!"Kenneth011" %in% rownames(feature_table_sens))
nrow(metadata_sens)      # should be 23, not 24
nrow(feature_table_sens) # should be 23, not 24

fit_sens <- maaslin3(
  input_data = feature_table_sens,
  input_metadata = metadata_sens,
  output = "./maaslin3_output_sens_excl_Kenneth011",
  formula = "~ tumor_status + site + (1 | mouse_id)",
  normalization = "NONE",
  transform = "LOG",
  median_comparison_abundance = FALSE,
  augment = TRUE,
  standardize = TRUE,
  reference = c("tumor_status,Mock", "site,Ileum"),
  small_random_effects = TRUE
)
