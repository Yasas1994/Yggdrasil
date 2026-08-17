library(dplyr)

# --- 0. Verify both output directories actually contain what we expect ---
stopifnot(dir.exists("./maaslin3_output"))
stopifnot(dir.exists("./maaslin3_output_sens_excl_Kenneth011"))

print(list.files("./maaslin3_output", recursive = TRUE))
print(list.files("./maaslin3_output_sens_excl_Kenneth011", recursive = TRUE))

# adjust these two paths if the filenames above differ from "all_results.tsv"
full_results_path <- "./maaslin3_output/all_results.tsv"
sens_results_path  <- "./maaslin3_output_sens_excl_Kenneth011/all_results.tsv"

stopifnot(file.exists(full_results_path))
stopifnot(file.exists(sens_results_path))

# --- 1. Load full-dataset results ---
full_results <- read.delim(full_results_path)

# --- 2. Flag and split hits (using .data pronoun to avoid name collision with the metadata data.frame) ---
tumor_hits <- full_results %>%
  filter(.data$metadata == "tumor_status", .data$model == "abundance") %>%
  arrange(qval_individual) %>%
  mutate(
    prevalence_frac = N_not_zero / N,
    flag_large_coef = abs(coef) > 5,
    flag_low_prevalence = prevalence_frac < 0.7,
    needs_scrutiny = flag_large_coef | flag_low_prevalence
  )

clean_hits <- tumor_hits %>%
  filter(qval_individual < 0.05, !needs_scrutiny)

flagged_hits <- tumor_hits %>%
  filter(qval_individual < 0.05, needs_scrutiny)

print(clean_hits %>% select(feature, coef, stderr, qval_individual, N_not_zero))
print(flagged_hits %>% select(feature, coef, stderr, qval_individual, N_not_zero))

# --- 3. Directional asymmetry check ---
direction_summary <- tumor_hits %>%
  filter(qval_individual < 0.05, !needs_scrutiny) %>%
  summarise(n_up = sum(coef > 0), n_down = sum(coef < 0))

print(direction_summary)
print(binom.test(direction_summary$n_up, direction_summary$n_up + direction_summary$n_down, p = 0.5))

# --- 4. Load sensitivity (excl. Kenneth011) results and join ---
sens_results <- read.delim(sens_results_path) %>%
  filter(.data$metadata == "tumor_status", .data$model == "abundance") %>%
  select(feature, coef_sens = coef, qval_sens = qval_individual)

robustness_check <- clean_hits %>%
  left_join(sens_results, by = "feature") %>%
  mutate(
    survived = qval_sens < 0.05,
    coef_direction_stable = sign(coef) == sign(coef_sens)
  )

# how many clean_hits features had NO match in the sensitivity results at all
# (would indicate a join problem, e.g. feature naming mismatch between runs)
n_unmatched <- sum(is.na(robustness_check$coef_sens))
print(paste("Unmatched features (NA after join):", n_unmatched))

print(robustness_check %>% filter(!survived | !coef_direction_stable) %>%
        select(feature, coef, qval_individual, coef_sens, qval_sens))

robustness_check %>% summarise(
  n_total = n(),
  n_survived = sum(survived, na.rm = TRUE),
  n_direction_flipped = sum(!coef_direction_stable, na.rm = TRUE)
)

# Of 1,834 vOTUs tested, 63 showed significant differential abundance by tumor status 
# after FDR correction (q<0.05), independent of gut site. Of these, 56 were elevated 
# and 7 depleted in tumor-bearing mice (binomial test, p≈[your value]). Results were
# robust to exclusion of the highest-depth sample: 55/63 (87%) remained significant
# with unchanged effect direction; the remaining 8 showed stable coefficient estimates 
# but q-values shifting marginally above threshold, consistent with reduced power at
# n=23 rather than loss of signal. Three additional vOTUs met the significance 
# threshold but were excluded from primary reporting due to low prevalence or 
# coefficient instability suggestive of estimation artifacts (see Methods/Supplementary).

final_votu_list <- clean_hits %>%
  left_join(sens_results, by = "feature") %>%
  mutate(
    survived_sensitivity = qval_sens < 0.05,
    coef_direction_stable = sign(coef) == sign(coef_sens),
    status = case_when(
      survived_sensitivity & coef_direction_stable ~ "robust",
      !survived_sensitivity & coef_direction_stable ~ "marginal_borderline",
      TRUE ~ "unstable"
    ),
    direction = if_else(coef > 0, "up_in_tumor", "down_in_tumor")
  ) %>%
  select(
    feature, coef, qval_individual, N_not_zero,
    coef_sens, qval_sens, status, direction
  ) %>%
  arrange(status, qval_individual)

print(final_votu_list, n = 63)

# quick tally to confirm counts match what we already established
final_votu_list %>% count(status, direction)

final_votu_list
write.csv(final_votu_list, "./final_votu_list_tumor_associated.csv", row.names = FALSE)

# plain list of representative IDs, useful for pulling sequences / running host prediction
writeLines(final_votu_list$feature, "./final_votu_ids.txt")

# if you want just the robust set for the primary analysis, and the marginal set separately
robust_ids <- final_votu_list %>% filter(status == "robust") %>% pull(feature)
marginal_ids <- final_votu_list %>% filter(status == "marginal_borderline") %>% pull(feature)

writeLines(robust_ids, "./robust_votu_ids.txt")
writeLines(marginal_ids, "./marginal_votu_ids.txt")
