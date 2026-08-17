library(glmmTMB)
library(dplyr)
library(broom.mixed)

df <- read.csv("./prop_temperate.csv") %>%
  rename(mouse_id = `Mouse.Number.Control`,
         site = `Feces.Cecum.Control.CRC`,
         tumor_status = `Mock.Tumor`) %>%
  mutate(
    mouse_id = factor(mouse_id),
    site = factor(site, levels = c("Ileum", "Cecum")),
    tumor_status = factor(tumor_status, levels = c("Mock", "Tumor"))
  )

# Smithson & Verkuilen boundary adjustment (standard for beta regression)
n <- nrow(df)
df <- df %>% mutate(p_adj = (p_temperate * (n - 1) + 0.5) / n)

# --- full dataset ---
model_full <- glmmTMB(
  p_adj ~ tumor_status * site + (1 | mouse_id),
  data = df,
  family = beta_family()
)

# --- sensitivity: drop Kenneth011 ---
df_sens <- df %>% filter(sample != "Kenneth011")

model_sens <- glmmTMB(
  p_adj ~ tumor_status * site + (1 | mouse_id),
  data = df_sens,
  family = beta_family()
)

# side-by-side comparison of the interaction term specifically
comparison <- bind_rows(
  tidy(model_full, effects = "fixed", conf.int = TRUE) %>% mutate(model = "full"),
  tidy(model_sens, effects = "fixed", conf.int = TRUE) %>% mutate(model = "excl_Kenneth011")
) %>%
  filter(grepl(":", term)) %>%
  select(model, term, estimate, conf.low, conf.high, p.value)

print(comparison)


model_cecum <- glmmTMB(p_adj ~ tumor_status, data = df %>% filter(site == "Cecum"), family = beta_family())
model_ileum <- glmmTMB(p_adj ~ tumor_status, data = df %>% filter(site == "Ileum"), family = beta_family())

summary(model_cecum)

# or, consistent with the earlier table format
tidy(model_cecum, effects = "fixed", conf.int = TRUE)

tidy(model_ileum, effects = "fixed", conf.int = TRUE)

model_additive <- glmmTMB(
  p_adj ~ tumor_status + site + (1 | mouse_id),
  data = df,
  family = beta_family()
)

summary(model_additive)

tidy(model_additive, effects = "fixed", conf.int = TRUE)

model_pooled_naive <- glmmTMB(
  p_adj ~ tumor_status + (1 | mouse_id),
  data = df,
  family = beta_family()
)

tidy(model_pooled_naive, effects = "fixed", conf.int = TRUE)
