# 02_EDA.R
# Author: Pramod BR
# Date:   2026-07-24
# Exploratory plots to understand the data before modelling.

library(dplyr)
library(tidyr)
library(ggplot2)
source("scripts/00_plotting_theme.R")

sets <- readRDS("data/processed/analysis_sets.rds")
obs  <- sets$pkpd %>% filter(evid == 0)

# Mean +/- SD summary at each nominal time, for an overlay on the profiles.
summ <- function(df) {
  df %>% group_by(time) %>%
    summarise(m = mean(dv, na.rm = TRUE), sd = sd(dv, na.rm = TRUE),
              .groups = "drop")
}

# 1. Concentration profiles: faint per-subject lines + mean trend + SD ribbon.
cp <- obs %>% filter(dvid == "cp")
p_pk <- ggplot(cp, aes(time, dv)) +
  geom_line(aes(group = id), colour = pk_col, alpha = 0.18) +
  geom_point(colour = pk_col, alpha = 0.20, size = 1) +
  geom_ribbon(data = summ(cp), aes(time, ymin = pmax(m - sd, 0), ymax = m + sd),
              inherit.aes = FALSE, fill = pk_col, alpha = 0.15) +
  geom_line(data = summ(cp), aes(time, m), inherit.aes = FALSE,
            colour = pk_col, linewidth = 1.1) +
  labs(title = "Warfarin plasma concentration",
       subtitle = "Individual profiles with population mean ± SD",
       x = "Time (h)", y = "Concentration (mg/L)",
       caption = "nlmixr2data::warfarin — 32 subjects, single oral dose")
save_fig(p_pk, "figures/02_pk_profiles.png")

# 2. Effect (PCA) profiles.
pca <- obs %>% filter(dvid == "pca")
p_pd <- ggplot(pca, aes(time, dv)) +
  geom_line(aes(group = id), colour = pd_col, alpha = 0.18) +
  geom_point(colour = pd_col, alpha = 0.20, size = 1) +
  geom_ribbon(data = summ(pca), aes(time, ymin = m - sd, ymax = m + sd),
              inherit.aes = FALSE, fill = pd_col, alpha = 0.15) +
  geom_line(data = summ(pca), aes(time, m), inherit.aes = FALSE,
            colour = pd_col, linewidth = 1.1) +
  labs(title = "Prothrombin complex activity (effect)",
       subtitle = "Effect falls then slowly recovers — note the delay",
       x = "Time (h)", y = "PCA (%)",
       caption = "Lower PCA = greater anticoagulant effect")
save_fig(p_pd, "figures/02_pd_profiles.png")

# 3. Covariate distributions: histogram + density, clean facets.
covs <- sets$pkpd %>% distinct(id, wt, age, sex)
cov_long <- covs %>%
  pivot_longer(c(wt, age), names_to = "cov", values_to = "value") %>%
  mutate(cov = recode(cov, wt = "Weight (kg)", age = "Age (years)"))
p_cov <- ggplot(cov_long, aes(value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12,
                 fill = pk_col, alpha = 0.35, colour = "white") +
  geom_density(colour = ink, linewidth = 0.9) +
  facet_wrap(~ cov, scales = "free") +
  labs(title = "Covariate distributions",
       subtitle = "Demographics across the 32 subjects",
       x = NULL, y = "Density")
save_fig(p_cov, "figures/02_covariates.png", height = 4)

# 4. Hysteresis: effect vs concentration, coloured by time.
#    A loop (not a single line) => effect lags concentration => use an
#    indirect-response (turnover) PD model.
wide <- obs %>%
  filter(!is.na(dv), dvid %in% c("cp", "pca")) %>%
  group_by(id, time, dvid) %>%
  summarise(dv = mean(dv), .groups = "drop") %>%
  pivot_wider(names_from = dvid, values_from = dv) %>%
  filter(!is.na(cp), !is.na(pca)) %>%
  arrange(id, time)

p_hyst <- wide %>%
  filter(id == unique(id)[1]) %>%
  ggplot(aes(cp, pca, colour = time)) +
  geom_path(linewidth = 1, arrow = grid::arrow(length = unit(0.18, "cm"),
                                                type = "closed")) +
  geom_point(size = 2.6) +
  scale_colour_viridis_c(option = "magma", end = 0.9) +
  labs(title = "Concentration–effect hysteresis",
       subtitle = "Subject 1 — the loop is why we need a delayed PD model",
       x = "Concentration (mg/L)", y = "PCA (%)", colour = "Time (h)")
save_fig(p_hyst, "figures/02_hysteresis.png", width = 6.5, height = 5.2)

message("02 done: EDA figures written to figures/")
