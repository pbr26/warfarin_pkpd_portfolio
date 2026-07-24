# 07_Simulation.R
# Author: Pramod BR
# Date:   2026-07-24
# Use the fitted model to simulate dosing scenarios with rxode2.
# NOTE: source data is single-dose, so multi-dose results are illustrative.

library(nlmixr2)
library(rxode2)
library(ggplot2)
library(tidyr)
library(dplyr)
source("scripts/00_plotting_theme.R")

pkpd_fit <- readRDS("models/pkpd_fit.rds")

# Build an rxode2 model matching the fitted structure, using the fitted
# typical-value (fixed-effect) parameter estimates.
sim_model <- rxode2({
  cl <- tcl * (WT / 70) ^ 0.75
  v  <- tv
  ka <- tka
  d/dt(depot)  <- -ka * depot
  d/dt(center) <-  ka * depot - (cl / v) * center
  cp <- center / v
  inh <- 1 - (imax * cp) / (ic50 + cp)
  pca(0) <- kin / kout
  d/dt(pca) <- kin * inh - kout * pca
})

# Pull typical values from the fit (back-transformed).
fe  <- fixef(pkpd_fit)
tka <- exp(fe[["tka"]]); tcl <- exp(fe[["tcl"]]); tv <- exp(fe[["tv"]])
kin <- exp(fe[["tkin"]]); kout <- exp(fe[["tkout"]]); ic50 <- exp(fe[["tic50"]])
imax <- plogis(fe[["timax"]])

params <- c(tka = tka, tcl = tcl, tv = tv,
            kin = kin, kout = kout, ic50 = ic50, imax = imax, WT = 70)

# Scenario A: single 100 mg dose (matches the study).
evA <- et(amt = 100, cmt = "depot") %>% et(seq(0, 240, by = 1))
simA <- rxSolve(sim_model, params, evA)

# Scenario B: 5 mg daily maintenance for 10 days (illustrative extrapolation).
evB <- et(amt = 5, cmt = "depot", ii = 24, addl = 9) %>% et(seq(0, 336, by = 1))
simB <- rxSolve(sim_model, params, evB)

target <- 25  # illustrative "therapeutic" PCA level

plot_df <- bind_rows(
  transform(as.data.frame(simA), scenario = "A: single 100 mg"),
  transform(as.data.frame(simB), scenario = "B: 5 mg daily x10 (illustrative)")
)

# First time PCA drops to/below target, per scenario (for annotation).
ttp <- plot_df %>%
  group_by(scenario) %>%
  filter(pca <= target) %>%
  summarise(t = ifelse(n() > 0, min(time), NA_real_), .groups = "drop")

# --- Plot 1: effect (PCA) over time, with therapeutic threshold -------------
p_sim <- ggplot(plot_df, aes(time, pca)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = target,
           fill = pd_col, alpha = 0.06) +
  geom_hline(yintercept = target, colour = pd_col, linetype = "dashed") +
  geom_line(colour = pd_col, linewidth = 1) +
  geom_vline(data = ttp, aes(xintercept = t), colour = accent,
             linetype = "dotted", na.rm = TRUE) +
  facet_wrap(~ scenario, scales = "free_x") +
  annotate("text", x = Inf, y = target, label = paste0("  target ", target, "%"),
           hjust = 1, vjust = -0.5, size = 3.2, colour = pd_col) +
  labs(title = "Simulated anticoagulant effect over time",
       subtitle = "Shaded band = therapeutic zone; dotted line = time to reach it",
       x = "Time (h)", y = "PCA (%)",
       caption = "Typical-value simulation. Scenario B extrapolates beyond single-dose data.")
save_fig(p_sim, "figures/07_dosing_scenarios.png", width = 9.5, height = 4.6)

# --- Plot 2: concentration and effect together (scenario B) -----------------
long_B <- as.data.frame(simB) %>%
  select(time, Concentration = cp, PCA = pca) %>%
  pivot_longer(-time, names_to = "series", values_to = "value")
p_pkpd <- ggplot(long_B, aes(time, value, colour = series)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(Concentration = pk_col, PCA = pd_col),
                      guide = "none") +
  labs(title = "Concentration drives effect (Scenario B)",
       subtitle = "Repeated dosing builds concentration; effect follows with a delay",
       x = "Time (h)", y = NULL)
save_fig(p_pkpd, "figures/07_conc_effect.png", width = 8, height = 5.4)

message("Time to reach PCA<=", target, "%:")
print(ttp)

saveRDS(plot_df, "data/processed/simulation_results.rds")
message("07 done: simulation figures written to figures/")
