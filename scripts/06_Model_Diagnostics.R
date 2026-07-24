# 06_Model_Diagnostics.R
# Author: Pramod BR
# Date:   2026-07-24
# Check whether the joint model fits well: GOF plots, residuals, VPC, shrinkage.

library(nlmixr2)
library(ggplot2)
library(dplyr)
source("scripts/00_plotting_theme.R")

pkpd_fit <- readRDS("models/pkpd_fit.rds")

# SAEM fits don't carry CWRES by default; add them (computed via FOCEi).
if (!"CWRES" %in% names(as.data.frame(pkpd_fit))) {
  pkpd_fit <- addCwres(pkpd_fit)
}

d <- as.data.frame(pkpd_fit) %>%
  mutate(CMT = recode(as.character(CMT),
                      cp = "PK: concentration", pca = "PD: effect (PCA)"))

# --- 1. Observed vs predicted (population and individual) --------------------
gof <- function(xvar, title) {
  ggplot(d, aes(.data[[xvar]], DV)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55",
                linetype = "dashed") +
    geom_point(aes(colour = CMT), alpha = 0.45, size = 1.6) +
    geom_smooth(method = "loess", se = FALSE, colour = ink,
                linewidth = 0.8, formula = y ~ x) +
    facet_wrap(~ CMT, scales = "free") +
    scale_colour_manual(values = c("PK: concentration" = pk_col,
                                   "PD: effect (PCA)" = pd_col),
                        guide = "none") +
    labs(title = title, x = xvar, y = "Observed (DV)")
}

p_pred  <- gof("PRED",  "Observed vs population prediction") +
  labs(subtitle = "Points should scatter evenly around the dashed identity line")
save_fig(p_pred, "figures/06_dv_vs_pred.png", width = 9, height = 4.4)

p_ipred <- gof("IPRED", "Observed vs individual prediction") +
  labs(subtitle = "Individual predictions should hug the identity line closely")
save_fig(p_ipred, "figures/06_dv_vs_ipred.png", width = 9, height = 4.4)

# --- 2. Conditional weighted residuals --------------------------------------
cwres_plot <- function(xvar, title) {
  ggplot(d, aes(.data[[xvar]], CWRES)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_hline(yintercept = c(-2, 2), colour = "grey80", linetype = "dotted") +
    geom_point(aes(colour = CMT), alpha = 0.45, size = 1.6) +
    geom_smooth(method = "loess", se = FALSE, colour = accent,
                linewidth = 1, formula = y ~ x) +
    scale_colour_manual(values = c("PK: concentration" = pk_col,
                                   "PD: effect (PCA)" = pd_col),
                        name = NULL) +
    labs(title = title, y = "CWRES")
}

p_cwres_t <- cwres_plot("TIME", "Conditional weighted residuals vs time") +
  labs(subtitle = "Should sit around 0 with no trend (yellow = loess)",
       x = "Time (h)")
save_fig(p_cwres_t, "figures/06_cwres_time.png", width = 8, height = 4.4)

p_cwres_p <- cwres_plot("PRED", "Conditional weighted residuals vs prediction") +
  labs(subtitle = "Should sit around 0 with no trend", x = "Population prediction")
save_fig(p_cwres_p, "figures/06_cwres_pred.png", width = 8, height = 4.4)

# --- 3. Eta shrinkage --------------------------------------------------------
print(pkpd_fit$shrink)

# --- 4. Visual predictive check ---------------------------------------------
vpc_res <- tryCatch({
  v <- vpcPlot(pkpd_fit, n = 300, show = list(obs_dv = TRUE)) +
    theme_pk() +
    labs(title = "Visual predictive check",
         subtitle = "Observed data should fall within the simulated intervals")
  save_fig(v, "figures/06_vpc.png", width = 8, height = 5)
  "VPC written"
}, error = function(e) paste("VPC skipped:", conditionMessage(e)))
message(vpc_res)

message("06 done: diagnostics written to figures/")
