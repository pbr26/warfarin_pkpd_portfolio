# 04_PD_Model.R
# Author: Pramod BR
# Date:   2026-07-24
# Fit the PD model: how warfarin concentration drives the effect (PCA).
# Warfarin inhibits synthesis of clotting factors -> an indirect-response
# (turnover) model with inhibition of the production rate (kin).

library(nlmixr2)
library(dplyr)

sets <- readRDS("data/processed/analysis_sets.rds")

# In a SEQUENTIAL PD fit the model needs the drug concentration (cp) at every
# record. The raw data is long format (cp and pca on separate rows), so we
# build a cp column by linearly interpolating each subject's observed
# concentrations onto all of their PD record times.
conc <- sets$pkpd %>%
  filter(dvid == "cp", evid == 0) %>%
  select(id, time, cp = dv)

pd_data <- sets$pd %>%
  group_by(id) %>%
  group_modify(function(df, key) {
    ci <- conc[conc$id == key$id, ]
    df$cp <- if (nrow(ci) >= 2) {
      approx(ci$time, ci$cp, xout = df$time, rule = 2)$y
    } else {
      0
    }
    df
  }) %>%
  ungroup()
pd_data$cp[is.na(pd_data$cp)] <- 0

# Turnover model with Imax inhibition on production.
#   baseline PCA = kin/kout
#   drug reduces kin via  (1 - Imax*cp/(IC50+cp))
pd_turnover <- function() {
  ini({
    tkin  <- log(1)      # zero-order production rate
    tkout <- log(0.05)   # first-order loss rate
    tic50 <- log(1)      # concentration giving half-maximal inhibition
    timax <- logit(0.9)  # maximum inhibition (0-1)
    eta.kout ~ 0.1
    pca.err <- 5
  })
  model({
    kin  <- exp(tkin)
    kout <- exp(tkout + eta.kout)
    ic50 <- exp(tic50)
    imax <- expit(timax)
    # cp is supplied per record from the observed/PK-predicted concentration.
    inh  <- 1 - (imax * cp) / (ic50 + cp)
    pca(0) <- kin / kout
    d/dt(pca) <- kin * inh - kout * pca
    pca ~ add(pca.err)
  })
}

# Sequential PD fit: uses the interpolated cp column built above.
# (In 05 we fit PK and PD jointly; here we isolate the PD structure.)
pd_fit <- nlmixr2(pd_turnover, pd_data, est = "saem",
                  control = saemControl(print = 0))

print(pd_fit)
saveRDS(pd_fit, "models/pd_fit.rds")
message("04 done: PD model saved to models/pd_fit.rds")
