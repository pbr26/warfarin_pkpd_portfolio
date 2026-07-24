# 05_PKPD_Model.R
# Author: Pramod BR
# Date:   2026-07-24
# Joint PK/PD model: fit concentration (cp) and effect (pca) simultaneously,
# and add a covariate (weight on clearance and volume).

library(nlmixr2)
library(dplyr)

sets      <- readRDS("data/processed/analysis_sets.rds")
pkpd_data <- sets$pkpd

pkpd_model <- function() {
  ini({
    # PK
    tka <- log(0.5)
    tcl <- log(0.13)
    tv  <- log(8)
    wt_cl <- 0.75      # allometric-style weight effect on CL
    wt_v  <- 1.0       # weight effect on V
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.v  ~ 0.1
    # PD
    tkin  <- log(1)
    tkout <- log(0.05)
    tic50 <- log(1)
    timax <- logit(0.9)
    eta.kout ~ 0.1
    # residual error
    add.err  <- 0.1
    prop.err <- 0.1
    pca.err  <- 5
  })
  model({
    # PK with weight normalised to 70 kg
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl) * (wt / 70) ^ wt_cl
    v  <- exp(tv  + eta.v)  * (wt / 70) ^ wt_v
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v

    # PD driven by model-predicted cp
    kin  <- exp(tkin)
    kout <- exp(tkout + eta.kout)
    ic50 <- exp(tic50)
    imax <- expit(timax)
    inh  <- 1 - (imax * cp) / (ic50 + cp)
    pca(0) <- kin / kout
    d/dt(pca) <- kin * inh - kout * pca

    # two endpoints
    cp  ~ add(add.err) + prop(prop.err)
    pca ~ add(pca.err)
  })
}

pkpd_fit <- nlmixr2(pkpd_model, pkpd_data, est = "saem",
                    control = saemControl(print = 0))

print(pkpd_fit)
saveRDS(pkpd_fit, "models/pkpd_fit.rds")
message("05 done: joint PK/PD model saved to models/pkpd_fit.rds")
