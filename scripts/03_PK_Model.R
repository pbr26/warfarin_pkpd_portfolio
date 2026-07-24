# 03_PK_Model.R
# Author: Pramod BR
# Date:   2026-07-24
# Fit the population PK model: how the body absorbs and clears warfarin.

library(nlmixr2)
library(dplyr)

sets    <- readRDS("data/processed/analysis_sets.rds")
pk_data <- sets$pk %>% filter(dvid == "cp" | evid != 0)

# One-compartment, first-order oral absorption.
pk_1cmt <- function() {
  ini({
    tka  <- log(0.5)    # absorption rate (ka)
    tcl  <- log(0.13)   # clearance (CL)
    tv   <- log(8)      # volume (V)
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.v  ~ 0.1
    add.err  <- 0.1
    prop.err <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl / v) * center
    cp <- center / v
    cp ~ add(add.err) + prop(prop.err)
  })
}

# Two-compartment alternative, to justify the structural choice.
pk_2cmt <- function() {
  ini({
    tka <- log(0.5); tcl <- log(0.13); tv <- log(8)
    tq  <- log(0.5); tv2 <- log(4)
    eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1
    add.err <- 0.1; prop.err <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    q  <- exp(tq)
    v2 <- exp(tv2)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - (cl/v)*center - (q/v)*center + (q/v2)*periph
    d/dt(periph) <-  (q/v)*center - (q/v2)*periph
    cp <- center / v
    cp ~ add(add.err) + prop(prop.err)
  })
}

fit_1 <- nlmixr2(pk_1cmt, pk_data, est = "saem",
                 control = saemControl(print = 0))
fit_2 <- nlmixr2(pk_2cmt, pk_data, est = "saem",
                 control = saemControl(print = 0))

# Compare by objective function / AIC (lower = better, penalising complexity).
cat("1-cmt AIC:", AIC(fit_1), "  2-cmt AIC:", AIC(fit_2), "\n")
pk_fit <- if (AIC(fit_1) <= AIC(fit_2)) fit_1 else fit_2

print(pk_fit)
saveRDS(pk_fit, "models/pk_fit.rds")
message("03 done: PK model saved to models/pk_fit.rds")
