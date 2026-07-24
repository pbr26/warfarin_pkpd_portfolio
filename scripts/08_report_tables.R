# 08_report_tables.R
# Author: Pramod BR
# Date:   2026-07-24
# Export fitted-model results as markdown tables so the published report can
# embed them WITHOUT re-fitting (fast, reliable render on Posit Connect Cloud).
# Run this locally after 03-07, then commit reports/results/ and figures/.

library(nlmixr2)

dir.create("reports/results", recursive = TRUE, showWarnings = FALSE)

# Turn a fit's parameter table into a tidy data frame.
tidy_params <- function(fit) {
  pf <- as.data.frame(fit$parFixedDf)
  pf <- cbind(Parameter = rownames(pf), pf)
  rownames(pf) <- NULL
  num <- vapply(pf, is.numeric, logical(1))
  pf[num] <- lapply(pf[num], function(x) signif(x, 3))
  pf
}

# Write a data frame as a GitHub-style markdown table with a caption.
write_md <- function(df, file, caption) {
  md <- knitr::kable(df, format = "pipe", caption = caption, row.names = FALSE)
  writeLines(md, file)
  message("wrote ", file)
}

pk_fit   <- readRDS("models/pk_fit.rds")
pd_fit   <- readRDS("models/pd_fit.rds")
pkpd_fit <- readRDS("models/pkpd_fit.rds")

write_md(tidy_params(pk_fit),   "reports/results/pk_params.md",
         "PK model — population parameter estimates")
write_md(tidy_params(pd_fit),   "reports/results/pd_params.md",
         "PD model — population parameter estimates")
write_md(tidy_params(pkpd_fit), "reports/results/pkpd_params.md",
         "Joint PK/PD model — population parameter estimates")

# Time-to-target from the simulation step.
sim <- readRDS("data/processed/simulation_results.rds")
target <- 25
tt <- do.call(rbind, lapply(split(sim, sim$scenario), function(s) {
  hit <- s$time[s$pca <= target]
  data.frame(Scenario = unique(s$scenario),
             `Time to PCA <= 25% (h)` =
               if (length(hit)) min(hit) else NA_real_,
             check.names = FALSE)
}))
write_md(tt, "reports/results/time_to_target.md",
         "Time to reach the illustrative therapeutic target (PCA <= 25%)")

message("08 done: report tables written to reports/results/")
