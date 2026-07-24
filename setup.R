# setup.R — install and load everything the project needs.
# Author: Pramod BR
# Date:   2026-07-24
# Run this once before the pipeline: source("setup.R")

pkgs <- c(
  "nlmixr2",      # model fitting (pulls in rxode2, nlmixr2data, etc.)
  "nlmixr2data",  # the warfarin dataset
  "rxode2",       # ODE simulation
  "ggplot2",      # plots
  "dplyr",        # data wrangling
  "tidyr",        # data reshaping
  "xpose",        # goodness-of-fit diagnostics
  "xpose.nlmixr2" # xpose bridge for nlmixr2 fits
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(pkgs, library, character.only = TRUE))

# A place to keep helper settings used across scripts.
set.seed(1234)
theme_set(theme_bw())

message("Setup complete. Packages loaded.")
