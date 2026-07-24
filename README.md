# Warfarin PK/PD Portfolio

[![View live report](https://img.shields.io/badge/View-Live%20Report-2C6E9B?style=for-the-badge)](https://019f94c6-630e-cde7-a624-e413f869e4f2.share.connect.posit.cloud/)

A complete, reproducible population pharmacokinetic/pharmacodynamic (PK/PD)
analysis of warfarin, built in R with **nlmixr2**.

**Live report:** <https://019f94c6-630e-cde7-a624-e413f869e4f2.share.connect.posit.cloud/>

**Goal:** demonstrate a full PK/PD workflow end to end — from raw data to a
fitted population model, diagnostics, and dosing simulations.

## Data

Uses the built-in `warfarin` dataset from the **`nlmixr2data`** package
(O'Reilly warfarin study: 32 subjects, single oral dose, plasma concentration
`cp` and prothrombin complex activity `pca`, with weight/age/sex covariates).
No download required — it installs with nlmixr2.

## Tooling

- **nlmixr2** — nonlinear mixed-effects model fitting (SAEM / FOCEi)
- **rxode2** — ODE simulation engine (dosing what-ifs)
- **ggplot2 / xpose** — plots and goodness-of-fit diagnostics
- **Quarto** — the narrative report (`analysis.qmd`)

## Pipeline

Run the scripts in order (each saves its outputs for the next):

| Script | Does |
|--------|------|
| `01_data_cleaning.R` | Load, inspect, label PK vs PD, save processed data |
| `02_EDA.R` | Exploratory plots: profiles, covariates, hysteresis |
| `03_PK_Model.R` | Fit the PK model (absorption, clearance, volume) |
| `04_PD_Model.R` | Fit the PD turnover / indirect-response model |
| `05_PKPD_Model.R` | Fit combined PK/PD with covariates |
| `06_Model_Diagnostics.R` | GOF plots, VPC, shrinkage |
| `07_Simulation.R` | Dosing-regimen simulations, time-to-target |
| `08_report_tables.R` | Export result tables for the published report |

## How to reproduce

```r
# 1. Install packages (once)
source("setup.R")

# 2. Run the pipeline
source("scripts/01_data_cleaning.R")
source("scripts/02_EDA.R")
source("scripts/03_PK_Model.R")
source("scripts/04_PD_Model.R")
source("scripts/05_PKPD_Model.R")
source("scripts/06_Model_Diagnostics.R")
source("scripts/07_Simulation.R")
source("scripts/08_report_tables.R")

# 3. Render the report
quarto::quarto_render("analysis.qmd")
```

To publish the report online (Posit Connect Cloud), see `PUBLISHING.md`.

## Folder layout

```
warfarin_pkpd_portfolio/
├── data/raw/          # raw data snapshot
├── data/processed/    # cleaned, analysis-ready data
├── scripts/           # 01–07 pipeline
├── models/            # saved fitted model objects (.rds)
├── figures/           # generated plots
├── reports/           # rendered report
├── setup.R            # installs/loads packages
├── analysis.qmd       # narrative report
└── LEARNING_NOTES.md  # beginner walkthrough of every step
```

## Notes

- The PD readout is **PCA** (prothrombin complex activity), not INR. INR is
  discussed in the report but not fabricated.
- The source data is **single-dose**, so multi-dose simulations in script 07
  are labelled as illustrative extrapolation.

**Start with `LEARNING_NOTES.md`** if you are new to PK/PD — it explains the
why behind every step.
