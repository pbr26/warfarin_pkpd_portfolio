# 01_data_cleaning.R
# Author: Pramod BR
# Date:   2026-07-24
# Load the warfarin data, inspect it, label PK vs PD, save analysis-ready sets.

library(nlmixr2data)
library(dplyr)

data(warfarin)

# Snapshot the raw data so the project is self-contained.
write.csv(warfarin, "data/raw/warfarin_raw.csv", row.names = FALSE)

# Quick look
str(warfarin)
summary(warfarin)

# The dataset is long format. Key columns:
#   id, time, amt, dv, dvid ("cp" = concentration, "pca" = effect),
#   evid (0 = observation, 1 = dose), wt, age, sex.

clean <- warfarin %>%
  mutate(
    dvid = as.character(dvid),
    sex  = as.factor(sex)
  ) %>%
  arrange(id, time)

# Dosing records (shared by both PK and PD models).
dose_records <- clean %>% filter(evid != 0)

# PK analysis set: dosing rows + concentration observations.
pk_data <- clean %>%
  filter(evid != 0 | dvid == "cp")

# PD analysis set: dosing rows + pca observations.
pd_data <- clean %>%
  filter(evid != 0 | dvid == "pca")

# Combined set for the joint PK/PD fit (both endpoints kept).
pkpd_data <- clean

# Basic data-quality checks
stopifnot(!any(is.na(pk_data$time)))
stopifnot(all(dose_records$amt > 0))

# Save processed sets
write.csv(pk_data,   "data/processed/pk_data.csv",   row.names = FALSE)
write.csv(pd_data,   "data/processed/pd_data.csv",   row.names = FALSE)
write.csv(pkpd_data, "data/processed/pkpd_data.csv", row.names = FALSE)
saveRDS(list(pk = pk_data, pd = pd_data, pkpd = pkpd_data),
        "data/processed/analysis_sets.rds")

message("01 done: ", nrow(pk_data), " PK rows, ", nrow(pd_data), " PD rows.")
