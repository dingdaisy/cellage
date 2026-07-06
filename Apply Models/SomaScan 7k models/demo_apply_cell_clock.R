################################################################################
# CELL AGING CLOCK APPLICATION (EXTERNAL COHORTS)
# Purpose: Apply the trained cell aging clocks to a new proteomic dataset
################################################################################

### 1. Load Pre-trained Clocks #################################################
# Load the trained cell aging clocks coefficients table (intercepts and weights)
# Extract this file from the supplementary files or repository.
coef_df <- read.csv("olink_clock_coefficients_min.csv", row.names = 1, check.names = FALSE)

### 2. Load Your Dataset ########################################################
# Load your full cohort's Olink proteomics dataset (e.g., WHI data)
# IMPORTANT: NPX values must be Z-scored.
df_cohort <- read.csv("path_to_your_proteomics_dataset.csv")

# Ensure required covariates exist in your dataframe
# df_cohort MUST contain:
# 1. `eid`: Participant ID, or other identifier
# 2. `actual_age`: Chronological age in years
# 3. `Sex`: Binary indicator (1 = Female, 0 = Male)
# 4. `Diagnosis`: E.g., "Healthy" or other disease status 

# Create output dataframe to hold biological ages
all_preds <- data.frame(eid = df_cohort$eid, actual_age = df_cohort$actual_age, Diagnosis = df_cohort$Diagnosis)

### 3. Predict Raw Biological Ages ##############################################
cell_types <- rownames(coef_df)

for (ct in cell_types) {
  # Extract the coefficients specifically for this cell type model
  cfs_row <- coef_df[ct, ]
  cfs_row <- cfs_row[, !is.na(cfs_row), drop = FALSE]
  
  # Separate intercept from weights
  interc <- as.numeric(cfs_row[["intercept"]])
  cfs_row[["intercept"]] <- NULL
  
  features <- names(cfs_row)
  
  # Pad missing proteins with 0 (Assuming NPX data is Z-scored, 0 = population mean)
  # This allows the clock to run even if a specific marker protein is missing in your specific platform
  for (f in features) {
    if (!(f %in% names(df_cohort))) df_cohort[[f]] <- 0
  }
  
  # Extract features matrix and calculate biological age: (NPX_Matrix * Weights) + Intercept
  all_mat <- as.matrix(df_cohort[, features, drop = FALSE])
  pred_col <- paste0(ct, "_predicted_age")
  all_preds[[pred_col]] <- (all_mat %*% as.numeric(cfs_row)) + interc
}

### 4. Baseline Fit and Age Gap Calculation ######################################
# Helper function to interpolate expected ages from the LOWESS curve
estimate_age_gaps_with_fit <- function(chronological_age, biological_age, lowess_fit) {
  predicted <- approx(x = lowess_fit$x, y = lowess_fit$y, xout = chronological_age)$y
  gaps <- biological_age - predicted
  list(gaps = gaps, predicted = predicted, fit = lowess_fit)
}

# Output dataframe for gaps
all_gaps_df <- data.frame(eid = all_preds$eid, actual_age = all_preds$actual_age)

for (ct in cell_types) {
  pred_col <- paste0(ct, "_predicted_age")
  
  # Isolate healthy controls to establish an accurate healthy aging baseline
  # (if you don't have enough healthy controls, use the whole cohort `all_preds`)
  hc_preds <- all_preds[all_preds$Diagnosis == "Healthy", ]
  
  if(nrow(hc_preds) < 5) {
     cat(sprintf("Warning: Insufficient healthy controls for %s. Fitting LOWESS on entire cohort.\n", ct))
     hc_preds <- all_preds
  }
  
  # Fit LOWESS strictly on the base group (Healthy Controls)
  all_fit <- lowess(x = hc_preds$actual_age, y = hc_preds[[pred_col]], f = 2/3)
  
  # 1. Calculate raw age gaps based on the LOWESS trajectory for all individuals
  all_gaps_res <- estimate_age_gaps_with_fit(all_preds$actual_age, all_preds[[pred_col]], all_fit)
  
  # 2. Normalize using the cohort's mean and standard deviation
  cohort_mean <- mean(all_gaps_res$gaps, na.rm = TRUE)
  cohort_sd <- sd(all_gaps_res$gaps, na.rm = TRUE)
  
  # Produce standardized Z-gaps (comparable across cell types)
  all_gaps_z <- (all_gaps_res$gaps - cohort_mean) / cohort_sd
  
  all_gaps_df[[paste0(ct, "_age_gap")]] <- all_gaps_res$gaps
  all_gaps_df[[paste0(ct, "_age_gap_z")]] <- all_gaps_z
}

# Save final gap outputs
write.csv(all_gaps_df, "external_cohort_cell_age_gaps.csv", row.names = FALSE)
cat("Successfully computed Cell Aging Clocks and Age Gaps for external cohort.\n")
