library(glmnet)
library(pROC)
library(PRROC)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)  

gnpc_hc = read.csv('./data/GNPC_HC_clinical.csv')
gnpc_all = read.csv('./data/GNPC_clinical_final.csv')

pred_results <- read.csv("./results/Rev_UPDATE_cell_type_predictions_update.csv")

# Calculate age gaps from using LOWESS fit on healthy controls
estimate_age_gaps_with_fit <- function(chronological_age, biological_age, lowess_fit) {
  predicted <- approx(x = lowess_fit$x, y = lowess_fit$y, xout = chronological_age)$y
  gaps <- biological_age - predicted
  list(gaps = gaps, predicted = predicted, fit = lowess_fit)
}

all_contributors <- unique(gnpc_all$contributor_code)

kadrc_hc <- gnpc_hc[gnpc_hc$contributor_code == "F", ]
kadrc_hc_sample_ids <- kadrc_hc$sample_id
kadrc_hc_indices <- which(pred_results$sample_id %in% kadrc_hc_sample_ids)

# Get cell type column names
cell_type_cols <- setdiff(colnames(pred_results), c("sample_id", "Chronological_Age"))

# Initialize list to store gaps for all contributors
all_contributor_gaps <- list()

# Process each cohort separately: fit LOWESS on cohort's HC
for (contributor in all_contributors) {
  cat("Processing contributor:", contributor, "\n")
  
  # Filter samples for this contributor
  contributor_data <- gnpc_all[gnpc_all$contributor_code == contributor, ]
  contributor_sample_ids <- contributor_data$sample_id
  contributor_indices <- which(pred_results$sample_id %in% contributor_sample_ids)
  
  # Get healthy controls for this contributor
  contributor_hc <- gnpc_hc[gnpc_hc$contributor_code == contributor, ]
  contributor_hc_sample_ids <- contributor_hc$sample_id
  contributor_hc_indices <- which(pred_results$sample_id %in% contributor_hc_sample_ids)
  
  # Check if we have enough healthy controls
  has_enough_hc <- length(contributor_hc_indices) >= 5
  
  if (has_enough_hc) {
    cat("  Using contributor's own healthy controls (n=", length(contributor_hc_indices), ")\n", sep="")
  } else {
    cat("  Insufficient healthy controls, falling back to KADRC\n")
  }
  
  # Initialize list to store gaps for this contributor
  gnpc_cell_gaps_contributor <- list()
  
  # Process each cell type
  for (cell_col in cell_type_cols) {
    # Extract the original cell type name from the column name
    cell_type <- gsub("_Age$", "", cell_col)
    cell_type <- gsub("\\.", " ", cell_type)
    
    use_kadrc <- FALSE
    
    # Try to fit LOWESS on cohort's own healthy controls
    if (has_enough_hc) {
      hc_chrono_age <- pred_results$Chronological_Age[contributor_hc_indices]
      hc_bio_age <- pred_results[[cell_col]][contributor_hc_indices]
      
      # Fit LOWESS on the contributor's healthy controls
      hc_fit <- try(lowess(x = hc_chrono_age, y = hc_bio_age, f = 2/3), silent = TRUE)
      
      if (inherits(hc_fit, "try-error")) {
        cat("  Warning: LOWESS fit failed for", cell_type, "using contributor's HC, falling back to KADRC\n")
        use_kadrc <- TRUE
      }
    } else {
      use_kadrc <- TRUE
    }
    
    # Fallback to KADRC healthy controls if needed
    if (use_kadrc) {
      kadrc_hc_chrono_age <- pred_results$Chronological_Age[kadrc_hc_indices]
      kadrc_hc_bio_age <- pred_results[[cell_col]][kadrc_hc_indices]
      
      if (length(kadrc_hc_chrono_age) < 5) {
        cat("  Warning: Too few KADRC HC samples for cell type", cell_type, ", skipping\n")
        next
      }
      
      hc_fit <- lowess(x = kadrc_hc_chrono_age, y = kadrc_hc_bio_age, f = 2/3)
    }
    
    # Calculate age gaps for all samples in this cohort
    contributor_chrono_age <- pred_results$Chronological_Age[contributor_indices]
    contributor_bio_age <- pred_results[[cell_col]][contributor_indices]
    
    # Calculate gaps using the selected healthy controls-based LOWESS fit
    res <- estimate_age_gaps_with_fit(
      chronological_age = contributor_chrono_age, 
      biological_age = contributor_bio_age, 
      lowess_fit = hc_fit
    )
    
    # Z-score normalize gaps within cohort
    res$gaps_z <- scale(res$gaps)
    res$sample_id <- pred_results$sample_id[contributor_indices]
    
    gnpc_cell_gaps_contributor[[cell_type]] <- res}

  all_contributor_gaps[[contributor]] <- gnpc_cell_gaps_contributor
}

# Combine gaps from all cohorts into final dataframes
gap_df_final <- data.frame(sample_id = pred_results$sample_id)
gap_df_final$Chronological_Age <- pred_results$Chronological_Age

# For each cell type, create a column in the final dataframe
all_cell_types <- unique(unlist(lapply(all_contributor_gaps, names)))

# Initialize columns for each cell type with NAs
for (cell_type in all_cell_types) {
  gap_df_final[[cell_type]] <- NA
}

gap_df_final_z_within = gap_df_final
gap_df_final_z_across = gap_df_final

for (contributor in names(all_contributor_gaps)) {
  print(contributor)
  for (cell_type in names(all_contributor_gaps[[contributor]])) {
    contributor_gaps <- all_contributor_gaps[[contributor]][[cell_type]]
    
    for (i in 1:length(contributor_gaps$sample_id)) {
      sample_id <- contributor_gaps$sample_id[i]
      idx <- which(gap_df_final$sample_id == sample_id)
      
      if (length(idx) > 0) {
        gap_df_final_z_within[idx, cell_type] <- contributor_gaps$gaps_z[i]
        gap_df_final_z_across[idx, cell_type] <- contributor_gaps$gaps[i]
      }
    }
  }
}

for (cell_type in all_cell_types) {
  if (cell_type %in% c("sample_id", "Chronological_Age")) {
    next
  }
  gap_df_final_z_across[[cell_type]] <- scale(gap_df_final_z_across[[cell_type]])
}

write.csv(gap_df_final_z_across, file = "./results/Rev_UPDATE_gnpc_cell_type_gaps_z_v3.3_by_COHORT_HC_across_cohorts_z.csv", row.names = FALSE)
write.csv(gap_df_final_z_within, file = "./results/Rev_UPDATE_gnpc_cell_type_gaps_z_v3.3_by_COHORT_HC_within_cohorts_z.csv", row.names = FALSE)


