library(glmnet)
library(pROC)
library(PRROC)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)  

########################
# Load clock models #
########################

clock = readRDS("./results/clock_v3.3_update.rds")
cell_type_df = read.csv("./../preprocessing/cell_type_mapping_update.csv")
cell_type_df$Somamers_list <- sapply(strsplit(cell_type_df$Somamers, ",\\s*"), function(x) x)

plmd = read.csv('./data/GNPC_processed_proteomics_clinical_final.csv')
pl_meta = plmd[,1:58]
pl_prot = plmd[,59:ncol(plmd)]
prot_std_new = pl_prot

# 2 is Female
sex_idx <- which(pl_meta$sex %in% c(1,2))
sex_F = as.numeric(pl_meta$sex[sex_idx] == 2)

pl_meta_sex_idx = pl_meta[sex_idx,]
#write.csv(pl_meta_sex_idx, "./results/gnpc_meta_sex_idx_update.csv", row.names = FALSE)

#######################
# Predict in GNPC #
#######################

prot_sex = cbind(sex_F, prot_std_new[sex_idx,])

predict_bootstrap_models <- function(models, x_test) {
  predictions <- matrix(NA, nrow = nrow(as.matrix(x_test)), ncol = length(models))
  for (i in 1:length(models)) {
    predictions[, i] <- predict(models[[i]], newx = as.matrix(x_test))
  }
  return(rowMeans(predictions))
}

sex_F = as.numeric(pl_meta$sex[sex_idx] == 2)
prot_std_new_fil = prot_std_new[sex_idx,]
chrono_age = pl_meta$age_at_visit[sex_idx]

cell_type_low_performance = c()
results <- list()
for (cell_type in names(clock)) {
 
    somamer_ids = cell_type_df[cell_type_df$Original.Cell.types == cell_type,]$Somamers_list[[1]]
    n_marker_total = length(somamer_ids)
    
    # Extract relevant protein data and combine with sex
    x_test <- cbind(
      Sex = sex_F,
      prot_std_new_fil[, somamer_ids]
    )
    
    # Get predictions using the bootstrap models for this cell type
    predictions <- predict_bootstrap_models(
      clock[[cell_type]]$models, 
      as.matrix(x_test)
    )
    
    # Store results for this cell type
    results[[cell_type]] <- list(
      predictions = predictions,
      evaluation_gnpc = cor(predictions, pl_meta$age[sex_idx]),
      n_marker = n_marker_total
    )
    
    cat("\nResults for", cell_type, ":\n")
    cat("Correlation:", round(results[[cell_type]]$evaluation_gnpc, 3), "\n")
    
    if (results[[cell_type]]$evaluation_gnpc < 0.1){
      cell_type_low_performance = c(cell_type_low_performance, cell_type)
      cat("\nLow Correlation", cell_type, ":\n")
      print(results[[cell_type]]$evaluation_gnpc)
    }
    
}

saveRDS(results, file = "results/Rev_UPDATE_clock_prediction_v3.3_update.rds")
results = readRDS("results/Rev_UPDATE_clock_prediction_v3.3_update.rds")

pred_results <- data.frame(
  sample_id = pl_meta$sample_id[sex_idx],
  Chronological_Age = pl_meta$age_at_visit[sex_idx]
)
for (cell_type in names(results)) {
  predictions <- results[[cell_type]]$predictions
  col_name <- paste0(make.names(cell_type), "_Age")
  pred_results[[col_name]] <- predictions
}
write.csv(pred_results, "./results/Rev_UPDATE_cell_type_predictions_update.csv", row.names = FALSE)

cell_type_performance <- list()
cell_type_markers_n = list()
for (cell_type in names(results)) {
  cell_type_performance[[cell_type]] <- results[[cell_type]]$evaluation_gnpc
  cell_type_markers_n[[cell_type]] <- results[[cell_type]]$n_marker
}
cell_type_performance_df <- data.frame(
  cell_type = names(cell_type_performance),
  evaluation_gnpc = unlist(cell_type_performance),
  n_marker = unlist(cell_type_markers_n),
  stringsAsFactors = FALSE
)
write.csv(cell_type_performance_df, "./results/Rev_UPDATE_cell_type_performance_gnpc_update.csv", row.names = FALSE)

