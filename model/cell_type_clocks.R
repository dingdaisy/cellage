library(glmnet)
library(pROC)
library(PRROC)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#######################
# Load proteomics data
#######################

kadrc_df = read.csv("./data/KADRC_SomaScan_Plasma_7K_ANML_clean_v2.csv", row.names = 1)

df_prot_kadrc = kadrc_df[, 21:ncol(kadrc_df)]
df_meta_kadrc = kadrc_df[, 1:20]

df_prot_kadrc_log = log10(df_prot_kadrc)
df_prot_kadrc_std = scale(df_prot_kadrc_log)

#######################
# Helper funcitons for training cellular aging clocks
#######################

train_bootstrap_models <- function(x_train, y_train, n_bootstrap = 10, 
                                   alpha_level) {
  bootstrap_models <- list()
  n_samples <- length(y_train)
  
  for (i in 1:n_bootstrap) {
    boot_indices <- sample(1:n_samples, n_samples, replace = TRUE)
    x_boot <- x_train[boot_indices, ]
    y_boot <- y_train[boot_indices]
    
    # Perform 10-fold cross-validation to select optimal lambda (regularization parameter)
    cv_fit <- cv.glmnet(as.matrix(x_boot), y_boot, alpha = alpha_level, nfolds = 10)
    selected_lambda = cv_fit$lambda.min  # Lambda that minimizes cross-validation error
    
    # Fit final model with selected lambda
    bootstrap_models[[i]] <- glmnet(as.matrix(x_boot), y_boot, 
                                    alpha = alpha_level,
                                    lambda=selected_lambda)
  }
  return(bootstrap_models)
}

predict_bootstrap_models <- function(models, x_test) {
  predictions <- matrix(NA, nrow = nrow(as.matrix(x_test)), ncol = length(models))
  for (i in 1:length(models)) {
    predictions[, i] <- predict(models[[i]], newx = as.matrix(x_test))
  }
  return(rowMeans(predictions))
}

evaluate_predictions <- function(predictions, true_values) {
  correlation <- cor(predictions, true_values)
  return(correlation)
}


#######################
# Fit cellular aging models
#######################

# Load cell type to protein mapping (created from Human Protein Atlas)
cell_type_df = read.csv("./../preprocessing/cell_type_mapping_update.csv")
cell_type_df$Somamers_list <- sapply(strsplit(cell_type_df$Somamers, ",\\s*"), function(x) x)

# Prepare data with sex adjustment (sex as a covariate)
prot_with_sex = cbind(Sex = df_meta_kadrc$Sex == "F", df_prot_kadrc_std)

# Models are trained only on healthy controls to capture normal aging
kadrc_hc_index = which(df_meta_kadrc$Diagnosis_group == "HC")
kadrc_prot_hc = df_prot_kadrc_std[kadrc_hc_index, ]
kadrc_prot_with_sex_hc = prot_with_sex[kadrc_hc_index, ]
age_hc = df_meta_kadrc$Age[kadrc_hc_index]
sex_hc = as.numeric(df_meta_kadrc$Sex[kadrc_hc_index] == "F")
cdr_hc = df_meta_kadrc$CDRGLOB[kadrc_hc_index]

kadrc_dc_index = which(df_meta_kadrc$Diagnosis_group != "HC")
kadrc_prot_dc = df_prot_kadrc_std[kadrc_dc_index, ]
kadrc_prot_with_sex_dc = prot_with_sex[kadrc_dc_index, ]
age_dc = df_meta_kadrc$Age[kadrc_dc_index]
sex_dc = as.numeric(df_meta_kadrc$Sex[kadrc_dc_index] == "F")
cdr_dc = df_meta_kadrc$CDRGLOB[kadrc_dc_index]

results <- list()

for (cell_type in cell_type_df$Original.Cell.types) {
  cat("Processing cell type:", cell_type, "\n")
  
  # Get cell-type-specific proteins (Somamers) for this cell type
  somamer_ids = cell_type_df[cell_type_df$Original.Cell.types == cell_type,]$Somamers_list[[1]]
  print(length(somamer_ids))
  
  x_train <- cbind(Sex = sex_hc, kadrc_prot_hc[, somamer_ids])
  
  bootstrap_models <- train_bootstrap_models(x_train, age_hc, alpha_level = 0.5, n_bootstrap=100)
  predictions_hc <- predict_bootstrap_models(bootstrap_models, x_train)
  
  results[[cell_type]] <- list(
    healthy_cohort_pred = predictions_hc,
    healthy_cohort_cor = evaluate_predictions(predictions_hc, age_hc),
    models = bootstrap_models
  )
  
  cat("\nResults for", cell_type, ":\n")
  cat("Correlation:", results[[cell_type]]$healthy_cohort_cor, "\n")
}

saveRDS(results, file = "results/clock_v3.3_update.rds")
