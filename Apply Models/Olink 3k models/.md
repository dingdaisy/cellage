# Olink Explore 3072 Cell and Lineage Aging Clocks

This repository contains the pre-trained models, metadata, performance metrics, and reference scripts corresponding to the **Olink Explore 3072 (Olink3k)** cell-type-specific and lineage-type proteomic aging clocks developed in our study. These models were utilized to generate the UK Biobank cohort results presented in the accompanying paper.

---

## File Contents

### Data & Metadata Files
* **`Olink_3k_protein_metadata.xlsx`**
    * **Description:** Contains comprehensive metadata for the ~3,000 Olink-measured proteins, including Assay Name, UniProt ID, Panel, Protein Name, and target Gene Name.
    * **Usage:** Use this file to map your cohort's proteomic features to the model coefficients. Your measured features must map accurately to these target IDs to apply the clocks.
* **`olink_clock_coefficients_min.xlsx`**
    * **Description:** Extracted from the study's Supplementary Tables, this file contains the pre-trained model weights (coefficients) and intercepts for all cell-type and lineage-type aging clocks. Each row corresponds to a specific model, and columns represent individual proteins or the model intercept.
* **`cell_type_performance_results_final.xlsx`**
    * **Description:** Provides detailed training and test performance metrics—including Pearson correlation ($r$), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE)—evaluated within our UK Biobank cohort. 
    * **Usage:** Use the `Source` column to filter between cell-type or lineage-type models and to prioritize the most performant cell lines for your downstream analyses.

### Scripts & Visualization
* **`demo_apply_cell_clock.R`**
    * **Description:** A reference R script demonstrating how to load the pre-trained weights and intercepts and apply them to an external validation cohort.
* **`clock_visualization_full_cohort.png`**
    * **Description:** Performance scatter plots (Predicted vs. Chronological Age) with computed Pearson correlations for 10 representative models across the full UK Biobank cohort. *Note: Sex-specific cell types are included here for reference, though they were not explicitly analyzed in the core study.*

---

## Key Implementation Notes

### 1. Model Filtering & Selection
While `olink_clock_coefficients_min.xlsx` contains weights for ~70 total models, our core study focused strictly on **~40 cell types and lineages of high interest**. To ensure statistical robustness, we generally excluded models that did not meet the following quality thresholds unless there was a strong prior biological justification:
* **Training Correlation:** $r \ge 0.25$
* **Test Correlation:** $r \ge 0.15$
* **Feature Size:** Minimum of 4 protein features

The specific models retained for analysis in our study are explicitly annotated inside the `olink_clock_coefficients_min.xlsx` sheet.

### 2. Cell-Type vs. Lineage-Type Models
Due to variable protein coverage of certain specific cell-mapped proteins in the Olink Explore 3072 assay, we introduced broader **lineage models** in the UK Biobank cohort. 
* Lineage models serve as a composite representation of grouped cell types.
* Because the feature selection process identifies unique proteomic signatures for each, **we highly recommend evaluating both cell-type and lineage-type models** in your parallel analyses.

### 3. Calculating Age Gaps in External Cohorts
External validation cohorts often exhibit unique structural variations, distinct age distributions, and batch effects. When calculating the "Age Gap" ($\text{Predicted Age} - \text{Chronological Age}$), keep the following in mind:
* **Healthy Baseline:** If your cohort contains a clearly defined "healthy control" subpopulation, we recommend fitting your reference aging curve (using a **LOWESS** regression) on those healthy individuals *only*. You can then project the full cohort onto this baseline to extract cleaner biological signals.
* **Methodology Reference:** For a comprehensive conceptual framework on handling proteomic age gaps, refer to the GitHub organ-aging tutorial by **Hamilton Oh**. Note that while his repository utilizes linear models for the baseline reference curve, this study opts for **LOWESS** smoothing to capture non-linear kinetics.
