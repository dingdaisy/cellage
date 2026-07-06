# SomaScan7k Cell Aging Clocks (Version 4.1 Assay)

This repository contains the pre-trained models, metadata, performance metrics, and implementation notes corresponding to the **SomaScan7k (Version 4.1 Assay)** cell-type-specific proteomic aging clocks developed in our study. These models were utilized for the GNPC and NSHD cohort results presented in the paper.

---

## Repository Contents

### Data & Metadata Files
*   **`Somalogic_7k_protein_metadata.xlsx`**
    *   **Description:** Contains comprehensive SomaScan metadata for the ~7,000 aptamers, including Gene Symbols, SeqID, and UniProt IDs.
    *   **Usage:** Use the specific **SeqIDs** in this file to map your proteomic features to the model coefficients accurately. 
*   **`Soma_clock_coefficients_min_.xlsx`**
    *   **Description:** Extracted from our Supplementary Tables, this file contains the pre-trained protein weights (coefficients) and intercepts for each of the cell-type-specific aging clocks. Each row corresponds to a specific cell-type model.
*   **`Soma_cell_type_performance_results_.xlsx`**
    *   **Description:** Summary performance table evaluating the Pearson correlation between predicted and actual chronological ages across our Healthy, Disease, and SADRC validation groups (SADRC is an external GNPC sub-cohort).
    *   **Usage:** The table is sorted in descending order of correlation magnitude within the healthy cohort. The total number of feature inputs for each clock is listed in the `n_proteins` column.

### Visualizations & Reference
*   **`visualize_soma_clocks.png`**
    *   **Description:** Visual performance plots (Predicted vs. Actual Age) for 10 representative example models evaluated in healthy individuals within the training cohort.

---

## Training Cohort Context
Our clocks were trained natively on healthy individuals from the **KADRC cohort**. When evaluating performance or applying these models to external data, keep our training demographic baseline in mind:
*   **Mean Age:** $73.2 \pm 10.8$ years
*   **Age Range:** 31 to 101 years

---

## Key Implementation & Preprocessing Notes

### 1. Mandatory Data Preprocessing Pipeline
To remain fully consistent with our cell aging clock training workflow, your input SomaScan **ANML** protein levels must undergo the following transformations *before* applying the model weights:
1.  **$\log_{10}$ Transformation:** Convert the raw ANML values to a logarithmic scale ($\log_{10}$).
2.  **Z-Score Normalization:** Center and scale the log-transformed levels so that each feature has a mean of 0 and a standard deviation of 1.

### 2. Sex Coefficient Encoding
Sex is explicitly included as a feature coefficient in these models. You must encode your cohort's sex covariate as a binary integer:
*   **Female:** `1`
*   **Male:** `0`

### 3. Model Filtering Criteria
While the coefficients sheet contains $>60$ total cell-type models, our core analysis was restricted to **~40 highly performant cell-type models**. To ensure statistical rigor, we filtered out models that failed to meet the following thresholds (unless anchored by a strong prior biological objective):
*   **Training Correlation:** $r \ge 0.25$
*   **Test Correlation:** $r \ge 0.15$
*   **Feature Size:** Minimum of 4 protein features (`n_proteins` $\ge 4$)

The specific ~40 models retained for our study are explicitly annotated inside the `Soma_clock_coefficients_min_.xlsx` sheet.

### Methodology References
For a complete step-by-step walkthrough on downstream analysis and matrix math pipelines, we highly recommend checking out Hamilton Oh's open-source workflows:
*   **SomaScan Organ Aging Tutorial:** Excellent baseline reference for calculating and handling proteomic age gaps (available at [github.com/hamiltonoh/organage](https://github.com/hamiltonoh/organage)).
*   **UK Biobank Organ Clock Tutorial:** Ideal for visualizing exactly how a multi-column coefficient/intercept table is programmatically applied to an external expression matrix to extract biological age estimates.
