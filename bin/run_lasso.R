#!/usr/bin/env Rscript

# Check and load dependencies
required_packages <- c("argparse", "dplyr", "tidyr", "caret", "glmnet")
invisible(lapply(required_packages, function(p) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) install.packages(p, repos="http://rstudio.com")
}))

library(argparse)
library(dplyr)
library(tidyr)
library(caret)
library(glmnet)

# 1. Parse command line arguments from Nextflow
parser <- ArgumentParser(description="Run Lasso Multinomial Regression on Imputed Legionella Data")
parser\$add_argument("-i", "--input", required=TRUE, help="Path to the imputed input CSV file")
parser\$add_argument("-m", "--model_output", default="lasso_model_summary.txt", help="Output filename for model performance")
parser\$add_argument("-t", "--table_output", default="short_reporting_table.csv", help="Output filename for the clean coefficients table")
args <- parser\$parse_args()

# 2. Read imputed dataset and clean factor variables
# Expects column names matching: Leg_subspecies, flaA, pilE, asd, MIP, mompS, proA, neuA.neuAh
df <- read.csv(args\$input, stringsAsFactors = FALSE)

# Clean subspecies string names into valid R variable names for Caret
df\(Leg_subspecies <- as.factor(make.names(df\)Leg_subspecies))

gene_cols <- c("Gene_flaA", "Gene_pilE", "Gene_asd", "Gene_mip", "Gene_mompS", "Gene_proA", "Gene_Gene_neuA_neuAh
")
df[gene_cols] <- lapply(df[gene_cols], as.factor)

# 3. Fit Lasso model using Stratified 5-fold Cross-Validation via Caret
train_control <- trainControl(
    method = "cv", 
    number = 5, 
    summaryFunction = multiClassSummary,
    classProbs = TRUE
)

lasso_model <- train(
    Leg_subspecies ~ Gene_flA + Gene_pilE + Gene_asd + Gene_MIP + Gene_mompS + Gene_proA + Gene_neuA_neuAh,
    data = df,
    method = "glmnet",
    metric = "Accuracy",
    tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, by = 0.01)),
    trControl = train_control
)

# 4. Save general model fit details and hyperparameters to file
sink(args\$model_output)
cat("=== GLMNET LASSO MULTINOMIAL MODEL TUNING OVERVIEW ===\n\n")
print(lasso_model)
cat("\n=== SELECTED HYPERPARAMETERS ===\n")
cat("Alpha : ", lasso_model\(bestTune\)alpha, "\n")
cat("Lambda: ", lasso_model\(bestTune\)lambda, "\n")
sink()

# 5. Extract sparse coefficients matrix side-by-side (Unified table block)
opt_lambda <- lasso_model\(bestTune\)lambda
raw_coefs <- coef(lasso_model\$finalModel, s = opt_lambda)

# Dynamically extract subspecies class names from the model object list
subspecies_classes <- names(raw_coefs)

extract_clean_df <- function(subsp_name) {
    matrix_data <- as.matrix(raw_coefs[[subsp_name]])
    df_out <- data.frame(
        Feature = rownames(matrix_data),
        Coefficient = round(matrix_data[, 1], 4),
        stringsAsFactors = FALSE
    )
    colnames(df_out) <- c("Feature", subsp_name)
    return(df_out)
}

# Sequentially join the subspecies frames dynamically
joined_table <- extract_clean_df(subspecies_classes[1])
for (subsp in subspecies_classes[-1]) {
    joined_table <- left_join(joined_table, extract_clean_df(subsp), by = "Feature")
}

# Simplify checking vectors for downstream case_when logic
raw_colnames <- colnames(joined_table)
colnames(joined_table) <- c("Feature", "fraseri", "pneumophila", "raphaeli")

# 6. Filter out background noise rows (the lasso zeros filter)
short_table <- joined_table %>%
    filter(!(fraseri == 0 & pneumophila == 0 & raphaeli == 0))

# 7. Apply biological interpretation logic labels
final_reporting_table <- short_table %>%
    mutate(`Decision / Biological Role` = case_when(
        Feature == "(Intercept)" ~ "Baseline Subspecies Prevalence",
        fraseri > 0 & pneumophila == 0 & raphaeli == 0 ~ "Critical Diagnostic Marker for subsp. fraseri",
        fraseri == 0 & pneumophila > 0 & raphaeli == 0 ~ "Critical Diagnostic Marker for subsp. pneumophila",
        fraseri == 0 & pneumophila == 0 & raphaeli > 0 ~ "Critical Diagnostic Marker for subsp. raphaeli",
        fraseri > 0 & pneumophila < 0 ~ "Critical Diagnostic Marker for subsp. fraseri / Negative Marker for subsp. pneumophila",
        TRUE ~ "Secondary Predictive Variant"
    ))

# Standardize final columns to match publication expectations
colnames(final_reporting_table)[1] <- "Genetic Predictor (Locus/Allele)"

# Save output to CSV file
write.csv(final_reporting_table, file = args\$table_output, row.names = FALSE)
