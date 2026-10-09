#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import cohen_kappa_score
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description="Run Python Lasso Multinomial Regression on Preprocessed Legionella Data")
    parser.add_argument("-i", "--input", required=True, help="Path to the preprocessed input TSV file")
    parser.add_argument("-m", "--model_output", default="lasso_model_summary.txt", help="Output text summary filename")
    parser.add_argument("-t", "--table_output", default="lasso_coefficients_reportingTable.csv", help="Output short coefficients table filename")
    parser.add_argument("-p", "--plot_output", default="lasso_coefficients_heatmap.png", help="Output plot filename (e.g. .png, .pdf)")
    args = parser.parse_args()

    # 1. Load your preprocessed master report file (expecting tab-separated)
    df = pd.read_csv(args.input, sep="\t")
    
    # Map the long subspecies names to the requested short names immediately
    subspecies_map = {
        "Legionella pneumophila subsp. fraseri": "subsp. fraseri",
        "Legionella pneumophila subsp. pneumophila": "subsp. pneumophila",
        "Legionella pneumophila subsp. raphaeli": "subsp. raphaeli"
    }
    df["subspeciesID_kraken"] = df["subspeciesID_kraken"].map(subspecies_map).fillna(df["subspeciesID_kraken"])

    # Renaming column for code consistency
    df = df.rename(columns={"subspeciesID_kraken": "lp_subspecies"})
    
    # Define your feature columns and target variable name
    gene_cols = ['Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh']
    target_col = 'lp_subspecies'

    # Safety structural check
    missing_cols = [c for c in [target_col] + gene_cols if c not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in input file: {missing_cols}")
        sys.exit(1)

    # 2. Extract features (X) and target (y)
    X = pd.get_dummies(df[gene_cols].astype(str), drop_first=False)
    y = df[target_col]

    # Save mapping for later feature tracking
    feature_names = X.columns.tolist()
    subspecies_classes = sorted(y.unique().tolist())

    # =========================================================================================
    # 3. Setup and Fit LogisticRegressionCV (Lasso L1 Multinomial) with Dynamic Folds
    # =========================================================================================
    # Calculate the size of the smallest sample group to prevent split value crashes
    min_class_size = y.value_counts().min()
    
    # Stratified cross-validation splits cannot be larger than the total members of your rarest class
    if min_class_size > 1:
        safe_cv_folds = min(5, min_class_size)
    else:
        print("Warning: One of your target groups has only 1 sample. Falling back to a 2-fold cross-validation setup.")
        safe_cv_folds = 2

    lasso_model = LogisticRegressionCV(
        cv=safe_cv_folds,
        penalty='l1',
        solver='saga',
        scoring='accuracy',
        max_iter=5000,
        random_state=42,
        Cs=20
    )
    lasso_model.fit(X, y)

    # =========================================================================================
    # 4. Compute Metrics for the Model Performance Metrics Report
    # =========================================================================================
    # Handle multi-class arrays securely by taking the mean across target categories
    raw_C = lasso_model.C_
    if isinstance(raw_C, (np.ndarray, list)):
        best_C = float(np.mean(raw_C))
    else:
        best_C = float(raw_C)
        
    equivalent_lambda = 1.0 / best_C
    alpha = 1.0 / best_C
    
    # Compute accuracy on the full training set
    y_pred = lasso_model.predict(X)
    accuracy = np.mean(y_pred == y)
    
    # Compute Cohen's Kappa score
    kappa = cohen_kappa_score(y, y_pred)

    with open(args.model_output, 'w', encoding='utf-8') as f:
        f.write("=== PYTHON SCIKIT-LEARN LASSO MULTINOMIAL MODEL TUNING OVERVIEW ===\n\n")
        f.write(f"Total Samples Analyzed : {len(df)}\n")
        f.write(f"Total Predictive Feats: {len(feature_names)}\n")
        f.write(f"Target Subspecies      : {', '.join(subspecies_classes)}\n\n")
        
        f.write("=== OPTIMIZED HYPERPARAMETERS & PERFORMANCE METRICS ===\n")
        f.write(f"Alpha                  : {alpha:.6f}\n")
        f.write(f"Lambda (R equivalent)  : {equivalent_lambda:.6f}\n")
        f.write(f"Optimal C Value        : {best_C:.6f}\n")
        f.write(f"Accuracy               : {accuracy * 100:.2f}%\n")
        f.write(f"Cohen's Kappa          : {kappa:.4f}\n")

    # 5. Extract Coefficients side-by-side & Build Short Diagnostic Table
    coef_matrix = lasso_model.coef_
    intercepts = lasso_model.intercept_

    rows = []
    intercept_row = {"Genetic Predictor": "(Intercept)", "Locus/Allele": "N/A"}
    for idx, subsp in enumerate(subspecies_classes):
        intercept_row[subsp] = round(intercepts[idx], 4)
    rows.append(intercept_row)

    # Add Feature variant coefficient rows
    for feat_idx, feat_name in enumerate(feature_names):
        parsed_gene = "Unknown_Locus"
        parsed_allele = feat_name
        
        for g in gene_cols:
            if feat_name.startswith(f"{g}_"):
                parsed_gene = g
                parsed_allele = feat_name[len(g)+1:]
                break

        row = {"Genetic Predictor": parsed_gene, "Locus/Allele": parsed_allele}
        for class_idx, subsp in enumerate(subspecies_classes):
            row[subsp] = round(coef_matrix[class_idx][feat_idx], 4)
        rows.append(row)

    coef_df = pd.DataFrame(rows)
    subsp_cols = subspecies_classes
    
    # 6. Filter out background noise rows (where ALL subspecies coefficients are exactly 0.0)
    short_table = coef_df[
        (coef_df["Genetic Predictor"] == "(Intercept)") |
        (coef_df[subsp_cols].abs().sum(axis=1) > 0.0)
    ].copy()

    # 7. Apply the Biological Decision engine rules dynamically
    def assign_biological_role(row):
        feat = row["Genetic Predictor"]
        if feat == "(Intercept)":
            return "Baseline Subspecies Prevalence"
        
        vals = {s: row[s] for s in subsp_cols}
        positive_classes = [s for s, v in vals.items() if v > 0.0]
        negative_classes = [s for s, v in vals.items() if v < 0.0]

        pos_str = ", ".join(positive_classes)
        neg_str = ", ".join(negative_classes)

        if len(positive_classes) == 1 and len(negative_classes) == 0:
            return f"Critical Diagnostic Marker for {pos_str}"
        elif len(positive_classes) == 1 and len(negative_classes) > 0:
            return f"Critical Diagnostic Marker for {pos_str} / Negative Exclusion Marker for {neg_str}"
        else:
            return "Secondary Predictive Variant"

    short_table["Decision / Biological Role"] = short_table.apply(assign_biological_role, axis=1)

    # 8. Clean up column data formatting
    short_table["Genetic Predictor"] = short_table["Genetic Predictor"].str.strip()
    short_table["Locus/Allele"] = short_table["Locus/Allele"].str.strip()

    # 9. Structure headers cleanly using your exact requested array configuration
    ordered_cols = [
        "Genetic Predictor",
        "Locus/Allele", 
        "subsp. fraseri", 
        "subsp. pneumophila", 
        "subsp. raphaeli", 
        "Decision / Biological Role"
    ]
    
    # Safety cross-check
    for col in ordered_cols:
        if col not in short_table.columns:
            short_table[col] = 0.0 if "subsp." in col else ""

    short_table = short_table[ordered_cols]

    # Export to file
    short_table.to_csv(args.table_output, index=False)

    # =========================================================================================
    # 10. Generate and Save Reporting Visualizations (With Safe Fallback Empty Plot)
    # =========================================================================================
    plot_df = short_table[short_table["Genetic Predictor"] != "(Intercept)"].copy()
    
    if not plot_df.empty:
        plot_df["Feature_Label"] = plot_df["Genetic Predictor"] + ": " + plot_df["Locus/Allele"]
        plot_df.set_index("Feature_Label", inplace=True)
        heatmap_data = plot_df[subsp_cols]
        
        plt.figure(figsize=(8, max(4, len(heatmap_data) * 0.4 + 2)))
        sns.heatmap(
            heatmap_data, 
            annot=True, 
            cmap="RdBu_r", 
            center=0.0, 
            cbar_kws={'label': 'LASSO Log-Odds Coefficient'},
            linewidths=0.5,
            fmt=".4f"
        )
        
        plt.title("LASSO Multinomial Regression: Diagnostic Allele Predictors", fontsize=12, pad=15)
        plt.ylabel("Genetic Predictor (Locus: Allele Variant)", fontsize=10)
        plt.xlabel("Target Subspecies", fontsize=10)
        plt.xticks(rotation=15)
        plt.tight_layout()
        
        plt.savefig(args.plot_output, dpi=300)
        plt.close()
        print(f"Successfully generated outputs: {args.model_output}, {args.table_output}, {args.plot_output}")
    else:
        # Fallback empty image generation block to satisfy Nextflow output requirements
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "No non-zero variants\nfound to plot.", 
                 ha='center', va='center', fontsize=12, color='darkred', weight='bold')
        plt.axis('off')
        plt.tight_layout()
        
        plt.savefig(args.plot_output, dpi=150)
        plt.close()
        print(f"Successfully generated text/table. Saved fallback empty plot layout to {args.plot_output}")

if __name__ == '__main__':
    main()
