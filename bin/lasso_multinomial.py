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
    parser = argparse.ArgumentParser(description="Universal Robust Lasso Multinomial Regression for Legionella Data")
    parser.add_argument("-i", "--input", required=True, help="Path to preprocessed input TSV file")
    parser.add_argument("-m", "--model_output", default="lasso_model_summary.txt", help="Output text summary filename")
    parser.add_argument("-t", "--table_output", default="lasso_coefficients_reportingTable.txt", help="Output short coefficients table filename")
    parser.add_argument("-p", "--plot_output", default="lasso_coefficients_heatmap.png", help="Output plot filename")
    args = parser.parse_args()

    # 1. Load input dataset
    df = pd.read_csv(args.input, sep="\t")
    
    # Standardize and map long names immediately
    subspecies_map = {
        "Legionella pneumophila subsp. fraseri": "subsp. fraseri",
        "Legionella pneumophila subsp. pneumophila": "subsp. pneumophila",
        "Legionella pneumophila subsp. raphaeli": "subsp. raphaeli",
        "Legionella pneumophila subsp. pascullei": "subsp. pascullei"
    }
    df["subspeciesID_kraken"] = df["subspeciesID_kraken"].map(subspecies_map).fillna(df["subspeciesID_kraken"])
    df = df.rename(columns={"subspeciesID_kraken": "lp_subspecies"})
    
    gene_cols = ['Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh']
    target_col = 'lp_subspecies'

    # Safety structural check
    missing_cols = [c for c in [target_col] + gene_cols if c not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in input file: {missing_cols}")
        sys.exit(1)

    # --- UNIVERSAL PIPELINE SAFEGUARDS ---
    # Count frequency of each active subspecies present in this unique batch run
    class_counts = df[target_col].value_counts()
    
    # Extract exactly what unique classes survived to represent this specific batch
    subspecies_classes = sorted(df[target_col].unique().tolist())
    num_classes = len(subspecies_classes)

    # CASE 1: Uniform Cohort or Insufficient Classes (Survival Mode)
    if num_classes <= 1:
        uniform_name = subspecies_classes[0] if num_classes == 1 else "Unknown/Uniform"
        print(f"Critical Warning: Only {num_classes} subspecies class remaining ({subspecies_classes}). Writing fallback files safely.")
        
        note_msg = f"Insufficient comparative variance. Entire batch dataset is uniform: {uniform_name}"
        
        # 1. Write to model summary file
        with open(args.model_output, 'w', encoding='utf-8') as f:
            f.write(f"=== PYTHON SCIKIT-LEARN LASSO MODEL SKIPPED ===\n")
            f.write(f"Reason: {note_msg}\n")
        
        # 2. Write dynamically configured headers to reporting table
        ordered_cols = ["Genetic Predictor", "Locus/Allele"] + subspecies_classes + ["Decision / Biological Role"]
        fallback_df = pd.DataFrame(columns=ordered_cols)
        fallback_df.to_csv(args.table_output, sep="\t", index=False)
        
        # 3. Append the clear note to the bottom of the reporting table
        with open(args.table_output, 'a', encoding='utf-8') as f:
            f.write(f"\n# Note: {note_msg}\n")
            
        sys.exit(0)

    # 2. Extract features (X) and target (y) without discarding any rare categories
    X = pd.get_dummies(df[gene_cols].astype(str), drop_first=False)
    y = df[target_col]
    feature_names = X.columns.tolist()

    min_class_size = class_counts.min()
    
    # 3. Dynamic Modeling Strategy Selector
    if min_class_size < 3:
        print(f"Warning: Smallest class size is {min_class_size} (< 3). Stratified CV is impossible.")
        print("Pipeline automatically switching to Resubstitution Mode to protect sparse arrays.")
        
        from sklearn.linear_model import LogisticRegression
        reporting_C = 1.0
        chosen_solver = 'liblinear' if num_classes == 2 else 'saga'
        
        lasso_model = LogisticRegression(
            penalty='l1', solver=chosen_solver, max_iter=5000, 
            random_state=42, C=reporting_C
        )
        lasso_model.fit(X, y)
        cv_status_message = "SKIPPED (Forced Resubstitution due to sparse class sizes < 3)"
        
    else:
        dynamic_cv = min(5, min_class_size - 1)
        if dynamic_cv < 2:
            dynamic_cv = 2
            
        print(f"Running Lasso model with {num_classes} classes. Dynamic CV set to: {dynamic_cv}")
        chosen_solver = 'liblinear' if num_classes == 2 else 'saga'
        
        lasso_model = LogisticRegressionCV(
            cv=dynamic_cv, penalty='l1', solver=chosen_solver,
            scoring='accuracy', max_iter=5000, random_state=42, Cs=20
        )
        lasso_model.fit(X, y)
        
        best_C = lasso_model.C_
        reporting_C = float(best_C[0]) if isinstance(best_C, np.ndarray) else float(best_C)
        cv_status_message = f"SUCCESS ({dynamic_cv}-fold cross-validation optimized)"                                                                       
    
    equivalent_lambda = 1.0 / reporting_C
    y_pred = lasso_model.predict(X)
    accuracy = np.mean(y_pred == y)
    kappa = cohen_kappa_score(y, y_pred)
    alpha = 1.0 / reporting_C

    # 4. Write the Model Summary
    with open(args.model_output, 'w', encoding='utf-8') as f:
        f.write("=== PYTHON SCIKIT-LEARN LASSO MULTINOMIAL MODEL TUNING OVERVIEW ===\n\n")
        f.write(f"Total Samples Analyzed : {len(df)}\n")
        f.write(f"Total Predictive Feats: {len(feature_names)}\n")
        f.write(f"Target Subspecies      : {', '.join(subspecies_classes)}\n\n")
        f.write("=== OPTIMIZED HYPERPARAMETERS & PERFORMANCE METRICS ===\n")
        f.write(f"Alpha                  : {alpha:.6f}\n")
        f.write(f"Lambda (R equivalent)  : {equivalent_lambda:.6f}\n")
        f.write(f"Optimal C Value        : {reporting_C:.6f}\n")
        f.write(f"Accuracy               : {accuracy * 100:.2f}%\n")
        f.write(f"Cohen's Kappa          : {kappa:.4f}\n")

    # 5. Extract Coefficients safely
    coef_matrix = lasso_model.coef_
    intercepts = lasso_model.intercept_

    if num_classes == 2:
        coef_matrix = np.vstack([-coef_matrix, coef_matrix])
        flat_intercept = intercepts[0] if isinstance(intercepts, np.ndarray) and intercepts.ndim > 0 else intercepts
        intercepts = np.array([-flat_intercept, flat_intercept])

    rows = []
    intercept_row = {"Genetic Predictor": "(Intercept)", "Locus/Allele": "N/A"}
    for idx, subsp in enumerate(subspecies_classes):
        val = intercepts[idx]
        if hasattr(val, "item"):
            val = val.item()
        intercept_row[subsp] = round(float(val), 4)
    rows.append(intercept_row)

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
            val_coef = coef_matrix[class_idx][feat_idx]
            if hasattr(val_coef, "item"):
                val_coef = val_coef.item()
            row[subsp] = round(float(val_coef), 4)
        rows.append(row)

    coef_df = pd.DataFrame(rows)

    # 6. Filter out background noise rows
    short_table = coef_df[
        (coef_df["Genetic Predictor"] == "(Intercept)") |
        (coef_df[subspecies_classes].abs().sum(axis=1) > 0.0)
    ].copy()

    # 7. Apply the Biological Decision engine rules dynamically
    def assign_biological_role(row):
        feat = row["Genetic Predictor"]
        if feat == "(Intercept)":
            return "Baseline Subspecies Prevalence"
        
        vals = {s: row[s] for s in subspecies_classes}
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
    short_table["Genetic Predictor"] = short_table["Genetic Predictor"].str.strip()
    short_table["Locus/Allele"] = short_table["Locus/Allele"].str.strip()

    # 9. Structure headers dynamically based ONLY on exactly what subspecies exist in this dataset run
    ordered_cols = ["Genetic Predictor", "Locus/Allele"] + subspecies_classes + ["Decision / Biological Role"]
    short_table = short_table[ordered_cols]

    # Save the dataframe
    short_table.to_csv(args.table_output, sep="\t", index=False)
    
    # Determine the status message based on data variance metrics
    if num_classes <= 1:
        note_msg = f"Insufficient comparative variance. Entire batch dataset is uniform: {subspecies_classes[0]}"
    elif kappa == 0.0:
        note_msg = "Warning: Model failed to find stable genetic predictors (Kappa=0.00). Coefficients reflect background baseline prevalence only."
    else:
        note_msg = "Model status: Completed successfully with valid feature selection iterations."

    # Append note to the bottom of the table file
    with open(args.table_output, 'a', encoding='utf-8') as f:
        f.write(f"\n# {note_msg}\n")

    # Append note to the bottom of the summary text file
    with open(args.model_output, 'a', encoding='utf-8') as f:
        f.write(f"\n=== SYSTEM NOTICE ===\n{note_msg}\n")

    # 10. Save Reporting Visualizations Safely
    plot_df = short_table[short_table["Genetic Predictor"] != "(Intercept)"].copy()
    if not plot_df.empty:
        plot_df["Feature_Label"] = plot_df["Genetic Predictor"] + ": " + plot_df["Locus/Allele"]
        plot_df.set_index("Feature_Label", inplace=True)
        heatmap_data = plot_df[subspecies_classes]
        
        plt.figure(figsize=(8, max(4, len(heatmap_data) * 0.4 + 2)))
        sns.heatmap(heatmap_data, annot=True, cmap="RdBu_r", center=0.0, cbar_kws={'label': 'LASSO Log-Odds Coefficient'}, linewidths=0.5, fmt=".4f")
        plt.xticks(rotation=15)
        plt.tight_layout()
        plt.savefig(args.plot_output, dpi=300)
        plt.close()
        print(f"[SUCCESS] Generated LASSO metrics outputs cleanly.")
    else:
        print("No non-zero variants found to plot.")

if __name__ == '__main__':
    main()
