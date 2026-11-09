#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report

def main():
    parser = argparse.ArgumentParser(description="Run Robust Dual-File kNN Classification inside Nextflow Pipeline")
    parser.add_argument("-i", "--input", required=True, help="Path to input scaled TSV file")
    parser.add_argument("-m", "--metrics", default="knn_performance_metrics.txt", help="Output standard metrics path")
    parser.add_argument("-s", "--summary", default="knn_performance_summary.txt", help="Output summary warning log path")
    parser.add_argument("-p", "--predictions", default="knn_predicted_subspecies.tsv", help="Output predictions path")
    args = parser.parse_args()

    # 1. Load data and filter/rename columns
    df = pd.read_csv(args.input, sep='\t')
    df = df.rename(columns={'subspeciesID_kraken': 'lp_subspecies'})

    gene_cols = ['Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh']
    target_col = 'lp_subspecies'

    # Filter dataset to rows containing clear classifications
    df = df.dropna(subset=[target_col] + gene_cols).copy()

    # --- UNIVERSAL PIPELINE SAFEGUARDS ---
    class_counts = df[target_col].value_counts()
    unique_classes = sorted(df[target_col].unique().tolist())
    num_classes = len(unique_classes)

    # Handle absolute uniform structures gracefully
    if num_classes <= 1:
        print(f"Critical Warning: Only {num_classes} subspecies class remaining. Writing fallback empty files.")
        error_msg = f"Insufficient comparative variance. Entire batch dataset is uniform: {unique_classes}\n"
        
        # Write clean fallback versions of both metrics sheets
        for path in [args.metrics, args.summary]:
            with open(path, 'w', encoding='utf-8') as f:
                f.write("======================================================================\n")
                f.write("         kNN CLASSIFIER TERMINATED GRACEFULLY                        \n")
                f.write("======================================================================\n\n")
                f.write(error_msg)
        
        df['knn_predicted_subspecies'] = df[target_col] if num_classes == 1 else "Unknown"
        df.to_csv(args.predictions, sep='\t', index=False)
        sys.exit(0)

    X = df[gene_cols]
    y = df[target_col]
    min_class_size = class_counts.min()

    # --- DYNAMIC SPLITTING & VALIDATION SWITCH ---
    forced_resubstitution = False
    resubstitution_reason = ""

    # Stratified 70/30 splits require at least 4 items per class to populate both partitions safely.
    if min_class_size < 4 or len(df) < 15:
        forced_resubstitution = True
        resubstitution_reason = (
            f"Extreme class imbalance or tiny sample cohort detected (Minimum class size = {min_class_size}, Total N = {len(df)}). "
            f"A stratified data split requires at least 4 samples per category to allocate across training/testing partitions "
            f"without causing structural calculation failures. Evaluation metrics reflect full-cohort training set resubstitution."
        )
        print(f"Warning: {resubstitution_reason}")
        
        X_train, X_test, y_train, y_test = X, X, y, y
        validation_strategy = "Full Cohort Resubstitution"
        # If resubstitution mode runs, we preserve the entire row frame array indices for predictions
        prediction_df = df.copy()
    else:
        # Standard robust partitioned data track when sufficient samples exist
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30, random_state=42, stratify=y)
        validation_strategy = "70% Training / 30% Testing Partitioned Split"
        # Extract matching dataframe subset tracking precisely back to the true validation block row index map
        prediction_df = df.loc[y_test.index].copy()

    # Dynamic K Neighbors parameter boundary calculations
    dynamic_k = min(5, len(X_train) - 1)
    if dynamic_k < 1:
        dynamic_k = 1

    # 3. Setup and fit kNN Model
    knn = KNeighborsClassifier(n_neighbors=dynamic_k)
    knn.fit(X_train, y_train)

    # Predict targets and lock prediction profiles into the final output dataframe 
    y_pred = knn.predict(X_test)
    prediction_df['knn_predicted_subspecies'] = y_pred
    prediction_df.to_csv(args.predictions, sep='\t', index=False)

    # 4. Compute metrics and structure layout tables dynamically
    acc = accuracy_score(y_test, y_pred) * 100
    test_labels = sorted(y_test.unique().tolist())

    cm = confusion_matrix(y_test, y_pred, labels=test_labels)
    report = classification_report(y_test, y_pred, target_names=test_labels)

    # Generate structured matrix alignments
    max_label_len = max(len(str(lbl)) for lbl in test_labels)
    col_width = max(max_label_len, 12)

    matrix_header = f"{'True / Pred':<{col_width}} | " + " | ".join(f"{lbl:<{col_width}}" for lbl in test_labels)
    divider = "-" * len(matrix_header)

    # ==================================================================
    # OUTPUT FILE 1: FILE SYSTEM REPORT (knn_performance_metrics.txt)
    # ==================================================================
    with open(args.metrics, 'w', encoding='utf-8') as f:
        f.write("======================================================================\n")
        f.write("         kNN SUBSPECIES CLASSIFIER PERFORMANCE SUMMARY                \n")
        f.write("======================================================================\n\n")
        f.write(f"Partitioning Strategy : {validation_strategy}\n")
        f.write(f"Overall Accuracy      : {acc:.2f}%\n\n")
        f.write("--- DETAILED EVALUATION METRICS (TEST SET) ---\n")
        f.write(report + "\n\n")
        f.write("--- ALIGNED CONFUSION MATRIX ---\n")
        f.write(matrix_header + "\n")
        f.write(divider + "\n")
        for i, row_label in enumerate(test_labels):
            row_str = " | ".join(f"{val:<{col_width}}" for val in cm[i])
            f.write(f"{row_label:<{col_width}} | {row_str}\n")
        f.write("======================================================================\n")

    # ==================================================================
    # OUTPUT FILE 2: FILE SYSTEM REPORT (knn_performance_summary.txt)
    # ==================================================================
    with open(args.summary, 'w', encoding='utf-8') as f:
        f.write("======================================================================\n")
        f.write("         kNN SUBSPECIES CLASSIFIER DIAGNOSTIC LOG SUMMARY            \n")
        f.write("======================================================================\n\n")
        f.write(f"Partitioning Strategy : {validation_strategy}\n")
        f.write(f"K-Neighbors Value     : {dynamic_k}\n")
        f.write(f"Overall Accuracy      : {acc:.2f}%\n\n")
        
        if forced_resubstitution:
            f.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
            f.write("WARNING: FORCED RESUBSTITUTION MODE ACTIVE\n")
            f.write(f"Reason: {resubstitution_reason}\n")
            f.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
        else:
            f.write("STATUS SUCCESS: Sufficient comparative variance verified across group channels.\n\n")
            
        f.write("--- DETAILED EVALUATION METRICS ---\n")
        f.write(report + "\n")
        f.write("======================================================================\n")

if __name__ == '__main__':
    main()
