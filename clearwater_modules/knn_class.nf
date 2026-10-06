process KNN_CLASSIFICATION {
    container null
    tag "kNN Classifier on ${scaled_csv.baseName}"
    label 'process_medium'

    // Automatically activates Python and Scikit-Learn tools
    //conda 'conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy'

    publishDir "${params.output}/knn_classification_results", mode: 'copy'

    input:
    path scaled_csv

    output:
    path "knn_performance_metrics.txt", emit: metrics
    path "knn_predicted_subspecies.tsv", emit: predictions

    script:
    """
    # Force the environment's python3 to process this entire block as raw Python text
    python3 << 'EOF'
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report

# 1. Load data and filter/rename columns
df = pd.read_csv('${scaled_csv}', sep='\t')
df = df.rename(columns={'subspeciesID_kraken': 'lp_subspecies'})

gene_cols = ['Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh']
target_col = 'lp_subspecies'

# Filter dataset to rows containing clear classifications
df = df.dropna(subset=[target_col] + gene_cols)

X = df[gene_cols]
y = df[target_col]

# 2. Partitions data (70% training, 30% testing)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30, random_state=42, stratify=y)

# 3. Setup and fit kNN Model (k=5 as a reliable baseline)
knn = KNeighborsClassifier(n_neighbors=5)
knn.fit(X_train, y_train)

# Predict targets on the testing partition
y_pred = knn.predict(X_test)
# Calculate whole-dataset projections for reporting
df['knn_predicted_subspecies'] = knn.predict(X)
df.to_csv('knn_predicted_subspecies.tsv', sep='\t', index=False)

# 4. Compute metrics and structure layout tables DYNAMICALLY
acc = accuracy_score(y_test, y_pred) * 100

# DYNAMIC FIX: Find unique classes that actually exist inside y_test
test_labels = sorted(y_test.unique().tolist())

# DYNAMIC FIX: Pass the exact test labels to both functions
cm = confusion_matrix(y_test, y_pred, labels=test_labels)
report = classification_report(y_test, y_pred, target_names=test_labels)

# Generate a beautiful, structured matrix table based on test_labels
max_label_len = max(len(str(lbl)) for lbl in test_labels)
col_width = max(max_label_len, 12)

matrix_header = f"{'True / Pred':<{col_width}} | " + " | ".join(f"{lbl:<{col_width}}" for lbl in test_labels)
divider = "-" * len(matrix_header)

# Write output file
with open('knn_performance_metrics.txt', 'w') as f:
    f.write('======================================================================\\n')
    f.write('         kNN SUBSPECIES CLASSIFIER PERFORMANCE SUMMARY                \\n')
    f.write('======================================================================\\n\\n')
    f.write(f'Partitioning Strategy : 70% Training / 30% Testing\\n')
    f.write(f'Overall Accuracy      : {acc:.2f}%\\n\\n')

    f.write('--- DETAILED EVALUATION METRICS (TEST SET) ---\\n')
    f.write(report)
    f.write('\\n')

    f.write('--- ALIGNED CONFUSION MATRIX ---\\n')
    f.write(f'{matrix_header}\\n')
    f.write(f'{divider}\\n')
    for i, row_label in enumerate(test_labels):
        row_str = " | ".join(f"{val:<{col_width}}" for val in cm[i])
        f.write(f'{row_label:<{col_width}} | {row_str}\\n')
    f.write('======================================================================\\n')
EOF
    """
}
