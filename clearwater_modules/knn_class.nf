process KNN_CLASSIFICATION {
    container null
    tag "kNN Classifier on ${scaled_csv.baseName}"
    label 'process_medium'

    // Automatically activates Python and Scikit-Learn tools if using conda
    //conda 'conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy'

    publishDir "${params.output}/knn_classification_results", mode: 'copy'

    input:
    path scaled_csv

    output:
    path "knn_performance_metrics.txt", emit: metrics
    path "knn_performance_summary.txt",   emit: diagnostic_summary
    path "knn_predicted_subspecies.tsv", emit: predictions
    
    script:
    """
    knn_classifier.py \\
        --input ${scaled_csv} \\
        --metrics knn_performance_metrics.txt \\
        --summary knn_performance_summary.txt \\
        --predictions knn_predicted_subspecies.tsv
    """
}
