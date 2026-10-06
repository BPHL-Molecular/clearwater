
process SKLEARN_LASSO {
    container null
    tag "Scikit-Learn Lasso Logistic Regression on ${imputed_csv.baseName}"
    label 'process_medium'

    // Configures Conda requirements automatically for this step
    // conda 'conda-forge::pandas conda-forge::scikit-learn conda-forge::numpy conda-forge::matplotlib conda-forge::seaborn'
    
    publishDir "${params.output}/lassoMregression_results", mode: 'copy'

    input:
    path imputed_csv

    output:
    path "lasso_model_summary.txt", emit: model_summary
    path "lasso_coefficients_reportingTable.csv", emit: short_table
    path "lasso_coefficients_heatmap.png", emit: coef_plot 

    script:
    """
    run_lasso.py \\
        --input ${imputed_csv} \\
        --model_output lasso_model_summary.txt \\
        --table_output lasso_coefficients_reportingTable.csv \\
        --plot_output lasso_coefficients_heatmap.png
    """
}
