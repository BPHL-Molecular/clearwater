process SKLEARN_LASSO {
    container null
    tag "Lasso Multinomial Regression on ${imputed_report.baseName}"
    label 'process_medium'

    publishDir "${params.output}/lassoMregression_results", mode: 'copy'

    input:
    path imputed_report

    output:
    path "lasso_model_summary.txt",          emit: summary
    path "lasso_coefficients_reportingTable.txt", emit: coefficients_table
    path "lasso_coefficients_heatmap.png",   emit: heatmap, optional: true

    script:
    """
    lasso_multinomial.py \\
        --input ${imputed_report} \\
        --model_output lasso_model_summary.txt \\
        --table_output lasso_coefficients_reportingTable.txt \\
        --plot_output lasso_coefficients_heatmap.png
    """
}