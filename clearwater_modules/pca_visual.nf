process PCA_VISUALIZATION {
    container null 
    tag "PCA Panel Analysis"
    label 'process_medium'

    // Recommended container or conda setup for reproducibility
    // container 'biocontainers/pandas:1.5.1' 
    // conda 'conda-forge::pandas conda-forge::scikit-learn conda-forge::matplotlib conda-forge::seaborn'

    publishDir "${params.output}/pca_plots", mode: 'copy'

    input:
    path full_dataset      // Expected columns: lp_subspecies + 7 Gene columns
    path knn_predictions   // Expected columns: knn_predicted_subspecies + 7 Gene columns

    output:
    path "subspecies_pca_comparison_panel.png", emit: pca_plot

    shell:
    '''
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA

    # Define the 7 core allele profiling columns
    gene_cols = ['Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS', 'Gene_proA', 'Gene_neuA_neuAh']

    # 1. Load both data matrices
    df_full = pd.read_csv('!{full_dataset}', sep='\t')
    df_knn  = pd.read_csv('!{knn_predictions}', sep='\t')

    # 2. Fit PCA on the FULL dataset to establish baseline variance space
    pca = PCA(n_components=2)
    pc_coords_full = pca.fit_transform(df_full[gene_cols])
    
    # 3. Project the kNN test dataset onto that EXACT same PCA space
    pc_coords_knn = pca.transform(df_knn[gene_cols])

    # Assign coordinates back to respective DataFrames
    df_full['PC1'] = pc_coords_full[:, 0]
    df_full['PC2'] = pc_coords_full[:, 1]
    
    df_knn['PC1'] = pc_coords_knn[:, 0]
    df_knn['PC2'] = pc_coords_knn[:, 1]

    # Calculate variance explained percentage based on the full matrix baseline
    var_exp = pca.explained_variance_ratio_ * 100

    # 4. Construct Dual-Panel Plot Layout with uniform axes limits
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharex=True, sharey=True)

    # Panel A: True Subspecies Mapping (Full Dataset)
    sns.scatterplot(
        ax=axes[0], data=df_full, x='PC1', y='PC2', 
        hue='lp_subspecies', palette='Set1', style='lp_subspecies', s=70, alpha=0.8
    )
    axes[0].set_title('A. Full Dataset Variance: True Subspecies ID', fontsize=12, pad=10)
    axes[0].set_xlabel(f'PC1 ({var_exp[0]:.2f}% Var)', fontsize=10)
    axes[0].set_ylabel(f'PC2 ({var_exp[1]:.2f}% Var)', fontsize=10)
    axes[0].legend(title='True Subspecies', loc='best')
    axes[0].grid(True, linestyle='--', alpha=0.5)

    # Panel B: kNN Predicted Subspecies Mapping (Test Split)
    sns.scatterplot(
        ax=axes[1], data=df_knn, x='PC1', y='PC2', 
        hue='knn_predicted_subspecies', palette='Set1', style='knn_predicted_subspecies', s=70, alpha=0.8
    )
    axes[1].set_title('B. Test Dataset Predictions: kNN Classification', fontsize=12, pad=10)
    axes[1].set_xlabel(f'PC1 ({var_exp[0]:.2f}% Var)', fontsize=10)
    axes[1].legend(title='kNN Predicted', loc='best')
    axes[1].grid(True, linestyle='--', alpha=0.5)

    plt.suptitle('Principal Component Analysis (PCA) Comparison of Legionella pneumophila 7-Gene Profiles', fontsize=14, y=0.98)
    plt.tight_layout()

    # Save final graphical panel asset
    plt.savefig('subspecies_pca_comparison_panel.png', dpi=300)
    plt.close()
    '''
}
