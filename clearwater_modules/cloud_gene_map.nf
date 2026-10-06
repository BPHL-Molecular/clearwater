process EXTRACT_CLOUD_MAP {
    tag { "Mapping Cloud Genes to Samples" }
    publishDir "${params.output}/pangenome_analysis", mode: 'copy'
    
    // Runs natively in your active terminal environment where pandas is installed
    container null 

    input:
    path roary_outdir

    output:
    path "cloud_genes_sample_map.txt", emit: cloud_map
    path "core_genes_list.txt",       emit: core_list

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np

    # 1. Load the master pangenome table from Roary's output directory
    spreadsheet_file = "${roary_outdir}/gene_presence_absence.csv"
    df = pd.read_csv(spreadsheet_file, low_memory=False)

    # 2. Extract Core Genes list for reference (No. isolates equals total strains)
    # Finding total strains dynamically from the row column counts
    sample_columns = df.columns[14:]
    total_strains = len(sample_columns)

    core_df = df[df['No. isolates'] == total_strains]
    core_df['Gene'].to_csv("core_genes_list.txt", index=False, header=False)

    # 3. Extract Cloud Genes (Present in exactly 1 isolate)
    cloud_df = df[df['No. isolates'] == 1]
    cloud_mapping = []

    for index, row in cloud_df.iterrows():
        gene_name = row['Gene']
        annotation = row['Annotation']
        
        # Check each sample column to see which strain holds this specific unique gene
        for sample in sample_columns:
            if pd.notna(row[sample]) and str(row[sample]).strip() != '':
                cloud_mapping.append({
                    'Gene': gene_name,
                    'Annotation': annotation,
                    'Found_In_Sample': sample
                })
                break

    # 4. Export the clean tracking map to a file
    cloud_mapping_df = pd.DataFrame(cloud_mapping)
    cloud_mapping_df.to_csv("cloud_genes_sample_map.txt", sep="\\t", index=False)
    """
}
