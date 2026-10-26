process EXTRACT_CLOUD_MAP {
    tag { "Mapping Cloud Genes to Samples" }
    publishDir "${params.output}/pangenome_analysis", mode: 'copy'
    
    container null 

    input:
    path roary_outdir
    path global_report
    path metadata_file

    output:
    path "cloud_genes_sample_map.txt", emit: cloud_map
    path "core_genes_list.txt",       emit: core_list

    script:
    // Safe ternary flag evaluation passed downstream into the Python script environment
    def has_metadata = (metadata_file.name != 'EMPTY_METADATA_PLACEHOLDER' && metadata_file.exists()) ? "True" : "False"
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np

    # 1. Load the master pangenome table
    spreadsheet_file = "${roary_outdir}/gene_presence_absence.csv"
    df = pd.read_csv(spreadsheet_file, low_memory=False)

    # 2. Load the global compilation report to get Subspecies mapping
    report_df = pd.read_csv("${global_report}", sep="\\t", low_memory=False)
    subspecies_map = dict(zip(report_df['sampleID'], report_df['subspeciesID_kraken']))

    # 3. Handle the Optional CSV Metadata cross-referencing
    has_meta_flag = ${has_metadata}
    metadata_map = {}
    
    if has_meta_flag:
        # Explicitly tailored to parse a standard Comma-Separated CSV structure
        meta_df = pd.read_csv("${metadata_file}", sep=",", low_memory=False)
        
        # Build dictionary mapping SampleID to metadata traits
        for _, meta_row in meta_df.iterrows():
            s_id = str(meta_row.get('SampleID', '')).strip()
            if s_id:
                metadata_map[s_id] = {
                    'Source': meta_row.get('Isolate_Sources', 'N/A'),
                    'County': meta_row.get('County', 'N/A')
                }

    # 4. Extract Core Genes list for reference
    sample_columns = df.columns[14:]
    total_strains = len(sample_columns)

    core_df = df[df['No. isolates'] == total_strains]
    core_df['Gene'].to_csv("core_genes_list.txt", index=False, header=False)

    # 5. Extract Cloud Genes and join all active arrays
    cloud_df = df[df['No. isolates'] == 1]
    cloud_mapping = []

    for index, row in cloud_df.iterrows():
        gene_name = row['Gene']
        annotation = row['Annotation']
        
        for sample in sample_columns:
            cleaned_sample = str(sample).strip()
            if pd.notna(row[sample]) and cleaned_sample != '':
                subspecies_id = subspecies_map.get(sample, 'Unknown')
                
                # Baseline tracking structure
                payload = {
                    'Gene': gene_name,
                    'Annotation': annotation,
                    'Found_In_Sample': sample,
                    'SubspeciesID': subspecies_id
                }
                
                # Append extra metadata descriptors safely if tracking exists
                if has_meta_flag:
                    sample_meta = metadata_map.get(cleaned_sample, {'Source': 'N/A', 'County': 'N/A'})
                    payload['Isolate_Sources'] = sample_meta['Source']
                    payload['County']          = sample_meta['County']
                
                cloud_mapping.append(payload)
                break

    # 6. Export the final cleanly stitched table (Outputs as a clean TSV report)
    cloud_mapping_df = pd.DataFrame(cloud_mapping)
    cloud_mapping_df.to_csv("cloud_genes_sample_map.txt", sep="\\t", index=False)
    """
}
