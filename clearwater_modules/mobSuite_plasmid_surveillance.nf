process RUN_MOB_TYPER {
    tag "Typing: ${sample_id}"
    
    // Saves raw individual reports to a subfolder to keep your output directory organized
    publishDir "${params.output}/plasmid_analysis/raw_outputs", mode: 'copy'
    container 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "${sample_id}_raw_mob.txt", emit: raw_txt

    script:
    """
    mob_typer \
        --infile ${fasta} \
        --out_file ${sample_id}_raw_mob.txt
    """
}

process MERGE_MOB_REPORTS {
    publishDir "${params.output}/plasmid_analysis", mode: 'copy'
    //container 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0'
    container null
    input:
    // Receives the combined flat array natively without breaking tuple syntax rules
    path all_files

    output:
    path "mobtyper_results.txt", emit: final_report

    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    import pandas as pd

    # Finds all raw mob reports staged natively in the root workspace folder
    files = glob.glob("*_raw_mob.txt")
    all_results = []

    if not files:
        raise FileNotFoundError("No raw mob_typer reports found to merge!")

    for f in files:
        file_name = os.path.basename(f)
        sample_id = file_name.replace("_raw_mob.txt", "")
        
        try:
            df = pd.read_csv(f, sep="\\t")
            if not df.empty:
                df['sample_id'] = sample_id
                all_results.append(df)
        except Exception as e:
            print(f"CRITICAL: Failed to parse execution file {f}. Error: {str(e)}")
            raise e

    final_columns = ["sample_id", "num_contigs", "size", "gc", "rep_type(s)", "relaxase_type(s)", "mpf_type", "predicted_mobility", "primary_cluster_id"]

    if all_results:
        master_df = pd.concat(all_results, ignore_index=True)
        
        for col in final_columns:
            if col not in master_df.columns:
                master_df[col] = "None"
        
        formatted_df = master_df[final_columns]
        
        formatted_df['primary_cluster_id'] = formatted_df['primary_cluster_id'].apply(
            lambda x: f"?? {x}" if pd.notna(x) and str(x).strip() != "" and "??" not in str(x) else x
        )

        # Dynamic clean check for the metadata sheet inside the workspace
        meta_files = [f for f in glob.glob("*.csv") if "EMPTY_METADATA_PLACEHOLDER" not in f]
        
        if meta_files:
            meta_path = meta_files[0] # Grab the staged CSV file name
            try:
                meta_df = pd.read_csv(meta_path, sep=",")
                
                if 'SampleID' in meta_df.columns:
                    cols_to_keep = ['SampleID']
                    if 'Isolate_Sources' in meta_df.columns: cols_to_keep.append('Isolate_Sources')
                    if 'County' in meta_df.columns: cols_to_keep.append('County')
                    
                    meta_df = meta_df[cols_to_keep]

                    formatted_df = pd.merge(
                        formatted_df, 
                        meta_df, 
                        left_on='sample_id', 
                        right_on='SampleID', 
                        how='left'
                    )
                    
                    if 'SampleID' in formatted_df.columns:
                        formatted_df.drop(columns=['SampleID'], inplace=True)
                    
                    if 'Isolate_Sources' in formatted_df.columns: 
                        formatted_df['Isolate_Sources'] = formatted_df['Isolate_Sources'].fillna('None')
                    if 'County' in formatted_df.columns: 
                        formatted_df['County'] = formatted_df['County'].fillna('None')
                else:
                    print("WARNING: Metadata file provided, but column 'SampleID' was missing. Skipping merge.")
            except Exception as meta_err:
                print(f"WARNING: Metadata flag was provided but merging failed: {str(meta_err)}")

    else:
        formatted_df = pd.DataFrame(columns=final_columns)

    formatted_df.to_csv("mobtyper_results.txt", sep="\\t", index=False)
    """
}
