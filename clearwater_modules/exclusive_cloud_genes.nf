process PLOT_EXCLUSIVE_CLOUD {
    tag { "Analyzing Exclusive Cloud Gene Profiles" }
    publishDir "${params.output}/pangenome_plots/exclusives", mode: 'copy'
    container null    
    //conda 'conda-forge::biopython=1.81 conda-forge::matplotlib=3.7.1 conda-forge::seaborn=0.12.2 conda-forge::pandas'

    input:
    path all_files

    output:
    path "*.png", emit: distribution_plots
    path "exclusive_*.txt", emit: gene_lists, optional: true

    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    # 1. Locate files flatly in sandbox workspace
    spreadsheet_files = glob.glob("*/gene_presence_absence.csv") or glob.glob("gene_presence_absence.csv")
    imputed_files = glob.glob("imputed_subspecies_report.txt")
    meta_csvs = [f for f in glob.glob("*.csv") if "gene_presence_absence" not in f and "EMPTY_METADATA_PLACEHOLDER" not in f]
    
    spreadsheet_file = spreadsheet_files[0]
    imputed_file = imputed_files[0]
    has_metadata = len(meta_csvs) > 0

    # 2. Extract Binary Cloud Gene Matrix
    roary = pd.read_csv(spreadsheet_file, low_memory=False)
    roary.set_index('Gene', inplace=True)
    roary.drop(list(roary.columns[:13]), axis=1, inplace=True)
    roary.replace('.{2,100}', 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)

    strains = list(roary.columns)
    num_strains = len(strains)
    
    gene_counts = roary.sum(axis=1)
    cloud_genes = gene_counts[gene_counts < (num_strains * 0.15)].index
    roary_cloud = roary.loc[cloud_genes]

    # 3. Parse and Map Groups (Subspecies always runs)
    df_imputed = pd.read_csv(imputed_file, sep="\\t")
    if 'sampleID' in df_imputed.columns:
        df_imputed.rename(columns={'sampleID': 'SampleID'}, inplace=True)
    elif 'sample_id' in df_imputed.columns:
        df_imputed.rename(columns={'sample_id': 'SampleID'}, inplace=True)
    df_imputed['SampleID'] = df_imputed['SampleID'].astype(str)

    subsp_map = {}
    for s in strains:
        match = df_imputed[df_imputed['SampleID'] == str(s)]
        val = match['subspeciesID_kraken'].values[0] if not match.empty and 'subspeciesID_kraken' in df_imputed.columns else "Unknown"
        subsp_map[s] = str(val)

    # 4. EXCLUSIVE PIE CHART 1: Subspecies Partition Allocation
    unique_subsp = sorted(list(set(subsp_map.values())))
    subsp_exclusive_counts = {sub: 0 for sub in unique_subsp}
    subsp_exclusive_lists = {sub: [] for sub in unique_subsp}
    shared_subsp_count = 0

    for gene in cloud_genes:
        present_strains = roary_cloud.columns[roary_cloud.loc[gene] == 1]
        active_groups = set([subsp_map[st] for st in present_strains])
        
        if len(active_groups) == 1:
            grp = list(active_groups)[0]
            subsp_exclusive_counts[grp] += 1
            subsp_exclusive_lists[grp].append(gene)
        else:
            shared_subsp_count += 1

    # Plot Subspecies Pie Chart
    plt.figure(figsize=(8, 8))
    labels = [f"Exclusive to {k}" for k in subsp_exclusive_counts.keys()] + ["Shared across Subsp"]
    sizes = list(subsp_exclusive_counts.values()) + [shared_subsp_count]
    
    plot_data = [(l, s) for l, s in zip(labels, sizes) if s > 0]
    if plot_data:
        p_labels, p_sizes = zip(*plot_data)
        plt.pie(p_sizes, labels=p_labels, autopct=lambda p: '{:.1f}%\\n({:.0f})'.format(p, p*sum(p_sizes)/100),
                startangle=140, colors=sns.color_palette("Pastel1", len(p_sizes)))
        plt.title(f"Distribution of Cloud Genes by Subspecies Specificity\\n(Total Cloud Genes: {len(cloud_genes)})", fontweight='bold', fontsize=12)
        plt.savefig("cloud_genes_by_subspecies_pie.png", dpi=300, bbox_inches='tight')
    plt.clf()

    for sub, glist in subsp_exclusive_lists.items():
        clean_name = sub.replace(" ", "_").replace(".", "").replace("/", "_")
        with open(f"exclusive_subspecies_{clean_name}.txt", "w") as f:
            f.write("\\n".join(glist) + "\\n")

    # 5. EXCLUSIVE PIE CHART 2: Isolate Source (With Strict Validation Guard)
    if has_metadata:
        df_meta = pd.read_csv(meta_csvs[0])
        
        # Guard 1: Verify the literal column string exists
        if 'SampleID' not in df_meta.columns:
            print("WARNING: Provided metadata file lacks a 'SampleID' column. Skipping source distribution evaluation.")
            has_metadata = False
        else:
            df_meta['SampleID'] = df_meta['SampleID'].astype(str)
            meta_samples = set(df_meta['SampleID'].unique())
            run_samples = set(strains)
            
            # Guard 2: Verify keys match active fasta files
            if len(meta_samples.intersection(run_samples)) == 0:
                print("WARNING: Metadata SampleID values do not match active assembly codes. Skipping source distribution evaluation.")
                has_metadata = False

    if has_metadata:
        source_map = {}
        for s in strains:
            match = df_meta[df_meta['SampleID'] == str(s)]
            raw_src = str(match['Isolate_Sources'].values[0]).lower() if not match.empty and 'Isolate_Sources' in df_meta.columns else "unknown"
            
            if any(x in raw_src for x in ['sputum', 'wash', 'aspirate', 'clinical', 'blood', 'patient']):
                source_map[s] = 'Clinical'
            elif any(x in raw_src for x in ['environmental', 'water', 'tower', 'soil', 'swab']):
                source_map[s] = 'Environmental'
            else:
                source_map[s] = 'Other_Unknown'

        unique_sources = sorted(list(set(source_map.values())))
        src_exclusive_counts = {src: 0 for src in unique_sources}
        src_exclusive_lists = {src: [] for src in unique_sources}
        shared_src_count = 0

        for gene in cloud_genes:
            present_strains = roary_cloud.columns[roary_cloud.loc[gene] == 1]
            active_groups = set([source_map[st] for st in present_strains])
            
            if len(active_groups) == 1:
                grp = list(active_groups)[0]
                src_exclusive_counts[grp] += 1
                src_exclusive_lists[grp].append(gene)
            else:
                shared_src_count += 1

        plt.figure(figsize=(8, 8))
        src_labels = [f"Exclusive to {k}" for k in src_exclusive_counts.keys()] + ["Shared across Sources"]
        src_sizes = list(src_exclusive_counts.values()) + [shared_src_count]
        
        src_plot_data = [(l, s) for l, s in zip(src_labels, src_sizes) if s > 0]
        if src_plot_data:
            ps_labels, ps_sizes = zip(*src_plot_data)
            plt.pie(ps_sizes, labels=ps_labels, autopct=lambda p: '{:.1f}%\\n({:.0f})'.format(p, p*sum(ps_sizes)/100),
                    startangle=90, colors=sns.color_palette("Set3", len(ps_sizes)))
            plt.title(f"Distribution of Cloud Genes by Epidemiological Origin\\n(Total Cloud Genes: {len(cloud_genes)})", fontweight='bold', fontsize=12)
            plt.savefig("cloud_genes_by_source_pie.png", dpi=300, bbox_inches='tight')
        plt.clf()

        for src, glist in src_exclusive_lists.items():
            clean_src_name = src.replace("/", "_").replace(" ", "_")
            with open(f"exclusive_source_{clean_src_name}.txt", "w") as f:
                f.write("\\n".join(glist) + "\\n")
    """
}
