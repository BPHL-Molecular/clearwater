process PLOT_CLOUD_PANGENOME {
    tag { "Plotting Dynamic Structured Cloud Matrix" }
    publishDir "${params.output}/pangenome_plots", mode: 'copy'
    container null 
    //conda 'conda-forge::biopython=1.81 conda-forge::matplotlib=3.7.1 conda-forge::seaborn=0.12.2 conda-forge::pandas'

    input:
    path all_files

    output:
    path "cloud_pangenome_annotated.png", emit: cloud_plot

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
    from Bio import Phylo
    from matplotlib.patches import Patch

    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    sns.set_style('white')

    # 1. Locate files flatly in sandbox workspace
    tree_files = glob.glob("*/accessory_binary_genes.fa.newick") or glob.glob("*.newick")
    spreadsheet_files = glob.glob("*/gene_presence_absence.csv") or glob.glob("gene_presence_absence.csv")
    imputed_files = glob.glob("imputed_subspecies_report.txt")
    
    if not tree_files or not spreadsheet_files or not imputed_files:
        raise FileNotFoundError("CRITICAL: Missing core files inside sandbox directory wrapper!")

    tree_file = tree_files[0]
    spreadsheet_file = spreadsheet_files[0]
    imputed_file = imputed_files[0]

    meta_csvs = [f for f in glob.glob("*.csv") if "gene_presence_absence" not in f and "EMPTY_METADATA_PLACEHOLDER" not in f]
    has_metadata = len(meta_csvs) > 0

    # Load phylogenetic layout
    t = Phylo.read(tree_file, 'newick')
    strains_ordered = [x.name for x in t.get_terminals()]
    num_strains = len(strains_ordered)

    # 2. Extract and Filter Cloud Genes (< 15% presence across strains)
    roary = pd.read_csv(spreadsheet_file, low_memory=False)
    roary.set_index('Gene', inplace=True)
    roary.drop(list(roary.columns[:13]), axis=1, inplace=True)

    roary.replace('.{2,100}', 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)

    gene_counts = roary.sum(axis=1)
    cloud_genes = gene_counts[gene_counts < (num_strains * 0.15)].index

    if len(cloud_genes) == 0:
        raise ValueError("CRITICAL: No cloud genes found matching frequency boundaries (<15%)!")

    roary_cloud = roary.loc[cloud_genes, strains_ordered]

    # 3. Parse Subspecies Data with Exact Column Mapping
    df_imputed = pd.read_csv(imputed_file, sep="\\t")
    if 'sampleID' in df_imputed.columns:
        df_imputed.rename(columns={'sampleID': 'SampleID'}, inplace=True)
    elif 'sample_id' in df_imputed.columns:
        df_imputed.rename(columns={'sample_id': 'SampleID'}, inplace=True)
        
    if 'SampleID' not in df_imputed.columns:
        raise KeyError(f"CRITICAL: Could not find identification key in imputed file! Found: {list(df_imputed.columns)}")

    df_imputed['SampleID'] = df_imputed['SampleID'].astype(str)

    subspecies_vector = []
    for strain in strains_ordered:
        sub_match = df_imputed[df_imputed['SampleID'] == str(strain)]
        subspecies = sub_match['subspeciesID_kraken'].values[0] if not sub_match.empty and 'subspeciesID_kraken' in df_imputed.columns else "Unknown"
        subspecies_vector.append(str(subspecies))

    unique_subsp = sorted(list(set(subspecies_vector)))
    subsp_ints = np.array([unique_subsp.index(s) for s in subspecies_vector]).reshape(-1, 1)

    # 4. Parse Metadata Columns (With Strict SampleID Validation Guard)
    if has_metadata:
        df_meta = pd.read_csv(meta_csvs[0])
        
        # Guard 1: Verify the literal column string exists
        if 'SampleID' not in df_meta.columns:
            print("WARNING: Provided metadata file lacks a 'SampleID' column. Skipping source mapping ribbon.")
            has_metadata = False
        else:
            df_meta['SampleID'] = df_meta['SampleID'].astype(str)
            meta_samples = set(df_meta['SampleID'].unique())
            run_samples = set(strains_ordered)
            
            # Guard 2: Verify keys match active fasta files
            if len(meta_samples.intersection(run_samples)) == 0:
                print("WARNING: Metadata SampleID values do not match active assembly codes. Skipping source mapping ribbon.")
                has_metadata = False

    if has_metadata:
        sources_vector = []
        for strain in strains_ordered:
            meta_match = df_meta[df_meta['SampleID'] == str(strain)]
            source = meta_match['Isolate_Sources'].values[0] if not meta_match.empty and 'Isolate_Sources' in df_meta.columns else "Unknown"
            sources_vector.append(str(source))
            
        unique_sources = sorted(list(set(sources_vector)))
        source_ints = np.array([unique_sources.index(s) for s in sources_vector]).reshape(-1, 1)

    # 5. Create the Multi-Panel Grid Visualization
    fig = plt.figure(figsize=(24, 12))
    
    if has_metadata:
        ax_tree = plt.subplot2grid((1, 50), (0, 0), colspan=7, facecolor='white')
        ax_matrix = plt.subplot2grid((1, 50), (0, 8), colspan=30)
        ax_source = plt.subplot2grid((1, 50), (0, 39), colspan=2)
        ax_subsp = plt.subplot2grid((1, 50), (0, 42), colspan=2)
    else:
        ax_tree = plt.subplot2grid((1, 48), (0, 0), colspan=7, facecolor='white')
        ax_matrix = plt.subplot2grid((1, 48), (0, 8), colspan=30)
        ax_subsp = plt.subplot2grid((1, 48), (0, 39), colspan=2)

    # A. Draw aligned phylogenetic guide tree
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])
    Phylo.draw(t, axes=ax_tree, show_confidence=False, label_func=lambda x: None,
               xticks=([],), yticks=([],), ylabel=('',), xlabel=('',),
               xlim=(-mdist*0.1, mdist+mdist*0.1), axis=('off',), do_show=False)
               
    ax_tree.set_title('Tree\\n(%d strains)' % num_strains, pad=15, fontweight='bold', fontsize=12)

    # B. Draw Main Heatmap Matrix
    ax_matrix.matshow(roary_cloud.T, cmap=plt.cm.Blues, vmin=0, vmax=1, aspect='auto', interpolation='none')
    ax_matrix.set_yticks([])
    ax_matrix.set_xticks([])
    ax_matrix.set_title('Cloud Matrix\\n(%d gene clusters)' % len(cloud_genes), pad=15, fontweight='bold', fontsize=12)

    # C. Draw Source Ribbon conditionally
    legend_elements = []
    if has_metadata:
        cmap_src = plt.cm.get_cmap('Set3', len(unique_sources))
        ax_source.matshow(source_ints, cmap=cmap_src, aspect='auto', interpolation='none')
        ax_source.set_yticks([])
        ax_source.set_xticks([])
        ax_source.set_title('Isolate\\nSource', pad=15, fontweight='bold', fontsize=10)
        
        for i, name in enumerate(unique_sources):
            legend_elements.append(Patch(facecolor=cmap_src(i), edgecolor='gray', label=f"Source: {name}"))

    # D. Draw Subspecies Metadata Strip Ribbon
    cmap_sub = plt.cm.get_cmap('Pastel1', len(unique_subsp))
    ax_subsp.matshow(subsp_ints, cmap=cmap_sub, aspect='auto', interpolation='none')
    ax_subsp.set_xticks([])
    ax_subsp.set_title('Subspecies\\n(Kraken)', pad=15, fontweight='bold', fontsize=10)

    ax_subsp.set_yticks(range(num_strains))
    ax_subsp.set_yticklabels(strains_ordered, fontsize=11, fontweight='bold')
    ax_subsp.yaxis.tick_right()

    for i, name in enumerate(unique_subsp):
        legend_elements.append(Patch(facecolor=cmap_sub(i), edgecolor='gray', label=f"Subsp: {name}"))

    fig.legend(handles=legend_elements, loc='lower center', ncol=min(6, len(legend_elements)), 
               bbox_to_anchor=(0.5, -0.05), fontsize=10, frameon=True)

    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('cloud_pangenome_annotated.png', dpi=300, bbox_inches='tight')
    plt.clf()
    """
}
