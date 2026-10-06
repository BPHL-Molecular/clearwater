process PLOT_PANGENOME {
    tag { "Generating Pangenome Plots" }
    publishDir "${params.output}/pangenome_plots", mode: 'copy'
    
    container null

    input:
    path roary_outdir

    output:
    path "*.png", emit: pangenome_plots

    script:
    """
#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
from Bio import Phylo

sns.set_style('white')

# Define inputs based on Nextflow directory mapping
tree_file = "${roary_outdir}/accessory_binary_genes.fa.newick"
spreadsheet_file = "${roary_outdir}/gene_presence_absence.csv"

# Read tree and calculate distance limits
t = Phylo.read(tree_file, 'newick')
mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

# Load and clean roary spreadsheet data matrix
roary = pd.read_csv(spreadsheet_file, low_memory=False)
roary.set_index('Gene', inplace=True)
roary.drop(list(roary.columns[:13]), axis=1, inplace=True)

# Convert values into binary format
roary.replace('.{2,100}', 1, regex=True, inplace=True)
roary.replace(np.nan, 0, regex=True, inplace=True)

# Calculate frequencies before sorting columns
num_strains = roary.shape[1]
gene_counts = roary.sum(axis=1)

# Sort the matrix index rows by strain frequencies
idx = gene_counts.sort_values(ascending=False).index
roary_sorted = roary.loc[idx]

# 1. Pangenome Frequency Histogram Plot
plt.figure(figsize=(7, 5))
plt.hist(gene_counts, num_strains, histtype="stepfilled", alpha=.7)
plt.xlabel('No. of genomes')
plt.ylabel('No. of gene clusters')
sns.despine(left=True, bottom=True)
plt.savefig('pangenome_frequency.png', dpi=300)
plt.clf()

# Align matrix data explicitly with terminal tips of the tree
roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]

# 2. Complete Core/Accessory Matrix Heatmap Plot
with sns.axes_style('whitegrid'):
    fig = plt.figure(figsize=(18, 12))
    
    ax1 = plt.subplot2grid((1,40), (0, 10), colspan=30)
    ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues, vmin=0, vmax=1, aspect='auto', interpolation='none')
    ax1.set_yticks([])
    ax1.set_xticks([])

    # Calculate exact column boundary positions based on frequencies
    sorted_counts = gene_counts.loc[idx].values
    total_genes = len(sorted_counts)
    
    core_len = len(np.where((sorted_counts >= num_strains * 0.99) & (sorted_counts <= num_strains))[0])
    softcore_len = len(np.where((sorted_counts >= num_strains * 0.95) & (sorted_counts < num_strains * 0.99))[0])
    shell_len = len(np.where((sorted_counts >= num_strains * 0.15) & (sorted_counts < num_strains * 0.95))[0])
    cloud_len = len(np.where(sorted_counts < num_strains * 0.15)[0])

    x0 = 0.0
    x1 = core_len / total_genes
    x2 = (core_len + softcore_len) / total_genes
    x3 = (core_len + softcore_len + shell_len) / total_genes
    x4 = 1.0

    # Draw continuous baseline
    ax1.plot([0, total_genes], [num_strains, num_strains], color='black', lw=1.5, clip_on=False)

    # Draw vertical ticks at boundaries
    for boundary in [x0, x1, x2, x3, x4]:
        ax1.plot([boundary, boundary], [-0.02, 0.0], color='gray', lw=1.5, clip_on=False, transform=ax1.get_xaxis_transform())

    # Define staggered configurations to eliminate text overlaps
    sections = [
        (x0, x1, f'Core\\n({core_len})', '#0d47a1', -0.05),
        (x1, x2, f'Soft-core\\n({softcore_len})', '#1976d2', -0.13),
        (x2, x3, f'Shell\\n({shell_len})', '#e65100', -0.05),
        (x3, x4, f'Cloud\\n({cloud_len})', '#1b5e20', -0.05)
    ]

    for start, end, label, color, y_offset in sections:
        if start == end: 
            continue
        ax1.plot([start, end], [-0.01, -0.01], color=color, lw=4, clip_on=False, transform=ax1.get_xaxis_transform())
        ax1.text((start + end) / 2, y_offset, label, transform=ax1.transAxes, ha='center', va='top', fontsize=11, color=color, fontweight='bold')

    try:
        ax = plt.subplot2grid((1,40), (0, 0), colspan=10, facecolor='white')
    except AttributeError:
        ax = plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

    fig.subplots_adjust(wspace=0, hspace=0)
    ax1.set_title('Roary matrix\\n(%d gene clusters)' % total_genes, pad=15)

    Phylo.draw(t, axes=ax, show_confidence=False, label_func=lambda x: None,
               xticks=([],), yticks=([],), ylabel=('',), xlabel=('',),
               xlim=(-mdist*0.1, mdist+mdist*0.1), axis=('off',),
               title=('Tree\\n(%d strains)' % num_strains,), do_show=False)
    
    plt.savefig('pangenome_matrix.png', dpi=300, bbox_inches='tight')
    plt.clf()

# 3. Pangenome Distribution Partition Pie Chart Plot
plt.figure(figsize=(10, 10))
core = core_len
softcore = softcore_len
shell = shell_len
cloud = cloud_len

def my_autopct(pct):
    val = int(round(pct*total_genes/100.0))
    return '{v:d}'.format(v=val)

plt.pie([core, softcore, shell, cloud],
        labels=['core\\n(%d <= strains <= %d)' % (num_strains*.99, num_strains),
                'soft-core\\n(%d <= strains < %d)' % (num_strains*.95, num_strains*.99),
                'shell\\n(%d <= strains < %d)' % (num_strains*.15, num_strains*.95),
                'cloud\\n(strains < %d)' % (num_strains*.15)],
        explode=[0.1, 0.05, 0.02, 0], radius=0.9,
        colors=[(0, 0, 1, float(x)/total_genes) for x in (core, softcore, shell, cloud)],
        autopct=my_autopct,
        textprops={'fontsize': 11})

plt.savefig('pangenome_pie.png', dpi=300, bbox_inches='tight')
plt.clf()
    """
}
