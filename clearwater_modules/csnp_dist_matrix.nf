process calculate_snp_matrix {
    container null
    tag { "Calculating Pairwise Core-SNP Distance Matrix" }
    publishDir "${params.output}/transmission_analysis", mode: 'copy'
    
    //container null // Runs natively inside your active bio-python environment

    input:
    path run_lpn_dir // Ingests your "run_lpnSamref2" directory safely

    output:
    path "pairwise_snp_distance_matrix.tsv", emit: matrix_tsv

    script:
    """
    #!/usr/bin/env python3
    from Bio import AlignIO
    import numpy as np
    import pandas as pd

    # 1. Access the sequence matrix inside the kSNP4 run directory output
    fasta_file = "${run_lpn_dir}/core_SNPs_matrix.fasta"
    alignment = AlignIO.read(fasta_file, "fasta")
    
    sample_names = [record.id for record in alignment]
    num_samples = len(sample_names)
    
    # Initialize an empty square distance matrix array
    dist_matrix = np.zeros((num_samples, num_samples), dtype=int)

    # 2. Iterate pairwise through string sequences to count mismatches
    for i in range(num_samples):
        seq_i = str(alignment[i].seq).upper()
        for j in range(i + 1, num_samples):
            seq_j = str(alignment[j].seq).upper()
            
            # Count physical mutations, ignoring gaps or missing ambiguous codes (N)
            mismatches = sum(1 for c1, c2 in zip(seq_i, seq_j) if c1 != c2 and c1 != 'N' and c2 != 'N' and c1 != '-' and c2 != '-')
            
            dist_matrix[i, j] = mismatches
            dist_matrix[j, i] = mismatches

    # 3. Export as a structured, clean tabular document
    matrix_df = pd.DataFrame(dist_matrix, index=sample_names, columns=sample_names)
    matrix_df.to_csv("pairwise_snp_distance_matrix.tsv", sep="\\t")
    """
}
