process legsta {
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'
    container 'docker://staphb/legsta:latest'
    input:
    tuple path(elgato_folder), path(fasta_file)

    output:
    // Pass both the generated tsv AND the fasta file to the next process
    tuple path("${fasta_file.baseName}_leg.tsv"), path(elgato_folder), path(fasta_file)

    script:
    """
    #apptainer exec docker://staphb/legsta \\
    legsta ${fasta_file} > ${fasta_file.baseName}_leg.tsv
    """
}