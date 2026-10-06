
nextflow.enable.dsl=2

process kraken_subs {
    publishDir { "${params.output}/${fasta_file.baseName}/kraken_out" }, mode: 'copy'
    container 'docker://staphb/kraken2:latest'
    input:
    // Receives the database directory file object
    path(krakenSub_db, stageAs: 'kraken_db/*')
    
    // Receives the tsv from legsta and the original fasta file
    tuple path(leg_tsv), path(elgato_folder), path(fasta_file)

    output:
    // Outputs a tuple containing the report, the raw output, and the fasta file for the next steps
    tuple path("${fasta_file.baseName}.report"), path("${fasta_file.baseName}_kraken.out"), path(elgato_folder), path(fasta_file)
    // ${params.krakenSub_db} to access in script the krakenSub_db without passing though the input section
    script:
    """
    # Run kraken2 purely on local files staged by Nextflow
    kraken2 \\
        --db ${krakenSub_db} \\
        --threads ${task.cpus ?: 1} \\
        --use-names \\
        --report ${fasta_file.baseName}.report \\
        --output ${fasta_file.baseName}_kraken.out \\
        ${fasta_file}
    """
}