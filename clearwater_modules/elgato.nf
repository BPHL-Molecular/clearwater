//nextflow.enable.dsl=2

process elgato4mlst {
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'
    container 'docker://staphb/elgato:latest'
    input:
    path fasta_file // Receives the file object directly

    output:
    // Emits the generated folder AND the original fasta file safely
    tuple path("elgato_out"), path(fasta_file)

    script:
    """
    mkdir -p elgato_out
    #apptainer exec docker://staphb/elgato \\
    el_gato.py \\
        --assembly ${fasta_file} \\
        --out elgato_out \\
        --sample ${fasta_file.baseName} \\
        --overwrite
    """
}