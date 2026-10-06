nextflow.enable.dsl=2

process prokka {
    // Automatically copies the prokka results directory back to your output folder
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'
    container 'docker://staphb/prokka:latest'
    input:
    // 1. Unpack the 5-item bundle flowing out of pyProc1
    //tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs),path(elgato_folder), path(fasta_file)
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs), path(elgato_folder), path(fasta_file) // same order as previous process
    output:
    // 2. Pass everything forward and add the generated prokka directory to the bundle
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs), path(elgato_folder), path(fasta_file), path("prokka_out")

    script:
    """
    # 3. Force Java to bypass shared /tmp performance locks
    #tells the Java virtual machine to completely skip creating or checking files inside /tmp/hsperfdata_*.
    export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

    # 4. Read variables directly from the local file without path splitting
    genus=\$(cat ${pyoutputs} | cut -d "," -f 1)
    species=\$(cat ${pyoutputs} | cut -d "," -f 2)

    # 5. Run prokka inside the local workspace, outputting to a clean local directory
    prokka \\
        --genus \${genus} \\
        --species \${species} \\
        --strain ${fasta_file.baseName} \\
        --outdir prokka_out \\
        --prefix ${fasta_file.baseName} \\
        --force \\
        --compliant \\
        --locustag \${genus} \\
        ${fasta_file}
    """
}
 
