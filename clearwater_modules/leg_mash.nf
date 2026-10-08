nextflow.enable.dsl=2

process proc_mash {
    // Automatically copies your final mash results into the sample folder
    publishDir { "${params.output}/${fasta_file.baseName}/mash_out" }, mode: 'copy'
    container 'docker://staphb/mash:2.3'
    input:
    // Unpack everything sent down from the kraken_subs process
    tuple path(kraken_report), path(kraken_out),path(elgato_folder), path(fasta_file)

    output:
    // Bundle ALL files generated so far so downstream scripts have everything
    tuple path(kraken_report), path(kraken_out), path(elgato_folder), path("${fasta_file.baseName}_distances_top10.tab"), path(fasta_file)

    script:
    """
    # 1. Generate the sketch locally (automatically names it with .msh extension)
    mash sketch -o ${fasta_file.baseName}_sketch ${fasta_file}

    # 2. Calculate distances using the local sketch
    mash dist /db/RefSeqSketchesDefaults.msh ${fasta_file.baseName}_sketch.msh > ${fasta_file.baseName}_distances.tab

    # 3. Sort and isolate the top 10 matches
    sort -gk3 ${fasta_file.baseName}_distances.tab | head > ${fasta_file.baseName}_distances_top10.tab
    """
}

//apptainer exec docker://staphb/ispcr
//apptainer exec docker://staphb/any2fasta
//apptainer exec docker://staphb/legsta