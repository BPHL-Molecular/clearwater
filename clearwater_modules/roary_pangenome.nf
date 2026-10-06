process PANGENOME_ANALYSIS {
    // Automatically copies the entire pangenome output folder back to your output master folder
    publishDir { "${params.output}/pangenome_analysis" }, mode: 'copy'
    container 'docker://staphb/roary:latest'
    
    cpus 4 // Adjust based on your infrastructure allocations

    input:
    // Receives a list of ALL individual gff files gathered from across your dataset
    path gff_files

    output:
    // Captures the complete results folder
    path "roary_out", emit: pgrun_dir

    script:
    """
    # Nextflow automatically links all 10 .gff files into the current folder.
    # We pass '*.gff' so Roary catches every single one of them.
    # -f: Output directory name
    # -p: CPU thread allocation
    # -e: Create a core gene alignment using MAFFT
    # -v: Verbose logging output
    roary -f roary_out -p ${task.cpus} -e --mafft ${gff_files}
    """
}