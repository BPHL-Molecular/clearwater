process checkm_qc {
    publishDir "${params.output}/checkm_qc", mode: 'copy'
    //cpus 8
    //memory '32 GB'
    container 'docker://nanozoo/checkm'  
    input:
    // Accepts the collected list of all fasta files at once
    path all_fastas

    output:
    path "checkm_results", emit: results_dir
    path "checkm_summary.txt", emit: summary

    script:
    """
    # 1. Organize fasta files into a true physical folder to bypass symlinks
    mkdir -p raw_assemblies
    cp -L --remove-destination *.fasta raw_assemblies/
    
    # 2. Run CheckM Lineage Workflow inside the Apptainer container wrapper
    # --no-home prevents cluster environment leaks
    #apptainer exec docker://nanozoo/checkm \\
        checkm lineage_wf \\
        --threads ${task.cpus} \\
        --extension fasta \\
        raw_assemblies/ \\
        checkm_results/
    # 3. FIX: CheckM2 natively creates a clean table called quality_report.tsv.
    # We copy it to match your expected file name so your sbatch script works unchanged.
    #apptainer exec docker://nanozoo/checkm \\
        checkm qa \\
        checkm_results/lineage.ms \\
        checkm_results/ \\
        --tab_table \\
        --file checkm_summary.txt #outputs a clean, tab-separated table containing precise percentages for Completeness and Contamination
    """
}