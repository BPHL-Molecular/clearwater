process runFastANI {
    tag "${fasta_file.baseName}"
    container 'docker://staphb/fastani:latest'

    // these completely ignore nextflow.config selectors and force Slurm's hand!
    cpus 4
    memory { 16.GB * task.attempt }
    time '2h'
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    
    input:
    // 1. The sample tuple (7 items coming out of pyProc3)
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_final), path(elgato_folder), path(fasta_file), path(prokka_out)
    // 2. The standalone reference directory channel
    path ref_dir 

    output:
    // Pass the 7 items + the generated fastani_raw.out forward (8 items total for pyProc6)
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_final), path("fastani_raw.out"), path(elgato_folder), path(fasta_file), path(prokka_out)
    // was under script block: def cpus = task.cpus ?: 2 and above """
    script:
    
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Added -L flag so find follows Nextflow's symbolic links
    find -L "${ref_dir}" -type f \\( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \\) > references.txt

    touch fastani_raw.out

    if [ -s references.txt ]; then
        fastANI \\
            --query "${fasta_file}" \\
            --rl references.txt \\
            --output fastani_raw.out \\
            --threads ${task.cpus} \\
            --fragLen 500 \\
            --minFraction 0.05
    fi
    """
}
