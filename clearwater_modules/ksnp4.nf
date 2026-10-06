process ksnp4 {
    cpus 4
    memory { task.attempt == 1 ? '16 GB' : (32.GB * (task.attempt - 1)) }
    time '4h'
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 2

    container 'docker://staphb/ksnp4:latest'
    containerOptions { '--env PATH="/usr/local/bin:/ksnp4/bin:$PATH"' }

    publishDir "${params.output}/ksnp_out", mode: 'copy'

    input:
    path all_fastas 

    output:
    path "run_lpnSamref2", emit: run_dir
    path "Kchooser_lpnSamref2_input.report"
    path "lpnSamref2_input.in"
    path "run_lpnSamref2Log"

    script:
    def cpu_count = task.cpus // Keep this as Groovy since Nextflow manages cpus
    """
    # 1. Setup environment
    export HOME=\$(pwd)  
    mkdir -p local_scratch
    export TEMP=\$(pwd)/local_scratch
    export TMPDIR=\$(pwd)/local_scratch

    # 2. Isolate and cleanly copy inputs
    mkdir -p cleaned_genomes
    for f in ${all_fastas}; do
        BASE_NAME=\$(basename "\$f" | sed -r 's/\\.(fasta|fa|fna)\$//')
        CLEAN_NAME=\$(echo "\${BASE_NAME}" | sed 's/[^a-zA-Z0-9_-]/_/g')
        cp -L "\$f" "cleaned_genomes/\${CLEAN_NAME}.fasta"
    done

    # 3. Capture absolute path
    ABS_CLEANED_DIR=\$(pwd)/cleaned_genomes

    # 4. Generate the 2-column input file
    MakeKSNP4infile -indir "\${ABS_CLEANED_DIR}" -outfile lpnSamref2_input.in

    # 5. Run Kchooser natively
    Kchooser4 -in lpnSamref2_input.in > Kchooser_lpnSamref2_input.report

    # --- NEW BASH STEP ---
    # Extract the K value using a Bash grep/sed command.
    # (Adjust the regex 'The optimum value of K is' to match your actual file output!)
    OPTIMUM_K_VALUE=\$(grep "The optimum value of K is" Kchooser_lpnSamref2_input.report | awk '{print \$NF}')
    
    # Fallback to 17 if Bash failed to parse a number
    if [ -z "\$OPTIMUM_K_VALUE" ]; then
        OPTIMUM_K_VALUE=17
    fi
    # ---------------------

    # 6. Run kSNP4 using the Bash variable \$BASH_K_VALUE
    kSNP4 \\
        -in lpnSamref2_input.in \\
        -outdir run_lpnSamref2 \\
        -k "\$OPTIMUM_K_VALUE" \\
        -CPU ${cpu_count} \\
        -ML \\
        -core \\
        -nocache \\
        2>&1 | tee run_lpnSamref2Log
    """
}
