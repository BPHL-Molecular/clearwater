process pyProc3 {
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'
    container null

    input:
    // FIXED: Combined everything into a single input tuple matching the channel structure
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs), path(elgato_folder), path(fasta_file), path(prokka_out), path(checkm_results_dir)

    output:
    // Pass forward the standard sample tuple (dropping the checkm folder since pyProc7 doesn't need it)
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path("pyoutputs_updated.txt"), path(elgato_folder), path(fasta_file), path(prokka_out)

    script:
    """
    #!/usr/bin/env python3
    import ast
    import sys
    import os

    sample_id = "${fasta_file.baseName}"
    stats_file_path = os.path.join("${checkm_results_dir}", "storage", "bin_stats_ext.tsv")

    completeness  = "N/A"
    contamination = "N/A"
    genome_size   = "N/A"
    n_contigs     = "N/A"

    if os.path.exists(stats_file_path):
        with open(stats_file_path, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split(None, 1)
                if len(parts) < 2:
                    continue
                
                if parts[0] == sample_id:
                    try:
                        stats = ast.literal_eval(parts[1])
                        completeness  = str(stats.get('Completeness', 'N/A'))
                        contamination = str(stats.get('Contamination', 'N/A'))
                        genome_size   = str(stats.get('Genome size', 'N/A'))
                        n_contigs     = str(stats.get('# contigs', 'N/A'))
                        break
                    except Exception as e:
                        print(f"Error parsing dictionary data for {sample_id}: {e}", file=sys.stderr)

    existing_data = ""
    if os.path.exists("${pyoutputs}"):
        with open("${pyoutputs}", "r") as f:
            existing_data = f.read().strip()

    updated_row = f"{existing_data},{completeness},{contamination},{genome_size},{n_contigs}"

    with open("pyoutputs_updated.txt", "w") as out:
        out.write(updated_row + "\\n")
    """
}

