process pyProc7 {
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'
    container null

    input:
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_ani), path(elgato_folder), path(fasta_file), path(prokka_out)

    output:
    tuple path("report.txt"), path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_ani), path(elgato_folder), path(fasta_file), path(prokka_out)

    script:
    """
    #!/usr/bin/env python3

    sample_id = "${fasta_file.baseName}"

    results = []
    with open("${pyoutputs_ani}", "r") as aline:
        for line in aline:
            if not line.strip():
                continue
            cells = line.rstrip().split(",")

            while len(cells) < 22:
                cells.append("N/A")

            results = [
                sample_id,
                cells[0] + '_' + cells[1], 
                cells[3],                  
                cells[2],                  
                cells[14],                 
                cells[15],                 
                cells[16],                 
                cells[17],                 
                cells[4],                  
                cells[5],                  
                cells[18],                 # 🌟 Safely holds best_ref
                cells[19],                 # 🌟 Safely holds best_ani
                cells[6],                  
                cells[7],                  
                cells[8],                  
                cells[9],                  
                cells[10],                 
                cells[11],                 
                cells[12],                 
                cells[13]                  
            ]

    with open("report.txt", "w") as report:
        header = [
            'sampleID', 'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance',
            'Assembly_Completeness', 'Contamination_level', 'Genome_size', 'Ncontig',
            'subspeciesID_kraken', 'krakenSubsp_percent', 'Best_Ref_Genome', 'Highest_ANI_Percent',
            'SeqType', 'Gene_flaA', 'Gene_pilE', 'Gene_asd', 'Gene_mip', 'Gene_mompS','Gene_proA', 'Gene_neuA_neuAh'
        ]

        report.write("\\t".join(header) + "\\n")
        if results:
            report.write("\\t".join(results) + "\\n")
        else:
            empty_results = [sample_id] + ["N/A"] * (len(header) - 1)
            report.write("\\t".join(empty_results) + "\\n")
    """
}
