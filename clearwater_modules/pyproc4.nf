process pyProc4 {
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'

    // Force Nextflow NOT to run this process in a container
    container null
    
    input:
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_updated), path(elgato_folder), path(fasta_file), path(prokka_out)

    output:
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path("pyoutputs_updated2.txt"), path(elgato_folder), path(fasta_file), path(prokka_out)

    script:
    """
    #!/usr/bin/env python3
    import os

    filepath3 = "${kraken_report}"
    #input_file = "${params.output}/${fasta_file.baseName}/${pyoutputs_updated}"
    #CORRECT: Read directly from the dynamic Nextflow input file variable staged locally
    input_file = "${pyoutputs_updated}"
    output_file = "pyoutputs_updated2.txt"

    # Copy the existing contents to our new output file first
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write(infile.read())

    # Open the new file in append mode to add the kraken results
    with open(output_file, "a") as f:
        try:
            found_s1 = False
            with open(filepath3, 'r') as ksreport:
                for line in ksreport:
                    fields = line.strip().split('\t')
                    if len(fields) < 6:
                        continue

                    tx_level = fields[3].strip()
                    if tx_level == 'S1':
                        percentage = float(fields[0].strip())
                        taxonomic_name = fields[5].strip()
                        f.write("," + str(taxonomic_name) + "," + str(percentage))
                        found_s1 = True
                        break

            if not found_s1:
                percentage = "No Percentage"
                taxonomic_name = "No subspecies found"
                f.write("," + str(taxonomic_name) + "," + str(percentage))

        except FileNotFoundError:
            percentage = 'check sample'
            taxonomic_name = 'no report file'
            f.write("," + str(taxonomic_name) + "," + str(percentage))
        except Exception as e:
            percentage = f'error occurred {e}'
            taxonomic_name = 'error occurred'
            f.write("," + str(taxonomic_name) + "," + str(percentage))
    """
}
