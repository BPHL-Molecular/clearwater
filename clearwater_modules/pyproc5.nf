#!/usr/bin/env nextflow

process pyProc5 {
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'

    // Force Nextflow NOT to run this process in a container
    container null
    
    input:
    // 1. Unpack the 7-item bundle flowing out of pyProc4
    //tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs), path(prokka_out), path(elgato_folder), path(fasta_file)
    // Clear variable naming to match exactly what came from pyProc4
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_updated2), path(elgato_folder), path(fasta_file), path(prokka_out)
    output:
    //Forward the 7-item bundle to pyProc7, output dir should be quoted ""
    //tuple path(kraken_report), path(kraken_out), path(mash_top10), p0ath(pyoutputs), path(prokka_out), path(elgato_folder), path(fasta_file)
    // Generate the final cumulative report file
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path("pyoutputs_final.txt"), path(elgato_folder), path(fasta_file), path(prokka_out)

    script:
    """
    #!/usr/bin/env python3
    import os

    # 3. Because Nextflow stages the folder locally using its name in the tuple,
    # we can point straight to it in the execution path:
    filepath3 = os.path.join("${elgato_folder}", "possible_mlsts.txt")
    
    # Use the matching variable name here
    input_file = "${pyoutputs_updated2}"
    output_file = "pyoutputs_final.txt"

    # Copy records forward
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write(infile.read())
 

    # 4. Use python's 'with' statement to cleanly read and append text
    with open(output_file, "a") as f:
        try:
            with open(filepath3, "r") as mlstreport:
                next(mlstreport, None) # Skip the first line
                second_line = next(mlstreport, None) # Read the second line
                if second_line is not None:
                    var1 = second_line.rstrip().split("\t") 
                    seqtype = var1[1]
                    flA = var1[2]
                    pilE = var1[3]
                    asd = var1[4]
                    mip = var1[5]
                    mompS = var1[6]
                    proA = var1[7]
                    neuA_neuAh = var1[8]
                    f.write(","+str(seqtype)+","+str(flA)+","+str(pilE)+","+str(asd)+","+str(mip)+","+str(mompS)+","+str(proA)+","+str(neuA_neuAh))
                else:
                    var3 = "Gene absent"
                    f.write("," + str(var3))

        except FileNotFoundError:
            var3 = "check sample"
            f.write("," + str(var3))

        except Exception as e:
            var3 = f"error occurred {e}"
            f.write("," + str(var3))
    """
}
