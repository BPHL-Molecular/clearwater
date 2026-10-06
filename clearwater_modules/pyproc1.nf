process pyProc1 {
    // Automatically copies the partial report to the sample directory
    publishDir { "${params.output}/${fasta_file.baseName}" }, mode: 'copy'
    // Force Nextflow NOT to run this process in a container
    container null
    
    input:
    // Unpack the bundle coming from proc_mash
    tuple path(kraken_report), path(kraken_out),path(elgato_folder), path(mash_top10), path(fasta_file)

    output:
    // Pass everything forward AND add the new pyoutputs.txt report to the bundle
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path("pyoutputs.txt"), path(elgato_folder), path(fasta_file)

    script:
    """
    #!/usr/bin/env python3
    import sys
    import os
    import re

    # 1. Use the Nextflow local variable directly
    filepath1 = "${mash_top10}"

    with open(filepath1, 'r') as mash:
        top_hit = mash.readline()
        if top_hit:
            top_hit = str(top_hit)
            #gn = re.sub('.*-\\.-', '', top_hit)
            gn = re.sub('.*-\\\\.-', '', top_hit)
            cells = gn.split()
            acell = cells[0]
            acell = acell.lstrip("_")
            #agn = re.split('^([^_]*_[^_]*)(_|\\.).*\$', acell)[1]
            agn = re.split('^([^_]*_[^_]*)(_|\\\\.).*\$', acell)[1]
            genus = agn.split('_')[0]
            species = agn.split('_')[1]
            distance = top_hit.split()[2]
            accession = top_hit.split("-")[5]
            
            # 2. Write straight to the local workspace file
            with open("pyoutputs.txt", "w") as f:
                f.write(f"{genus},{species},{distance},{accession}")
    """
}
