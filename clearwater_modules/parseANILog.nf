process pyProc6 {
    container null
    tag "${fasta_file.baseName}"

    input:
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path(pyoutputs_final), path(fastani_raw), path(elgato_folder), path(fasta_file), path(prokka_out)

    output:
    tuple path(kraken_report), path(kraken_out), path(mash_top10), path("pyoutputs_ani.txt"), path(elgato_folder), path(fasta_file), path(prokka_out)

    script:
    """
#!/usr/bin/env python3
import os
import sys

best_ref = "N/A"
best_ani = "N/A"

raw_log = "${fastani_raw}"

if os.path.exists(raw_log) and os.path.getsize(raw_log) > 0:
    highest_score = 0.0
    with open(raw_log, "r") as fani:
        for line in fani:
            parts = line.strip().split()
            if len(parts) >= 3:
                current_score = float(parts[2])
                if current_score > highest_score:
                    highest_score = current_score
                    base_filename = os.path.basename(parts[1])
                    best_ref = base_filename.rsplit('.', 1)[0]
    if highest_score > 0.0:
        best_ani = f"{highest_score}%"

with open("${pyoutputs_final}", "r") as infile:
    cells = infile.read().strip()

with open("pyoutputs_ani.txt", "w") as out:
    out.write(cells + "," + best_ref + "," + best_ani + "\\n")
    """
}

    
