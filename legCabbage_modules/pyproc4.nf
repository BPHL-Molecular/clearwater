
process pyProc4 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    #!/usr/bin/env python3

    import subprocess
    import json
    import itertools
    import fileinput
    import os

    items = "${mypath}".strip().split("/")
    #filepath1="${mypath}/"+items[-1]+"_assembly/prokka/"+items[-1]+".txt"
    #filepath2="${mypath}/"+items[-1]+"_assembly/"+items[-1]+".mlst"
    filepath3="${mypath}/kraken_out/"+items[-1]+".report" #// path to subspecies report
    filepath4="${mypath}/report.txt"
    #uppath="/".join(items[:-2])
    
    f = open("${pyoutputs}", "a") # file to write report results
    try:
        found_s1 = False
        with open(filepath3, 'r') as ksreport:
            for line in ksreport:
                # Kraken report lines are tab-delimited
                fields = line.strip().split('\t')
                
                if len(fields) < 6:
                    continue # Skip empty or incomplete lines
                
                # Tx Rank code is the 4th field (index 3)
                tx_level = fields[3].strip()
                
                if tx_level == 'S1':
                    percentage = float(fields[0].strip())
                    taxonomic_name = fields[5].strip()
                    f.write(","+str(taxonomic_name)+","+str(percentage))                      
                    found_s1 = True
                    break # Stop after the first occurrence
                
        if not found_s1:
            percentage = "No Percentage"
            taxonomic_name = "No subspecies found"
            f.write(","+str(taxonomic_name)+","+str(percentage))
    
    except FileNotFoundError:
        
                             percentage ='check sample',
                             taxonomic_name = 'no report file'
                             f.write(","+str(taxonomic_name)+","+str(percentage))
       
    except Exception as e:
                             percentage = 'error occurred {e}',
                             taxonomic_name = 'error occurred'
                             f.write(","+str(taxonomic_name)+","+str(percentage))


    f.close
        
    """
}