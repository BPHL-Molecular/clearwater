#!/usr/bin/env nextflow

process pyProc5 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
    script:    
    """

    #!/usr/bin/env python3
    import subprocess
    import json
    import itertools
    import fileinput
    import os

    #items = "${mypath}".strip().split("/")
    #filepath1="${mypath}/"+items[-1]+"_assembly/prokka/"+items[-1]+".txt"
    #filepath2="${mypath}/"+items[-1]+"_assembly/"+items[-1]+".mlst"

    filepath3="${mypath}/possible_mlsts.txt"
    filepath4="${mypath}/report.txt"
    #uppath="/".join(items[:-2])
    
    f = open("${pyoutputs}", "a") # file to write report results
    try:
        with open(filepath3, "r") as mlstreport:
            next(mlstreport, None)# Skip the first line
            second_line = next(mlstreport, None) # Read the second line(allele number)
            if second_line is not None:
                var1 = second_line.rstrip().split("\t") # convert to a list and remove \n
                #var2 = tuple(var1[2:]) # omit first item
                #var3 = ",".join(var2) # remove space
                #var4 = f"({var3})"
                seqtype = var1[1]
                flA = var1[2]
                pilE = var1[3]
                asd = var1[4]
                mip = var1[5]
                mompS = var1[6]
                proA = var1[7]
                neuA_neuAH = var1[8] 
                f.write(","+str(seqtype)+","+str(flA)+","+str(pilE)+","+str(asd)+","+str(mip)+","+str(mompS)+","+str(proA)+","+str(neuA_neuAH))
            else:
                var3 = "Gene absent"
                f.write("," +str(var3))
                
    except FileNotFoundError:
                         var3 ="check sample"
                         f.write(","+str(var3))
                         
    except Exception as e:
                         var3 = "error occurred {e}"
                         f.write(","+str(var3))
                         

    f.close
    
    
    """
}