//nextflow.enable.dsl=2

process elgato4mlst {
   input:
      val x
   output:
      //stdout
      //path "${params.output}/elgato_out/${x.baseName}" //"${params.output}/${x}" //"${workflow.workDir}", emit: outputpath1
      val "${params.output}/${x}", emit: outputpath1
   script:
   """
   mkdir -p ${params.output}/${x}
   
   #cp ${params.input}/${x}.fasta ${params.output}/${x}/      
   #apptainer exec docker://staphb/elgato
   cp ${params.input}/${x}.fasta ${params.output}/${x}/ 
   el_gato.py --assembly ${params.input}/${x}.fasta --out ${params.output}/${x} --sample ${x} --overwrite
  
   """
}