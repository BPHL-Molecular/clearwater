nextflow.enable.dsl=2

//legsta sanibel_dir/output_all_20250421-0501/${samp}/${samp}_assembly/${samp}.fasta > legsta_out/${samp}_leg.tsv

process legsta {
   input:
      val mypath
   output:
      //stdout
      val "${mypath}" //val "${workflow.workDir}/", emit: outputpath1
      //${params.output}/${x}/{${x}_leg.tsv
      //${params.input}/*fasta
      // mkdir -p ${params.output}/${x}
      //mv ${params.output}/${x}_leg.tsv ${params.output}/${x}/
   """
   samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
   
   legsta ${params.input}/\${samplename}.fasta > ${mypath}/\${samplename}_leg.tsv 
   
   """
}

//apptainer exec docker://staphb/ispcr
//apptainer exec docker://staphb/any2fasta
//apptainer exec docker://staphb/legsta