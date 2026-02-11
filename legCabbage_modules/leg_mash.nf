nextflow.enable.dsl=2

//legsta sanibel_dir/output_all_20250421-0501/${samp}/${samp}_assembly/${samp}.fasta > legsta_out/${samp}_leg.tsv

process proc_mash {
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
   mkdir -p ${mypath}/mash_out
   #cp ${params.input}/\${samplename}.fasta  ${mypath}/mash_out/\${samplename}.fasta 
   mash sketch -o ${mypath}/mash_out/\${samplename}_sketch ${params.input}/\${samplename}.fasta
   mash dist /db/RefSeqSketchesDefaults.msh ${mypath}/mash_out/\${samplename}_sketch.msh > ${mypath}/mash_out/\${samplename}_distances.tab
   sort -gk3 ${mypath}/mash_out/\${samplename}_distances.tab | head > ${mypath}/mash_out/\${samplename}_distances_top10.tab   
   """
}

//apptainer exec docker://staphb/ispcr
//apptainer exec docker://staphb/any2fasta
//apptainer exec docker://staphb/legsta