nextflow.enable.dsl=2

process kraken_subs {
    input:
        val mypath
        //path pyoutputs
    output:
        //stdout
        val mypath
        //path pyoutputs
    // backslash (`\`) to prevent Nextflow from trying to interpret them as Nextflow variables    
    """
    samplename=\$(echo "${mypath}" | rev | cut -d "/" -f 1 | rev)
    #mkdir -p ${params.output}/\${samplename}
    #mkdir -p ${params.output}/\${samplename}/kraken_out
    mkdir ${mypath}/kraken_out
    kraken2 --db /blue/bphl-florida/t.bazile1/legionella_pn/legion_pipeline/krakendb_subs/leg_subsp_customedDB/ --threads $task.cpus --use-names --report ${mypath}/kraken_out/\${samplename}.report --output ${mypath}/kraken_out/\${samplename}_kraken.out ${params.input}/\${samplename}.fasta
        
    """
}