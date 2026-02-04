#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.input=""
//params.output ="/blue/bphl-florida/t.bazile1/legionella_pn/legion_pipeline/out"

myDir = file("$params.input")
def sampLst =[]
def sampNames =[]

myDir.eachFileMatch ~/.*.fasta/, {sampLst << it.name}
sampLst.sort()

sampLst.each{
   def x = it.minus(".fasta")
   sampNames.add(x)
   }
   
A = Channel.fromList(sampNames)
A.view()

//C = Channel.fromPath("${params.input}/*", type: 'file')
                          //.collect()
//C.view()
//D = Channel.fromPath("${params.input}")
//D.view()
//C = Channel.of(Clist)

//B= Channel.value($params.input/"*.fasta")
include { elgato4mlst } from './legCabbage_modules/elgato.nf'
include { legsta } from './legCabbage_modules/legsta.nf'
include { kraken_subs } from './legCabbage_modules/kraken_s.nf'
include {proc_mash} from './legCabbage_modules/leg_mash.nf'
include {pyProc1} from './legCabbage_modules/pyproc1.nf'
include {prokka} from './legCabbage_modules/prokka.nf'
include {pyProc4} from './legCabbage_modules/pyproc4.nf'
include {pyProc5} from './legCabbage_modules/pyproc5.nf'
include {pyProc7} from './legCabbage_modules/pyproc7.nf' //final report
//include { proc_gff } from './legCabbage_modules/proc_gff.nf'

workflow {
    //A |(legsta & elgato4mlst & kraken_subs)
    elgato4mlst(A)|legsta|kraken_subs|proc_mash|pyProc1|prokka|pyProc4|pyProc5|pyProc7 //everything in one package dir 
    //elgato4mlst(A) 
    //A|elgato4mlst|view
    //legsta(myDir)|view
    //proc_gff(myDir) | view
    //elgato4mlst.out.view()
}

// feeding two process with same input
// A |(legsta&elgato4mlst)
