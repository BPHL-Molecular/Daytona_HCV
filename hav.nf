#!/usr/bin/env nextflow

/*
Note:
Before running the script, please set the parameters in the config file params.yaml
*/

//Step1:input data files
nextflow.enable.dsl=2
def L001R1Lst = []
def sampleNames = []
myDir = file("$params.input")

myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
L001R1Lst.sort()
L001R1Lst.each{
   def x = it.minus("_1.fastq.gz")
     //println x
   sampleNames.add(x)
}
//println L001R1Lst
//println sampleNames


//Step2: process the inputed data
A = Channel.fromList(sampleNames)
//A.view()

include { quality } from './modules/quality.nf'
include { nofrag } from './modules/nofrag_hav.nf'
include { frag } from './modules/frag_hav.nf'
//include { primer } from './modules/primer.nf'
include { unprimer } from './modules/unprimer.nf'
include { assembly } from './modules/assembly_hav.nf'
include { pystats } from './modules/pystats.nf'

workflow {
    if("${params.frag}" == "frag"){
       quality(A) | frag | unprimer | assembly | pystats | view
    }
    else{
       quality(A) | nofrag | unprimer | assembly | pystats | view
    }
}


