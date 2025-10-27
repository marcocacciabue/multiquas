#!/usr/bin/env nextflow

process BBDUK {
     tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/trimming", mode: 'copy'
    memory "4G"
    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
  
    tuple  val("${sample_id}"),path("${read1.simpleName}_trimmed.fq"), path("${read2.simpleName}_trimmed.fq"), emit: trimmed_reads
    


    script:

    """
   
    bbduk.sh in1=${read1} out1='${read1.simpleName}_trimmed.fq' in2=${read2} out2='${read2.simpleName}_trimmed.fq' ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl

    """
}