#!/usr/bin/env nextflow

process MERGE_RECONSTRUCTED_QURE {
          label 'process_low'

    tag "$sample_id"
    container "cacciabue/multiquas:developing"

    input:
    tuple   val(sample_id), path(haplotypes), path(adjusted), val(contig_name) 
  
    output:
    tuple val("$sample_id"),
          path("${sample_id}_qure_haplotypes.fasta"), 
                path("${sample_id}_qure_proportions.txt"),
                 val("${contig_name}")

    script:
    """
 
    
    cat ${haplotypes} >  ${sample_id}_qure_haplotypes.fasta
    
    
    cat ${adjusted} >  ${sample_id}_qure_proportions.txt
    
   
    """
}