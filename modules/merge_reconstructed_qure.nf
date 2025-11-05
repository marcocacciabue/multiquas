#!/usr/bin/env nextflow

process MERGE_RECONSTRUCTED_QURE {
    tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/haplotypes_qure/", mode: 'copy'    

    input:
    tuple path(haplotypes), path(adjusted), val(sample_id) 
  
    output:
    tuple path("${sample_id}_qure_haplotypes.fasta"), 
        val("$sample_id"),
        path("${sample_id}_qure_proportions.txt"),
        path("${sample_id}_qure_haplotypes_aligned.fasta")
    
    cpus 4

    script:
    """
 
    
    cat ${haplotypes} >  ${sample_id}_qure_haplotypes.fasta
    
    mafft --thread ${task.cpus} --auto --quiet ${sample_id}_qure_haplotypes.fasta > ${sample_id}_qure_haplotypes_aligned.fasta
    
    cat ${adjusted} >  ${sample_id}_qure_proportions.txt
    
   
    """
}