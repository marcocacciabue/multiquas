#!/usr/bin/env nextflow

process MERGE_RECONSTRUCTED {
  label 'process_low'

  container "cacciabue/multiquas:tools.v0.0.1"
  tag "${[sample_id, contig_name].join(':')}"
  containerOptions "--cpus=${task.cpus}"
  input:
    tuple val(sample_id),
          path(first_ref),
          path(general_ref),
          val(reconstructer),
          val(method),
          val(contig_name),
          path(haplotypes),
          path(adjusted)
  
  output:
    tuple val("$sample_id"),
          path("${first_ref.first()}"),
          path("${general_ref.first()}"),
          val("${reconstructer.first()}"),
          val("${method.first()}"),
          val("${contig_name.first()}"),
          path("${sample_id}_${reconstructer.first()}_${method.first()}_haplotypes.fasta"), 
          path("${sample_id}_${reconstructer.first()}_${method.first()}_proportions.txt")
  
  script:
 

    """
 
    
    cat ${haplotypes} >  '${sample_id}_${reconstructer.first()}_${method.first()}_haplotypes.fasta'
    
    
    cat ${adjusted} >  '${sample_id}_${reconstructer.first()}_${method.first()}_proportions.txt'
    
   
    """
}
