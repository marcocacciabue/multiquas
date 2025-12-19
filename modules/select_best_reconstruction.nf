#!/usr/bin/env nextflow

process SELECT_BEST_RECONSTRUCTION {
  label 'process_low'
  tag "$sample_id"
  container "cacciabue/multiquas:developing"
  publishDir "results/${sample_id}/best", mode: 'copy'
  
  
  input:
  tuple val(sample_id), 
  path(haplotypes),
  path(proportions),
  path(graphs),
  val(r_sq)
  
  output:
    tuple val("$sample_id"), 
  path("$haplotypes"),
  path("$proportions"),
  path("$graphs"),
  val("$r_sq")
  
  
  script:
  """
  echo $r_sq

  """
}
