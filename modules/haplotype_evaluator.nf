#!/usr/bin/env nextflow

process HAPLOTYPE_EVALUATOR {
  tag "$sample_id"
  container "cacciabue/multiquas:developing"
  publishDir "results/${sample_id}/haplotypes_qure", mode: 'copy'

  
  input:
    tuple val(sample_id), 
    path(haplotypes),
    path(proportions), 
    path(haplotypes_aligned), 
    path(variants_vcf)

  
  output:
    path("${sample_id}_qure_graphs.png")
  script:
    """
   Rscript /root/data/haplotype_evaluator.R ${haplotypes_aligned} ${proportions} ${variants_vcf} ${sample_id}_qure_graphs.png
    """
}