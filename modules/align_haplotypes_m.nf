#!/usr/bin/env nextflow

process ALIGN_HAPLOTYPES_M {
   label 'process_low'
  tag "$sample_id"
  container "cacciabue/multiquas:developing"
  publishDir "results/${sample_id}/haplotypes_qure/evaluator_multiple", mode: 'copy'

  
  input:
    tuple val(sample_id), 
    path(haplotypes),
    path(proportions), 
     val(contig_name),
    path(variants_vcf),
    path(first_ref)
    
  
  output:
     tuple val("${sample_id}"),
     path("${sample_id}_qure_haplotypes_aligned.fasta"),
     path("${proportions}"),
     path("${variants_vcf}"),
     val("{$contig_name}")
  
  script:
    """
    
    cat    ${first_ref} ${haplotypes}  > haplo_plus_ref.fasta    
    mafft --thread ${task.cpus} --auto --quiet haplo_plus_ref.fasta > ${sample_id}_qure_haplotypes_aligned.fasta

    #Rscript /root/data/haplotype_evaluator.R ${sample_id}_qure_haplotypes_aligned.fasta ${proportions} ${variants_vcf} ${sample_id}_qure_graphs.png
    """
}
