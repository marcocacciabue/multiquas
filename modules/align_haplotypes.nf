#!/usr/bin/env nextflow

process ALIGN_HAPLOTYPES {
  label 'process_low'

  container "cacciabue/multiquas:tools.v0.0.1"
    tag "${sample_id}"

  
  input:
    tuple val(sample_id),
          path(first_ref),
          path(general_ref), 
          val(reconstructer),
          val(method),
          val(contig_name),
          path(haplotypes),
          path(proportions),
          path(variants_vcf)
    
  
  output:
     tuple val("$sample_id"),
           val("${reconstructer}"),
           val("${method}"),
           path("haplotypes_aligned.fasta"),
           path("${proportions}"),
           path("${variants_vcf}"),
           val("{$contig_name}")
  script:
    """
    

    cat    ${first_ref} ${haplotypes}  > haplo_plus_ref.fasta
    
    mafft --adjustdirectionaccurately --thread ${task.cpus} --auto --quiet haplo_plus_ref.fasta > haplotypes_aligned.fasta

    """
}
