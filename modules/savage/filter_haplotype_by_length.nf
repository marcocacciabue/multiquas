#!/usr/bin/env nextflow

process FILTER_HAPLOTYPES_BY_LENGTH {
  label 'process_low'
  container "cacciabue/seqkit"
  tag "${sample_id}"
  containerOptions "--cpus=${task.cpus}"
  input:
    tuple val(sample_id),
  path(first_ref),
  path(general_ref),
  val(reconstructer),
  val(method),
  val(contig_name),
  path(stats),
  path(reconstructedVariants),
  path(s1_fq),
  path(s2_fq)
  
  
  output:
    tuple val("$sample_id"), 
  path("$first_ref"),
  path("$general_ref"),
  val("$reconstructer"),
  val("$method"),
  val("$contig_name"),
  path("$stats"),
  path("${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants_filtered.fasta"),
    path("$s1_fq"),
  path("$s2_fq"),emit: haplotypes
  
  script:
    """
    # TODO add check to see if larger than expected size sequences are found or not.
       seqkit seq  --min-len 100 ${reconstructedVariants} > ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants_filtered.fasta

    
    
    
    

    """
}