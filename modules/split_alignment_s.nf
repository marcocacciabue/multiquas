#!/usr/bin/env nextflow

process SPLIT_ALIGNMENT_S {
  label 'process_low'
  tag "$sample_id"
  container "cacciabue/multiquas:developing"
  publishDir "results/${sample_id}/aligment_s/single", mode: 'copy'
  
  input:
  tuple path(sorted_bam),
  path(sorted_bai),
  path(map_bam),
  path(map_bai),
  path(stats), 
  path(number_reads), 
  path(index),
  path(reference),
  val(sample_id)
  
  
  
  output:
  tuple path("single.s.1.fq"), 
  path("single.s.2.fq"), 
  path("single.o.fq"), 
  val("single"),
  val("$sample_id"), 
  path("${reference}"),
  path("$stats"),emit: splitted_reads
  cpus 4
  
  script:
    """
  

    samtools view $map_bam  -b | samtools sort -n - | samtools view -f2 -b | samtools bam2fq -1 single.s.1.fq -2 single.s.2.fq -0 single.o.fq -s single.s.fq -n -
    
    
    """
}
