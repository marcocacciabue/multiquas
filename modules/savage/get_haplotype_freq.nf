#!/usr/bin/env nextflow

process GET_HAPLOTYPE_FREQ {
  label 'process_medium'
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
          path("$reconstructedVariants"),
          path("reads_per_haplotype.txt"),emit: haplotypes
  
  script:
    """
    bowtie2-build ${reconstructedVariants} ref 
    bowtie2  -p ${task.cpus} -x ref -1 ${s1_fq}  -2 ${s2_fq} | samtools view -@ ${task.cpus} -bT ${reconstructedVariants} - | samtools sort -@ ${task.cpus} -m 2G - > sorted.bam
    samtools index -@ ${task.cpus} sorted.bam sorted.bai
    samtools idxstats sorted.bam > reads_per_haplotype.txt
    """
}
