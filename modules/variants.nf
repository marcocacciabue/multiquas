#!/usr/bin/env nextflow

process VARIANTS {
  
    label 'process_high'
    errorStrategy  { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries  3

    container "cacciabue/multiquas:tools.v0.0.1"
    
    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(first_ref),
          path(general_ref),
          path(sorted_bam),
          path(sorted_bai),
          path(map_bam),
          path(map_bai),
          path(stats), 
          path(number_reads), 
          path(index)

  
   
    output:
    tuple val("${sample_id}"), 
          path("${sample_id}_variants.vcf"), emit: variants
    script:
    """
     lofreq viterbi -f ${first_ref} -o ${sample_id}_map_viterbi.bam ${map_bam}
     samtools sort  -@ ${task.cpus} ${sample_id}_map_viterbi.bam > ${sample_id}_map_viterbi_sorted.bam
     samtools index -@ ${task.cpus} ${sample_id}_map_viterbi_sorted.bam ${sample_id}_map_viterbi_sorted.bai
     lofreq call-parallel --pp-threads ${task.cpus} --use-orphan  -f ${first_ref} ${sample_id}_map_viterbi_sorted.bam -o ${sample_id}_variants.vcf
    """
}
