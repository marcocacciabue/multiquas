#!/usr/bin/env nextflow

process VARIANTS {
    tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/single_alignment/variants", mode: 'copy'
    cpus 4
    memory "2G"
    input:
    tuple     path(sorted_bam),
    path(sorted_bai),
    path(map_bam),
    path(map_bai),
    path(stats), 
    path(number_reads), 
    path(index),
    path(reference),
    val(sample_id)
  
   
    output:
    tuple path("${sample_id}_variants.vcf"), val("${sample_id}"), emit: variants
    script:
    """
     lofreq viterbi -f ${reference} -o ${sample_id}_map_viterbi.bam ${map_bam}
     samtools sort  -@ ${task.cpus} ${sample_id}_map_viterbi.bam > ${sample_id}_map_viterbi_sorted.bam
  samtools index -@ ${task.cpus} ${sample_id}_map_viterbi_sorted.bam ${sample_id}_map_viterbi_sorted.bai
  lofreq call-parallel --pp-threads ${task.cpus} --use-orphan  -f ${reference} ${sample_id}_map_viterbi_sorted.bam -o ${sample_id}_variants.vcf
    """
}