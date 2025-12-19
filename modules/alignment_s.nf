#!/usr/bin/env nextflow

process ALIGNMENT_S {
    label 'process_medium'
    tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/aligment_s", mode: 'copy'
    cpus 4
    memory "2G"
    input:
    tuple val(sample_id), path(read1), path(read2),  path(first_ref)
   
    output:
    tuple  path("sorted.bam"),
    path("sorted.bai"),
    path("map.bam"),
    path("map.bai"),
    path("stats.txt"), 
    path("number_reads.txt"), 
    path("first_index.tar.gz"),
    path("${first_ref}"),
    val("$sample_id"), emit: single_bam

    script:
    """

    bowtie2-build ${first_ref} ref
  
    bowtie2  -p ${task.cpus} -x ref -1 ${read1}  -2 ${read2} | samtools view -@ ${task.cpus} -bT ${first_ref} - | samtools sort -@ ${task.cpus} -m 2G - > sorted.bam
    samtools index -@ ${task.cpus} sorted.bam sorted.bai
    samtools idxstats sorted.bam > stats.txt
    samtools view -c sorted.bam > number_reads.txt
    samtools view -@ ${task.cpus} -h -F 4 -b sorted.bam > map.bam
    samtools index -@ ${task.cpus} map.bam map.bai
    
    
    tar -czf 'first_index.tar.gz' *.bt2 

   #generate list of references in bam file

    samtools view -H map.bam | grep -P '^@SQ' | cut -f 2 -d ':' | cut -f 1 > references_list.txt
    echo "unmapped" >> references_list.txt
    """
}
