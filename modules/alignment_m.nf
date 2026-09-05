#!/usr/bin/env nextflow

process ALIGNMENT_M {
    label 'process_medium'

    container "cacciabue/multiquas:tools.v0.0.1"
        tag "${sample_id}"
    containerOptions "--cpus=${task.cpus}"
    input:
    tuple val(sample_id), 
          path(read1), 
          path(read2),  
          path(first_ref),
          path(general_ref)
   
    output:
    tuple val("$sample_id"),
          path("${first_ref}"),
          path("$general_ref"),
          path("sorted.bam"),
          path("sorted.bai"),
          path("map.bam"),
          path("map.bai"),
          path("stats.txt"), 
          path("number_reads.txt"), 
          path("${general_ref.simpleName}_index.tar.gz"), emit: multiple_bam
    tuple val("$sample_id"),     
          path("references_list.txt"),emit: reference_list
    script:
    """
    bowtie2-build $general_ref ref
  
    bowtie2 --no-discordant --no-mixed -p ${task.cpus} -x ref -1 ${read1}  -2 ${read2} | samtools view -@ ${task.cpus} -bT $general_ref - | samtools sort -@ ${task.cpus} -m 2G - > sorted.bam
    #bowtie2 -p ${task.cpus} -x ref -1 ${read1}  -2 ${read2} | samtools view -@ ${task.cpus} -bT $general_ref - | samtools sort -@ ${task.cpus} -m 2G - > sorted.bam

    samtools index -@ ${task.cpus} sorted.bam sorted.bai
    samtools idxstats sorted.bam > stats.txt
    samtools view -c sorted.bam > number_reads.txt
    samtools view -@ ${task.cpus} -h -F 4 -b sorted.bam > map.bam
    samtools index -@ ${task.cpus} map.bam map.bai
    
    
    tar -czf '${general_ref.simpleName}_index.tar.gz' *.bt2 

   #generate list of references in bam file

    samtools view -H map.bam | grep -P '^@SQ' | cut -f 2 -d ':' | cut -f 1 > references_list.txt
    echo "unmapped" >> references_list.txt
    """
}
