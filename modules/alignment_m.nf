#!/usr/bin/env nextflow

process ALIGNMENT_M {
    label 'process_medium'

    tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/aligment_m", mode: 'copy'
    cpus 4
    memory "2G"
    input:
    tuple val(sample_id), path(read1), path(read2),  path(input_ref)
   
    output:
    tuple  path("sorted.bam"),
    path("sorted.bai"),
    path("map.bam"),
    path("map.bai"),
    path("stats.txt"), 
    path("number_reads.txt"), 
    path("${input_ref.simpleName}_index.tar.gz"),
    path("$input_ref"),
    val("$sample_id"), emit: multiple_bam
    path("references_list.txt"),emit: reference_list
    script:
    """
    bowtie2-build $input_ref ref
  
    bowtie2 --no-discordant --no-mixed -p ${task.cpus} -x ref -1 ${read1}  -2 ${read2} | samtools view -@ ${task.cpus} -bT $input_ref - | samtools sort -@ ${task.cpus} -m 2G - > sorted.bam
    samtools index -@ ${task.cpus} sorted.bam sorted.bai
    samtools idxstats sorted.bam > stats.txt
    samtools view -c sorted.bam > number_reads.txt
    samtools view -@ ${task.cpus} -h -F 4 -b sorted.bam > map.bam
    samtools index -@ ${task.cpus} map.bam map.bai
    
    
    tar -czf '${input_ref.simpleName}_index.tar.gz' *.bt2 

   #generate list of references in bam file

    samtools view -H map.bam | grep -P '^@SQ' | cut -f 2 -d ':' | cut -f 1 > references_list.txt
    echo "unmapped" >> references_list.txt
    """
}
