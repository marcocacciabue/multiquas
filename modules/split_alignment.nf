#!/usr/bin/env nextflow

process SPLIT_ALIGNMENT {
  label 'process_medium'

  container "cacciabue/multiquas:tools.v0.0.1"
    tag "${sample_id}"
  containerOptions "--cpus=${task.cpus}"
  
  input:
    tuple val(contig_name),
          val(sample_id),
          path(sorted_bam),
          path(sorted_bai),
          path(map_bam),
          path(map_bai),
          path(stats), 
          path(number_reads), 
          path(index),
          path(general_ref)
          
    
  output:
    tuple val("${contig_name}"),
          val("$sample_id"),
          path("sorted_new.bam"), 
          path("sorted_new.bai"), 
          path("map_new.bam"), 
          path("map_new.bai"),
          path("${stats}"), 
          path("${number_reads}"),
          path("$index"),
          path("contig.fasta"),emit: splitted_bams

  script:
    """
    if [[ ${contig_name} = "unmapped" ]]
    then
    samtools faidx ${general_ref}  
    
    sed -n '1p' ${general_ref.simpleName}.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq ${general_ref} seq_names.txt > contig.fasta
    
    samtools view -b -f 12 -F 256 $sorted_bam -b | samtools sort -n - | samtools bam2fq -1 s.1.fq -2 s.2.fq -0 o.fq -s s.fq -n -
    
    bowtie2-build contig.fasta temp_ref
    
    bowtie2  -p ${task.cpus} -x temp_ref -1 s.1.fq  -2 s.2.fq | samtools view -@ ${task.cpus} -bT contig.fasta - | samtools sort -@ ${task.cpus} -m 2G - > sorted_new.bam
    
    samtools index -@ ${task.cpus} sorted_new.bam sorted_new.bai
    samtools view -@ ${task.cpus} -h -F 4 -b sorted_new.bam > map_new.bam
    samtools index -@ ${task.cpus} map_new.bam map_new.bai
    
    else
    echo ${contig_name}  >  contig.txt

    seqtk subseq ${general_ref} contig.txt > contig.fasta

    samtools view $map_bam ${contig_name}  -b > sorted_new.bam
    samtools index -@ ${task.cpus} sorted_new.bam sorted_new.bai

    samtools view -h -F 4 -b sorted_new.bam > map_new.bam
    samtools index -@ ${task.cpus} map_new.bam map_new.bai

    fi
    """
}
