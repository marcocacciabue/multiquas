#!/usr/bin/env nextflow

process PREPARE_ALIGNMENT {
  label 'process_low'
  
  container "cacciabue/multiquas:tools.v0.0.1"
  tag "${[sample_id, contig_name].join(':')}"
  
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
          path(index),
          val(contig_name)
    val(reconstructer)
    val(method)
  
  output:
    tuple val("$sample_id"),
          path("${sample_id}_${reconstructer}_${contig_name}_first_ref.fasta"),
          path("${sample_id}_${reconstructer}_${contig_name}_general_ref.fasta"),
          val("$reconstructer"),
          val("$method"),
          val("$contig_name"),
          path("sorted_new.bam"), 
          path("sorted_new.bai"), 
          path("map_new.bam"), 
          path("map_new.bai"),
          path("${stats}"), 
          path("${number_reads}"),
          path("$index"),
          path("contig.fasta")
  
  
  script:
    """
      if [[ ${contig_name} = "unmapped" && ${method} = "multiple" ]]
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
    
    fi
    if [[ ${contig_name} != "unmapped" && ${method} = "multiple" ]]
    then
    echo ${contig_name}  >  contig.txt

    seqtk subseq ${general_ref} contig.txt > contig.fasta

    samtools view $map_bam ${contig_name}  -b > sorted_new.bam
    samtools index -@ ${task.cpus} sorted_new.bam sorted_new.bai

    samtools view -h -F 4 -b sorted_new.bam > map_new.bam
    samtools index -@ ${task.cpus} map_new.bam map_new.bai

    fi
    if [[  ${method} = "single" ]]
    then
    
    echo ${contig_name}  >  contig.txt

    seqtk subseq ${general_ref} contig.txt > contig.fasta
    
  
    cp $map_bam  map_new.bam
    cp $sorted_bam  sorted_new.bam
    cp $map_bai  map_new.bai
    cp $sorted_bam  sorted_new.bai
    fi
    
    # added this step to avoid name collision down the road
    cp $first_ref  ${sample_id}_${reconstructer}_${contig_name}_first_ref.fasta
    cp $general_ref  ${sample_id}_${reconstructer}_${contig_name}_general_ref.fasta

    """
}
