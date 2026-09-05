#!/usr/bin/env nextflow

process GET_READS {
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
          path("${contig_name}.fasta"),
          path("$stats"),
          path("s.1.fq"), 
          path("s.2.fq"), 
          path("o.fq"), emit: splitted_reads

  
  script:
    """
  
     if [[ ${contig_name} = "unmapped" ]]
    then
    samtools faidx ${general_ref}  
    
    sed -n '1p' *.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq ${general_ref} seq_names.txt > ${contig_name}.fasta
    samtools view -b -f 12 -F 256 $sorted_bam -b | samtools sort -n - | samtools bam2fq -1 s.1.fq -2 s.2.fq -0 o.fq -s s.fq -n -
    else
     echo ${contig_name}  >  '${contig_name}.txt'

    seqtk subseq $general_ref '${contig_name}.txt' > "${contig_name}.fasta"
    samtools view $map_bam ${contig_name} -b | samtools sort -n - | samtools view -f2 -b | samtools bam2fq -1 s.1.fq -2 s.2.fq -0 o.fq -s s.fq -n -

    
    fi
    # added this step to avoid name collision down the road
    cp $first_ref  ${sample_id}_${reconstructer}_${contig_name}_first_ref.fasta
    cp $general_ref  ${sample_id}_${reconstructer}_${contig_name}_general_ref.fasta
    
    """
}
