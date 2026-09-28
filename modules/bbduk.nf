#!/usr/bin/env nextflow

process BBDUK {
    label 'process_low'
    tag "${sample_id}"
    container "cacciabue/multiquas:tools.v0.0.1"
    memory "4G"
    input:
    tuple val(sample_id), 
          path(read1), 
          path(read2),
          path(input_ref)

    output:
  
    tuple val("${sample_id}"),
          path("${sample_id}_R1_trimmed.fq"), 
          path("${sample_id}_R2_trimmed.fq"), 
          path("first.fasta"), 
          path("general_ref.fasta"),emit: trimmed_reads
 

    script:

    """
    cp ${input_ref} general_ref.fasta
    samtools faidx general_ref.fasta  
    
    sed -n '1p' general_ref.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq  general_ref.fasta seq_names.txt > first.fasta
    bbduk.sh in1=${read1} out1='${sample_id}_R1_trimmed.fq' in2=${read2} out2='${sample_id}_R2_trimmed.fq' ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl ignorebadquality

    """
}
