#!/usr/bin/env nextflow

process SPLIT_ALIGNMENT {
     tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/multiple_alignment/${contig}", mode: 'copy'

    input:
    tuple val(contig),
    path(sorted_bam),
    path(sorted_bai),
    path(map_bam),
    path(map_bai),
    path(stats), 
    path(number_reads), 
    path(index),
    path(reference),
    val(sample_id)
  

    
    output:
      tuple path("${contig}.s.1.fq"), 
      path("${contig}.s.2.fq"), 
      path("${contig}.o.fq"), 
      val("${contig}"),
      val("$sample_id"), 
      path("${contig}.fasta"),
      path("$stats"),emit: splittted_reads
    cpus 4

    script:
    """
    if [[ ${contig} = "unmapped" ]]
    then
    samtools faidx ${reference}  
    
    sed -n '1p' *.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq ${reference} seq_names.txt > ${contig}.fasta


    samtools view -b -f 12 -F 256 $sorted_bam -b | samtools sort -n - | samtools bam2fq -1 ${contig}.s.1.fq -2 ${contig}.s.2.fq -0 ${contig}.o.fq -s ${contig}.s.fq -n -
    else
    echo ${contig}  >  "${contig}.txt"

    seqtk subseq $reference "${contig}.txt" > "${contig}.fasta"

    samtools view $map_bam ${contig}  -b | samtools sort -n - | samtools view -f2 -b | samtools bam2fq -1 ${contig}.s.1.fq -2 ${contig}.s.2.fq -0 ${contig}.o.fq -s ${contig}.s.fq -n -
    
    fi
    """
}