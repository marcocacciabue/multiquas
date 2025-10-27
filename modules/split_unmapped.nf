#!/usr/bin/env nextflow

process SPLIT_UNMAPPED {
     tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/single_alignment/unmapped", mode: 'copy'
  

    input:
    tuple path(sorted_bam),
    path(sorted_bai),
    path(map_bam),
    path(map_bai),
    path(stats), 
    path(number_reads), 
    path(index),
    path(reference),
    val(sample_id)
    output:

    
    output:
      tuple path("unmapped.s.1.fq"), path("unmapped.s.2.fq"), path("unmapped.o.fq"), val("unmapped"), val("$sample_id"), path("first.fasta"), emit: merged_reads
    cpus 4

    script:
    def contig = "unmapped"
    """
    samtools faidx ${reference}  
    
    sed -n '1p' *.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq ${reference} seq_names.txt > first.fasta


    samtools view -b -f 12 -F 256 $sorted_bam -b | samtools sort -n - | samtools bam2fq -1 ${contig}.s.1.fq -2 ${contig}.s.2.fq -0 ${contig}.o.fq -s ${contig}.s.fq -n -

    """
}