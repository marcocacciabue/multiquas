#!/usr/bin/env nextflow

process EXTRACT_REFERENCE {
    tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/", mode: 'copy'

    cpus 4
    
    input:
    path(input_ref)
    
    output:
    path("first.fasta"), emit:first_ref
    script:
    """
    samtools faidx ${input_ref}  
    
    sed -n '1p' *.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq ${input_ref} seq_names.txt > first.fasta
    
 
    """
}