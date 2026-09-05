#!/usr/bin/env nextflow

process EXTRACT_REFERENCE {
    label 'process_low'
    container "cacciabue/multiquas:tools.v0.0.1"


    input:

    path input_ref

    output:
    path("first.fasta"), emit: first_ref
    path("general_ref.fasta"), emit: general_ref

    script:
    """
    cp ${input_ref} general_ref.fasta
    samtools faidx general_ref.fasta  
    
    sed -n '1p' general_ref.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq  general_ref.fasta seq_names.txt > first.fasta
    
 
    """
}
