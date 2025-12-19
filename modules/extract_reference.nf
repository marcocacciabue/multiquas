#!/usr/bin/env nextflow

process EXTRACT_REFERENCE {
    label 'process_low'
    container "cacciabue/multiquas:developing"

    input:

    path input_ref

    output:
    path ("first.fasta"), emit: first_ref
    path ("${input_ref}"), emit: general_ref

    script:
    """
    samtools faidx ${input_ref}  
    
    sed -n '1p' *.fasta.fai | cut -f 1 > seq_names.txt
    
    seqtk subseq ${input_ref} seq_names.txt > first.fasta
    
 
    """
}
