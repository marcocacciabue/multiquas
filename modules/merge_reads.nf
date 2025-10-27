#!/usr/bin/env nextflow

process MERGE_READS {
     tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/single_alignment/${contig_name}", mode: 'copy'

    input:
    tuple path(s_1_fq), path(s_2_fq), path(o_fq), val(contig_name), val(sample_id), path(contig), path(stats)
 
    output:
    
    tuple path("${contig_name}_all.fasta"), val("$sample_id"), path("$contig"), val("$contig_name"), path("$stats"), emit: reads
    cpus 4

    script:
    """
 
    pear -f $s_1_fq -r $s_2_fq -j ${task.cpus} -o out
    cat out.assembled.fastq out.unassembled.reverse.fastq out.unassembled.forward.fastq  | seqtk seq -a - >  '${contig_name}_all.fasta'



    """
}