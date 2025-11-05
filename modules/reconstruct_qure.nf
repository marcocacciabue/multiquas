#!/usr/bin/env nextflow

process RECONSTRUCT_QURE {
    tag "$sample_id"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/haplotypes_qure/${contig_name}", mode: 'copy'    
    memory '6G'
    errorStrategy { task.attempt <3 ? 'retry' : 'ignore' }
    maxRetries 3
    input:
    
    tuple path(reads), val(sample_id), path(ref), val(contig_name), path(stats)
    
    output:
    
    tuple path("${sample_id}_${contig_name}_temporary_reconstructedVariants.txt"), 
      path("${contig_name}_adjusted_proportions.txt"),
       val("$sample_id"), emit: haplotypes
    script:
   
    """
   if [[ ${task.attempt} -eq 1 ]]; then
        seqtk sample -s 11 ${reads} 0.1 > '${sample_id}_${contig_name}_temporary.fasta'
     java -cp /programs/QuRe_v0.99971/  QuRe  '${sample_id}_${contig_name}_temporary.fasta' ${ref} 1E-25 0.00035 10
   fi
 
   if [[ ${task.attempt} -eq 2 ]]; then
      seqtk sample -s 11 ${reads} 0.05 > '${sample_id}_${contig_name}_temporary.fasta'
     java -cp /programs/QuRe_v0.99971/  QuRe '${sample_id}_${contig_name}_temporary.fasta' ${ref} 1E-25 0.00035 10
   fi
   
   if [[ ${task.attempt} -eq 3 ]]; then
      seqtk sample -s 11 ${reads} 0.01 > '${sample_id}_${contig_name}_temporary.fasta'
     java -cp /programs/QuRe_v0.99971/  QuRe '${sample_id}_${contig_name}_temporary.fasta' ${ref} 1E-25 0.00035 10
   fi
   Rscript /root/data/adjust_prop.R  ${stats} "${sample_id}_${contig_name}_temporary_reconstructedVariants.txt" ${contig_name} ${contig_name}_adjusted_proportions.txt
    """
}