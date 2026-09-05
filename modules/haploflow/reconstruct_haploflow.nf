#!/usr/bin/env nextflow

process RECONSTRUCT_HAPLOFLOW {
    label 'process_high'
    container "cacciabue/multiquas:haploflow.v0.0.1"
    tag "${[sample_id, contig_name].join(':')}"

    containerOptions "--cpus=${task.cpus}"
    
    input:
    tuple val(sample_id),
          path(first_ref),
          path(general_ref),
          val(reconstructer),
          val(method),
          val(contig_name), 
          path(reference),
          path(reads),
          path(stats)
 
        


    output:
    tuple val("$sample_id"), 
          path("$first_ref"),
          path("$general_ref"),
          val("$reconstructer"),
          val("$method"),
          val("$contig_name"),
          path("$stats"),
          path("${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta"),
          path("${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt"),emit: haplotypes
   
    script:
        

    """


    conda run -n env_haploflow haploflow --read-file ${reads} --out out --filter 2000 --long 1
    
    mv out/contigs.fa ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta
 cp $stats ${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt


 
    """
    
}
