#!/usr/bin/env nextflow

process RECONSTRUCT_QURE {
    label 'process_high'
    container "cacciabue/multiquas:qure.v0.0.1"
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
          path("${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt"), emit: haplotypes
   
    script:
        def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")

    """
    #TODO add step to check if reads are enough.
   if [[ ${task.attempt} -eq 1 ]]; then
        seqtk sample -s 11 ${reads} 0.99 >  temporary.fasta
   fi
 
   if [[ ${task.attempt} -eq 2 ]]; then
      seqtk sample -s 11 ${reads} 0.5 >  temporary.fasta
   fi
   
   if [[ ${task.attempt} -eq 3 ]]; then
      seqtk sample -s 11 ${reads} 0.1 >  temporary.fasta
   fi
   java -Xmx${memory} -cp /programs/QuRe_v0.99971/ QuRe   temporary.fasta ${reference} 1E-25 0.00035 10
   mv ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta
 
    """
    stub: 
            def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")

    """
   
      seqtk sample -s 11 ${reads} 0.1 >  temporary.fasta
   
     java -Xmx${memory} -cp /programs/QuRe_v0.99971/ QuRe   temporary.fasta ${reference} 1E-25 0.00035 10
  mv temporary_reconstructedVariants.txt ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta
 cp $stats ${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt

"""
}
