#!/usr/bin/env nextflow

process RECONSTRUCT_QURE {
    label 'process_high'
    tag "$sample_id"
    container "cacciabue/multiquas:developing"  
    containerOptions "--cpus=${task.cpus}"
    
    input:
    
    
    tuple path(reads), val(sample_id), path(ref), val(contig_name), path(stats)
    
    output:
    
    tuple path("${contig_name}_temporary_reconstructedVariants.txt"),     
       val("$sample_id"),
       val("$contig_name"),
       path("$stats"), emit: haplotypes
    script:
        def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")

    """
   if [[ ${task.attempt} -eq 1 ]]; then
        seqtk sample -s 11 ${reads} 0.9 >  '${contig_name}_temporary.fasta'
   fi
 
   if [[ ${task.attempt} -eq 2 ]]; then
      seqtk sample -s 11 ${reads} 0.4 >  '${contig_name}_temporary.fasta'
   fi
   
   if [[ ${task.attempt} -eq 3 ]]; then
      seqtk sample -s 11 ${reads} 0.08 >  '${contig_name}_temporary.fasta'
   fi
   java -Xmx${memory} -cp /programs/QuRe_v0.99971/ QuRe   '${contig_name}_temporary.fasta' ${ref} 1E-25 0.00035 10
  
 
    """
    stub: 
            def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")

    """
   
      seqtk sample -s 11 ${reads} 0.08 >  '${contig_name}_temporary.fasta'
   
     java -Xmx${memory} -cp /programs/QuRe_v0.99971/ QuRe   '${contig_name}_temporary.fasta' ${ref} 1E-25 0.00035 10
  
 
    """
}