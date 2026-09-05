#!/usr/bin/env nextflow

process RECONSTRUCT_CLIQUE {
  label 'process_high'

  container "cacciabue/multiquas:clique.v0.0.1"  
  tag "${[sample_id, contig_name].join(':')}"
  containerOptions "--cpus=${task.cpus}"
  
  input:
    
    
  tuple val(sample_id),
        path(first_ref),
        path(general_ref),
        val(reconstructer),
        val(method),
        val(contig_name),
        path(sorted_bam),
        path(sorted_bai),
        path(map_bam),
        path(map_bai),
        path(stats), 
        path(number_reads), 
        path(first_index),
        path(reference)
   
   
        
  
  
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
    def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")
   
    """
    java -Xmx${memory} -jar /programs/clique/clique-snv.jar -threads ${task.cpus} -m snv-illumina -tf 0.1 -t 1000 -in ${map_bam} 
    cp snv_output/*.fasta ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta
    
     cp $stats ${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt

    """
 
}
