#!/usr/bin/env nextflow

process RECONSTRUCT_CLIQUE_S {
  label 'process_medium'
  tag "$sample_id"
  def reconstructer = "clique"
  container "cacciabue/multiquas:developing"  

  containerOptions "--cpus=${task.cpus}"
  
  input:
    
    
  tuple path(sorted_bam),
        path(sorted_bai),
        path(map_bam),
        path(map_bai),
        path(stats), 
        path(number_reads), 
        path(first_index),
        path(first_ref),
        val(sample_id)
  
  
  output:
    tuple val("$sample_id"),
          val("$reconstructer"),
          path("${sample_id}_reconstructed_variants_clique.fasta"),
          val("single"),
          path("$stats"), emit: haplotypes
  
  script:
    def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")
   
    """
    java -Xmx${memory} -jar /programs/clique-snv.jar -threads ${task.cpus} -m snv-illumina -tf 0.1 -t 1000 -in ${map_bam} 
    cp snv_output/map.fasta ${sample_id}_reconstructed_variants_clique.fasta
    """
 
}