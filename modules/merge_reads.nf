#!/usr/bin/env nextflow

process MERGE_READS {
      label 'process_medium'

  container "cacciabue/multiquas:tools.v0.0.1"
  tag "${[sample_id, contig_name].join(':')}"

  input:
    tuple val(sample_id),
          path(first_ref),
          path(general_ref),
          val(reconstructer),
          val(method),
          val(contig_name), 
          path(reference), 
          path(stats),
          path(s_1_fq), 
          path(s_2_fq), 
          path(o_fq) 
          
          
  
  output:
    
    tuple val("$sample_id"), 
          path("$first_ref"),
          path("$general_ref"),
          val("$reconstructer"), 
          val("$method"), 
          val("$contig_name"),
          path("$reference"), 
          path("${contig_name}_all.fasta"),
          path("$stats"), emit: reads
  
  script:
    """
 
    pear -f $s_1_fq -r $s_2_fq -j ${task.cpus} -o out
    cat out.assembled.fastq out.unassembled.reverse.fastq out.unassembled.forward.fastq | seqtk seq -a - >  '${contig_name}_all.fasta'



    """
}
