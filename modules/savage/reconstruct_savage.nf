#!/usr/bin/env nextflow

process RECONSTRUCT_SAVAGE {
  label 'process_high'
  container "cacciabue/multiquas:haploconduct.v0.0.1"
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
          path(stats),
          path(s1_fq),
          path(s2_fq),
          path(o_fq)
  
  
  
  
  output:
    tuple val("$sample_id"), 
          path("$first_ref"),
          path("$general_ref"),
          val("$reconstructer"),
          val("$method"),
          val("$contig_name"),
          path("$stats"),
          path("${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta"),
          path("$s1_fq"),
          path("$s2_fq"),emit: haplotypes
  
  script:
    
    
    """
    #only considers paired reads
    reads=\$(seqtk size $s1_fq | cut -f1)
    frac=\$(echo "scale=0; (2*\$reads/750)" | bc)
    
    conda run -n env_haploconduct haploconduct savage -s $o_fq -p1 $s1_fq  -p2 $s2_fq  --split \$frac

    cp contigs_stage_c.fasta ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta
   

 
    """
  
}
