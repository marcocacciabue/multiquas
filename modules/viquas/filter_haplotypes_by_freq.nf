#!/usr/bin/env nextflow

process FILTER_HAPLOTYPES_BY_FREQ {
  label 'process_high'
  
  container "cacciabue/multiquas:viquas.v0.0.1"  
  tag "${[sample_id, contig_name].join(':')}"
  containerOptions "--cpus=${task.cpus}"
  
  input:
    
    
    tuple val(sample_id),
          path(first_ref),
          path(general_ref),
          val(reconstructer),
          val(method),
          val(contig_name),
          path(stats), 
          path(reconstructed_variants),
          path(dummy_stats)
  
  
  
  
  output:
    tuple val("$sample_id"),
          path("$first_ref"),
          path("$general_ref"),
          val("$reconstructer"),
          val("$method"),
          val("$contig_name"),
          path("$stats"),
          path("${sample_id}_${reconstructer}_${method}_${contig_name}_filtered_reconstructed_variants.fasta"),
          path("$dummy_stats"),emit: haplotypes
  
  script:

  """
    ## Viquas produces too many haplotypes
    ## keep only those of freq >= 0.02
    ## this step should be a different process in case user wants to change threshold and not repeat the reconstruction step
    threshold=0.02
    samtools faidx $reconstructed_variants
    cat ${reconstructed_variants}.fai | cut -f 2 -d ':' | cut -f 1 | awk -F'_' -v t="\$threshold" '\$2+0 >= t' > haplo_list.txt 
    seqtk subseq ViQuaS-Spectrum.fa haplo_list.txt  > ${sample_id}_${reconstructer}_${method}_${contig_name}_filtered_reconstructed_variants.fasta
    

    """
  
}