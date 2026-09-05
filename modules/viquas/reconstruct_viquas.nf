#!/usr/bin/env nextflow

process RECONSTRUCT_VIQUAS {
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
        path(sorted_bam),
        path(sorted_bai),
        path(map_bam),
        path(map_bai),
        path(stats), 
        path(number_reads), 
        path(first_index),
        path(reference)
  val(max_reads)
   
        
  
  
  output:
    tuple val("$sample_id"),
          path("$first_ref"),
          path("$general_ref"),
          val("$reconstructer"),
          val("$method"),
          val("$contig_name"),
          path("$stats"),
          path("ViQuaS-Spectrum.fa"),
          path("${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt"),emit: haplotypes
  
  script:
    def memory = "${task.memory}".replaceAll("\\s","").replaceAll("B","")
   
    """
    #get number of mapped reads to all references (it assumes that only relevant reads are present)
    reads=\$(samtools idxstats ${map_bam} | cut -f3  | awk '{sum += \$1} END {print sum}')

    if [[ \$reads -gt $max_reads ]]
    then

    frac=\$(echo "scale=2; ($max_reads/\$reads)+12" | bc)

    #downsample bam file to avoid failure in process. 1 is seed. 
    samtools view -b -s \$frac ${map_bam} > map_reduced.bam
    samtools index map_reduced.bam
    Rscript /programs/viquas/ViQuaS.R ${reference} map_reduced.bam  
    else
    #use all reads
    Rscript /programs/viquas/ViQuaS.R ${reference} ${map_bam} 

    fi
   
    
    ## Viquas produces too many haplotypes
    ## keep only those of freq >= 0.02
    ## this step should be a different process in case user wants to change threshold and not repeat the reconstruction step
  ##  threshold=0.02
  ##  samtools faidx ViQuaS-Spectrum.fa
  ##  cat ViQuaS-Spectrum.fa.fai | cut -f 2 -d ':' | cut -f 1 | awk -F'_' -v t="\$threshold" '\$2+0 >= t' > haplo_list.txt 
    ##seqtk subseq ViQuaS-Spectrum.fa haplo_list.txt  > ${sample_id}_${reconstructer}_${method}_${contig_name}_reconstructed_variants.fasta

    
    
    cp $stats ${sample_id}_${reconstructer}_${method}_${contig_name}_stats_dummy.txt
    
    

    """
 
}
