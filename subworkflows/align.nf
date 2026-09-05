
include { ALIGNMENT_M                                      } from '../modules/alignment_m.nf'
include { ALIGNMENT_S                                      } from '../modules/alignment_s.nf'
include { VARIANTS                                         } from '../modules/variants.nf'

workflow GET_ALIGNMENT_S {
  
  take:
    reads
  
  main:
    
    ALIGNMENT_S(reads)
    ALIGNMENT_S_SAVE=ALIGNMENT_S.out.single_bam.map { sample_id,
    first_ref,
    general_ref,
    sorted_bam, 
    sorted_bai, 
    map_bam,
    map_bai,
    stats,
    number_reads,
    index  -> [sample_id:sample_id,
                 map_bam:map_bam,
                 map_bai:map_bai,
                 stats:stats]
    }.map { sample -> [sample.sample_id, sample] }
  
  
  emit:
    single_bam =  ALIGNMENT_S.out.single_bam
    reference_list = ALIGNMENT_S.out.reference_list
    save       =  ALIGNMENT_S_SAVE
}


workflow GET_ALIGNMENT_M {
  
  take:
    reads
  
  main:
    
    ALIGNMENT_M(reads)
  
  emit:
    multiple_bam   = ALIGNMENT_M.out.multiple_bam
    reference_list = ALIGNMENT_M.out.reference_list
}

workflow GET_VARIANTS {
  
  take:
    single_bam
  
  main:
    
    VARIANTS(single_bam)
  
  emit:
    variants = VARIANTS.out
}



