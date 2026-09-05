include { GET_READS             } from '../modules/get_reads.nf'
include { MERGE_READS           } from '../modules/merge_reads.nf'
include { RECONSTRUCT_QURE      } from '../modules/qure/reconstruct_qure.nf'
include { ADJUST_HAPLOTYPE_FREQ } from '../modules/adjust_haplotype_freq.nf'
include { ALIGN_HAPLOTYPES      } from '../modules/align_haplotypes.nf'
include { HAPLOTYPE_EVALUATOR   } from '../modules/haplotype_evaluator.nf'


workflow GET_QURE_S {
  
  take:
    single_bam
    variants
  
  main:
    

    GET_READS(single_bam,
                     "qure_s",
                     "single")
    MERGE_READS(GET_READS.out.splitted_reads)
    RECONSTRUCT_QURE(MERGE_READS.out.reads)
  
    ADJUST_HAPLOTYPE_FREQ(RECONSTRUCT_QURE.out)
  
  
    haplo_ch_single = ADJUST_HAPLOTYPE_FREQ.out
    .combine(variants,by:0)

  
  
  ALIGN_HAPLOTYPES(haplo_ch_single)
  HAPLOTYPE_EVALUATOR(ALIGN_HAPLOTYPES.out)
  
  
  emit:
    reconstructed_data =  HAPLOTYPE_EVALUATOR.out.reconstructed_data
   .map { sample_id, reconstructer, method,haplo, proportions,
    graphs,r_sq -> [sample_id:sample_id,
                    reconstructer:reconstructer,
                    method:method,
                    haplo:haplo,
                    proportions:proportions,
                    graphs:graphs,
                    r_sq:r_sq.text.toFloat()]
  }
  .map { sample -> [sample.sample_id, sample] }
  results            =  HAPLOTYPE_EVALUATOR.out.results
}
