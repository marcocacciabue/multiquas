include { PREPARE_ALIGNMENT    } from '../modules/prepare_alignment.nf'

include { RECONSTRUCT_VIQUAS   } from '../modules/viquas/reconstruct_viquas.nf'
include { FILTER_HAPLOTYPES_BY_FREQ   } from '../modules/viquas/filter_haplotypes_by_freq.nf'

include { ADJUST_HAPLOTYPE_FREQ } from '../modules/adjust_haplotype_freq.nf'
include { MERGE_RECONSTRUCTED                      } from '../modules/merge_reconstructed.nf'

include { ALIGN_HAPLOTYPES      } from '../modules/align_haplotypes.nf'
include { HAPLOTYPE_EVALUATOR   } from '../modules/haplotype_evaluator.nf'

workflow GET_VIQUAS_M {
  
  take:
    single_bam
  variants
  
  
  main:
    
    
    PREPARE_ALIGNMENT(single_bam,
                      "viquas_m",
                      "multiple")
  
  
  RECONSTRUCT_VIQUAS(PREPARE_ALIGNMENT.out,
                     5000)
  FILTER_HAPLOTYPES_BY_FREQ(RECONSTRUCT_VIQUAS.out)
  ADJUST_HAPLOTYPE_FREQ(FILTER_HAPLOTYPES_BY_FREQ.out)

  
  MERGE_RECONSTRUCTED(ADJUST_HAPLOTYPE_FREQ.out.haplotypes
                      .groupTuple(by:0))
  haplo_ch_clique_m = MERGE_RECONSTRUCTED.out
  .combine(variants,by:0)
  
  ALIGN_HAPLOTYPES(haplo_ch_clique_m)
  HAPLOTYPE_EVALUATOR(ALIGN_HAPLOTYPES.out)
  
  emit:
    reconstructed_data = HAPLOTYPE_EVALUATOR.out.reconstructed_data
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
  results            = HAPLOTYPE_EVALUATOR.out.results
}