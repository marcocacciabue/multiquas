include { PREPARE_ALIGNMENT    } from '../modules/prepare_alignment.nf'

include { RECONSTRUCT_CLIQUE    } from '../modules/clique/reconstruct_clique.nf'
include { ADJUST_HAPLOTYPE_FREQ } from '../modules/adjust_haplotype_freq.nf'
include { ALIGN_HAPLOTYPES      } from '../modules/align_haplotypes.nf'
include { HAPLOTYPE_EVALUATOR   } from '../modules/haplotype_evaluator.nf'

workflow GET_CLIQUE_S {
  
  take:
    single_bam
    variants

  
  main:
  

    PREPARE_ALIGNMENT(single_bam,
                     "clique_s",
                     "single")
    
    
    RECONSTRUCT_CLIQUE(PREPARE_ALIGNMENT.out)
  
    ADJUST_HAPLOTYPE_FREQ(RECONSTRUCT_CLIQUE.out)
  
  
  haplo_ch_clique = ADJUST_HAPLOTYPE_FREQ.out
  .combine(variants,by:0)

  
  
  ALIGN_HAPLOTYPES(haplo_ch_clique)
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