include { GET_READS          } from '../modules/get_reads.nf'
include { MERGE_READS              } from '../modules/merge_reads.nf'

include { RECONSTRUCT_SAVAGE    } from '../modules/savage/reconstruct_savage.nf'
include { FILTER_HAPLOTYPES_BY_LENGTH    } from '../modules/savage/filter_haplotype_by_length.nf'

include { GET_HAPLOTYPE_FREQ    } from '../modules/savage/get_haplotype_freq.nf'

include { ADJUST_HAPLOTYPE_FREQ } from '../modules/adjust_haplotype_freq.nf'
include { MERGE_RECONSTRUCTED                         } from '../modules/merge_reconstructed.nf'

include { ALIGN_HAPLOTYPES      } from '../modules/align_haplotypes.nf'
include { HAPLOTYPE_EVALUATOR   } from '../modules/haplotype_evaluator.nf'

workflow GET_SAVAGE_M {
  
  take:
    multiple_bam
  variants
  
  
  main:
    
    
    GET_READS(multiple_bam,
              "savage_m",
              "multiple")
  //skip merge reads because savage needs fq files.
  // MERGE_READS(GET_READS.out.splitted_reads)
  RECONSTRUCT_SAVAGE(GET_READS.out.splitted_reads)
  FILTER_HAPLOTYPES_BY_LENGTH(RECONSTRUCT_SAVAGE.out)
  GET_HAPLOTYPE_FREQ(FILTER_HAPLOTYPES_BY_LENGTH.out)
  ADJUST_HAPLOTYPE_FREQ(GET_HAPLOTYPE_FREQ.out)
  MERGE_RECONSTRUCTED(ADJUST_HAPLOTYPE_FREQ.out.haplotypes.groupTuple(by:0))

  
  haplo_ch_multiple = MERGE_RECONSTRUCTED.out
                        .combine(variants,by:0)
  
    ALIGN_HAPLOTYPES(haplo_ch_multiple)
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