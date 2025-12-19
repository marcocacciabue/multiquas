#!/usr/bin/env nextflow
params.version = 0.1

def Logo() {
    log.info"""

########################################################################################
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MULTIQUAS     -      Multiple reference quasispecies reconstruction pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################

                
             
             A Comprehensive Workflow for the Reconstruction of Viral Haplotypes 
                                from Short Read Data
   

                     #############################################
                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              MULTIQUAS  ~  version ${params.version}
                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     #############################################
""".stripIndent()
}

Logo()
// Primary input
params.input_csv = "samples.csv"
params.input_ref = "data/referencia/referencia.fasta"
// params.qure_time = '1m'
//params.input_ref2 = "data/referencia/referencia.fasta"

// Module INCLUDE statements
include { FASTQC                                           } from './modules/fastqc.nf'
include { BBDUK                                            } from './modules/bbduk.nf'
include { EXTRACT_REFERENCE                                } from './modules/extract_reference.nf'
include { ALIGNMENT_M                                      } from './modules/alignment_m.nf'
include { ALIGNMENT_S                                      } from './modules/alignment_s.nf'
include { SPLIT_ALIGNMENT_M                                } from './modules/split_alignment_m.nf'
include { SPLIT_ALIGNMENT_S                                } from './modules/split_alignment_s.nf'
include { MERGE_READS_M                                    } from './modules/merge_reads_m.nf'
include { MERGE_READS_S                                    } from './modules/merge_reads_s.nf'
include { RECONSTRUCT_QURE as RECONSTRUCT_QURE_M           } from './modules/reconstruct_qure.nf'
include { RECONSTRUCT_QURE as RECONSTRUCT_QURE_S           } from './modules/reconstruct_qure.nf'
include { MERGE_RECONSTRUCTED_QURE                         } from './modules/merge_reconstructed_qure.nf'
include { VARIANTS                                         } from './modules/variants.nf'
include { ADJUST_HAPLOTYPE_FREQ as ADJUST_HAPLOTYPE_FREQ_M } from './modules/adjust_haplotype_freq.nf'
include { ADJUST_HAPLOTYPE_FREQ as ADJUST_HAPLOTYPE_FREQ_S } from './modules/adjust_haplotype_freq.nf'
include { ALIGN_HAPLOTYPES_M                               } from './modules/align_haplotypes_m.nf'
include { ALIGN_HAPLOTYPES_S                               } from './modules/align_haplotypes_s.nf'
include { HAPLOTYPE_EVALUATOR_S                            } from './modules/haplotype_evaluator_s.nf'
include { HAPLOTYPE_EVALUATOR_M                            } from './modules/haplotype_evaluator_m.nf'
include { SELECT_BEST_RECONSTRUCTION                       } from './modules/select_best_reconstruction.nf'
include { RECONSTRUCT_CLIQUE_S                             } from './modules/clique/reconstruct_clique_s.nf'
include { ADJUST_HAPLOTYPE_FREQ_CLIQUE_S                   } from './modules/clique/adjust_haplotype_freq_clique_s.nf'
include { ALIGN_HAPLOTYPES_CLIQUE_S                        } from './modules/clique/align_haplotypes_clique_s.nf'
include { HAPLOTYPE_EVALUATOR_CLIQUE_S                     } from './modules/clique/haplotype_evaluator_clique_s.nf'



workflow GET_QURE_S {
  
  take:
  alignment_s
  reference
  variants
  
  main:
   SPLIT_ALIGNMENT_S(alignment_s)
   MERGE_READS_S(SPLIT_ALIGNMENT_S.out.splitted_reads)
   RECONSTRUCT_QURE_S(MERGE_READS_S.out.reads)

   ADJUST_HAPLOTYPE_FREQ_S(RECONSTRUCT_QURE_S.out)
   
       
   haplo_ch_single = ADJUST_HAPLOTYPE_FREQ_S.out
   .combine(variants,by:0)
   .combine(reference)
   

   ALIGN_HAPLOTYPES_S(haplo_ch_single)
   HAPLOTYPE_EVALUATOR_S(ALIGN_HAPLOTYPES_S.out)
   
   
  emit:
        reconstructed_data=  HAPLOTYPE_EVALUATOR_S.out.reconstructed_data

   results = HAPLOTYPE_EVALUATOR_S.out.results
}




workflow GET_QURE_M {
  
  take:
  alignment_m
  reference_list
  reference
  variants
  
  main:
  split_reads_ch = reference_list.first().splitCsv(strip:true)
	     .combine(alignment_m)

  
  SPLIT_ALIGNMENT_M(split_reads_ch)

  MERGE_READS_M(SPLIT_ALIGNMENT_M.out.splitted_reads)
  RECONSTRUCT_QURE_M(MERGE_READS_M.out.reads)
  ADJUST_HAPLOTYPE_FREQ_M(RECONSTRUCT_QURE_M.out)
  MERGE_RECONSTRUCTED_QURE(ADJUST_HAPLOTYPE_FREQ_M.out.haplotypes.groupTuple(by:0))

 
  haplo_ch_multiple = MERGE_RECONSTRUCTED_QURE.out
  .combine(variants,by:0)
  .combine(reference)


  ALIGN_HAPLOTYPES_M(haplo_ch_multiple)
  HAPLOTYPE_EVALUATOR_M(ALIGN_HAPLOTYPES_M.out)
   
   
  emit:
     reconstructed_data=  HAPLOTYPE_EVALUATOR_M.out.reconstructed_data
     results = HAPLOTYPE_EVALUATOR_M.out.results
}



workflow GET_CLIQUE_S {
  
  take:
  alignment_s
  reference
  variants
  
  main:

   RECONSTRUCT_CLIQUE_S(alignment_s)

   ADJUST_HAPLOTYPE_FREQ_CLIQUE_S(RECONSTRUCT_CLIQUE_S.out)
   
       
   haplo_ch_clique_s = ADJUST_HAPLOTYPE_FREQ_CLIQUE_S.out
   .combine(variants,by:0)
   .combine(reference)
   

   ALIGN_HAPLOTYPES_CLIQUE_S(haplo_ch_clique_s)
   HAPLOTYPE_EVALUATOR_CLIQUE_S(ALIGN_HAPLOTYPES_CLIQUE_S.out)
   
   
  emit:
      reconstructed_data =  HAPLOTYPE_EVALUATOR_CLIQUE_S.out.reconstructed_data

   results = HAPLOTYPE_EVALUATOR_CLIQUE_S.out.results
}
workflow {

    // Create input channel from the contents of a CSV file
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, file(row.fastq_1), file(row.fastq_2)]}
        
     ref_ch =Channel.fromPath(params.input_ref)
     EXTRACT_REFERENCE(ref_ch)
   

  // Adapter trimming and post-trimming QC
   BBDUK(read_ch)


   reads_trimmed_ch_single = BBDUK.out.combine( EXTRACT_REFERENCE.out.first_ref)
   reads_trimmed_ch_multiple = BBDUK.out.combine(EXTRACT_REFERENCE.out.general_ref)

   ALIGNMENT_S(reads_trimmed_ch_single)
   VARIANTS(ALIGNMENT_S.out.single_bam)

  GET_QURE_S(ALIGNMENT_S.out.single_bam,
                       EXTRACT_REFERENCE.out.first_ref,
                       VARIANTS.out)
   
   // GET_CLIQUE_S(ALIGNMENT_S.out.single_bam,
     //                    EXTRACT_REFERENCE.out.first_ref,
       //                  VARIANTS.out)
   

   
    // Run aligment of the reads over all references available
  ALIGNMENT_M(reads_trimmed_ch_multiple)
 
  GET_QURE_M(ALIGNMENT_M.out.multiple_bam,
 ALIGNMENT_M.out.reference_list,
                       EXTRACT_REFERENCE.out.first_ref,
                     VARIANTS.out)
  
  
  
  GET_QURE_M.out.results
  .mix(  GET_QURE_S.out.results)
    .collectFile(name: 'run_results.txt', 
   newLine: false,
   keepHeader:true,
   storeDir:'results')
  .subscribe { file ->
    println "Summary is saved to file: $file"
     println "${file.text}"
  }
  
  result_qure_ch_m=GET_QURE_M.out.reconstructed_data
  .map { sample_id, haplo, proportions,graphs,r_sq -> [sample_id:sample_id,
  haplo:haplo,
  proportions:proportions,
  graphs:graphs,
 r_sq:r_sq.text.toFloat()]
 }
  .map { sample -> [sample.sample_id, sample] }
 
  
  result_qure_ch_s=GET_QURE_S.out.reconstructed_data
  .map { sample_id, haplo, proportions, graphs, r_sq -> [sample_id:sample_id,
  haplo:haplo,
  proportions:proportions,
  graphs:graphs,
  r_sq:r_sq.text.toFloat()]
  }
  .map { sample -> [sample.sample_id, sample] }
  
 //  result_clique_ch_s=GET_CLIQUE_S.out.reconstructed_data
//  .map { sample_id, haplo, proportions, graphs, r_sq -> [sample_id:sample_id,
//  haplo:haplo,
//  proportions:proportions,
//  graphs:graphs,
//  r_sq:r_sq.text.toFloat()]
//  }
  //.map { sample -> [sample.sample_id, sample] }
   
 data_ch=result_qure_ch_s
 .join(result_qure_ch_m)
   .map{ 
     def getMAX = [b.r_sq,c.r_sq].max(),
     def elementMAX = (getMAX == b.r_sq) ? b : c,
     a,b,c -> 
     [elementMAX]
   }
   .flatten()
   .view()
 
  SELECT_BEST_RECONSTRUCTION(data_ch)
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
