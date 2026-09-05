#!/usr/bin/env nextflow
params.version = "v0.0.1"

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
                              by Marco Cacciabue
                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     #############################################
                     samplesheet:        ${params.input_csv}
                     qure single:        ${params.qure_s}  
                     qure multiple:      ${params.qure_m}
                     clique single:      ON
                     clique multiple:    ${params.clique_m} 
                     haploflow single:   ${params.haploflow_s} 
                     haploflow multiple: ${params.haploflow_m} 
                     savage single:      ${params.savage_s} 
                     savage multiple:    ${params.savage_m} 
                     viquas single:      ${params.viquas_s} 
                     viquas multiple:    ${params.viquas_m} 

                     summary:            ${params.summary}  
                   
                     """.stripIndent()
}


// Primary input
params.input_csv = "samples.csv"


// Reconstructer selection. Clique_s is always ON.
params.qure_s    = "ON"
params.qure_m    = "ON"
params.clique_m  = "ON"
params.haploflow_s = "ON"
params.haploflow_m = "ON"
params.savage_s    = "ON"
params.savage_m    = "ON"
params.viquas_s    = "ON"
params.viquas_m    = "ON"


// Create summary?
params.summary   = "ON"

// Workflows INCLUDE statements

include { CHECK       } from './workflows/check.nf'
include { PREPARE     } from './workflows/prepare.nf'
include { SAVE        } from './workflows/save.nf'
include { ALIGN       } from './workflows/align.nf'
include { RECONSTRUCT } from './workflows/reconstruct.nf'



workflow {
  main:
  Logo()
  CHECK()
  PREPARE()
  //TODO include specific commands from parameters file to the reconstruction programs
  ALIGN(PREPARE.out.reads)

  RECONSTRUCT(ALIGN.out.single_bam,
              ALIGN.out.multiple_bam,
              ALIGN.out.variants)
 
  SAVE(RECONSTRUCT.out.reconstruction_ch,
      RECONSTRUCT.out.reconstruction_ch_summary)

  
  workflow.onComplete = { 
    if (workflow.success) { println "✓ Pipeline ran successfully!" } else 
    { println "❌ Pipeline failed! Error: ${workflow.errorMessage}" }
    println "Pipeline completed at: $workflow.complete"
    println "Duration    : ${workflow.duration}"
    println "workDir     : ${workflow.workDir}"
    println "commandLine : ${workflow.commandLine}"
    println "exit status : ${workflow.exitStatus}"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
  }
  
 publish: 
   ALIGNMENT_S_SAVE         = ALIGN.out.save
   SAVE_ALL_RECONSTRUCTIONS = SAVE.out.save_all
   SAVE_BEST_RECONSTRUCTION = SAVE.out.data_all
}


output {
    ALIGNMENT_S_SAVE         {
      path { sample_id, x -> "${sample_id}/ALIGNMENT_S" }
    }
    SAVE_ALL_RECONSTRUCTIONS {
      path { sample_id, 
             reconstructer,
             x  -> "${sample_id}/${reconstructer}"      }
    }
    SAVE_BEST_RECONSTRUCTION {
      path { sample_id, x -> "${sample_id}/BEST"        }
    }
}
