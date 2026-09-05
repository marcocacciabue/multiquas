include { GET_ALL_RECONSTRUCTIONS                                       } from '../subworkflows/save.nf'
include { GET_SUMMARY                                       } from '../subworkflows/save.nf'
include { GET_BEST_RECONSTRUCTION                          } from '../subworkflows/get_best_reconstruction.nf'

workflow SAVE {
  take:
    reconstruction_ch
  reconstruction_ch_summary
  main:
    GET_ALL_RECONSTRUCTIONS(reconstruction_ch)
  GET_BEST_RECONSTRUCTION(reconstruction_ch)
  
  if (params.summary == "ON") {
    GET_SUMMARY(reconstruction_ch_summary)
  }
  emit:
    save_all  =   GET_ALL_RECONSTRUCTIONS.out.save_all
  data_all  =GET_BEST_RECONSTRUCTION.out.data_all
}