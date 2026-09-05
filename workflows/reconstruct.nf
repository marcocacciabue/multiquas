include { GET_QURE_S    } from '../subworkflows/get_qure_s.nf'
include { GET_QURE_M    } from '../subworkflows/get_qure_m.nf'
include { GET_CLIQUE_S  } from '../subworkflows/get_clique_s.nf'
include { GET_CLIQUE_M  } from '../subworkflows/get_clique_m.nf'
include { GET_HAPLOFLOW_S  } from '../subworkflows/get_haploflow_s.nf'
include { GET_HAPLOFLOW_M  } from '../subworkflows/get_haploflow_m.nf'
include { GET_SAVAGE_S  } from '../subworkflows/get_savage_s.nf'
include { GET_SAVAGE_M  } from '../subworkflows/get_savage_m.nf'
include { GET_VIQUAS_S } from '../subworkflows/get_viquas_s.nf'
include { GET_VIQUAS_M } from '../subworkflows/get_viquas_m.nf'


workflow RECONSTRUCT {
  take:
    single_bam
    multiple_bam
    variants

  
  main:
    GET_CLIQUE_S(single_bam,
                 variants)
    reconstruction_ch=GET_CLIQUE_S.out.reconstructed_data
    reconstruction_ch_summary= GET_CLIQUE_S.out.results
    
         if (params.viquas_s == "ON"){
    GET_VIQUAS_S(single_bam,
               variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_VIQUAS_S.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_VIQUAS_S.out.results)
}
     if (params.haploflow_s == "ON"){
    GET_HAPLOFLOW_S(single_bam,
               variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_HAPLOFLOW_S.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_HAPLOFLOW_S.out.results)
}
     if (params.savage_s == "ON"){
    GET_SAVAGE_S(single_bam,
               variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_SAVAGE_S.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_SAVAGE_S.out.results)
}
  
  if (params.qure_s == "ON"){
    GET_QURE_S(single_bam,
               variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_QURE_S.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_QURE_S.out.results)
}
  if (params.savage_m == "ON"){
    GET_SAVAGE_M(multiple_bam,
               variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_SAVAGE_M.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_SAVAGE_M.out.results)
  }
  if (params.qure_m == "ON"){
    GET_QURE_M(multiple_bam,
               variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_QURE_M.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_QURE_M.out.results)
  }
  
  if (params.clique_m == "ON"){
    GET_CLIQUE_M(multiple_bam,
                 variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_CLIQUE_M.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_CLIQUE_M.out.results)
  }
  
    if (params.haploflow_m == "ON"){
    GET_HAPLOFLOW_M(multiple_bam,
                 variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_HAPLOFLOW_M.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_HAPLOFLOW_M.out.results)
  }
     if (params.viquas_m == "ON"){
    GET_VIQUAS_M(multiple_bam,
                 variants)
    reconstruction_ch=reconstruction_ch
    .mix(GET_VIQUAS_M.out.reconstructed_data)
    reconstruction_ch_summary=reconstruction_ch_summary
    .mix(GET_VIQUAS_M.out.results)
  }
  emit:
    reconstruction_ch         = reconstruction_ch
    reconstruction_ch_summary = reconstruction_ch_summary
}
