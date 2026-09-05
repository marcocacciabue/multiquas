include { GET_ALIGNMENT_S                                       } from '../subworkflows/align.nf'
include { GET_ALIGNMENT_M                                       } from '../subworkflows/align.nf'
include { GET_VARIANTS                                          } from '../subworkflows/align.nf'



workflow ALIGN {
  take:
    reads
    
  main:
    GET_ALIGNMENT_S(reads)
    
    single_bam_ch=GET_ALIGNMENT_S.out.single_bam.combine(GET_ALIGNMENT_S.out.reference_list.splitCsv(strip:true).map{tup ->
      [tup[0]] + tup[1].flatten()
      },by:0)


    GET_VARIANTS(GET_ALIGNMENT_S.out.single_bam)
    // Run aligment of the reads over all references available
    GET_ALIGNMENT_M(reads)

    multiple_bam_ch=GET_ALIGNMENT_M.out.multiple_bam.combine(GET_ALIGNMENT_M.out.reference_list.splitCsv(strip:true).map{tup ->
      [tup[0]] + tup[1].flatten()
      },by:0)

  emit:
    single_bam     =  single_bam_ch
    save           =  GET_ALIGNMENT_S.out.save
    variants       =  GET_VARIANTS.out.variants
    multiple_bam   =  multiple_bam_ch
}