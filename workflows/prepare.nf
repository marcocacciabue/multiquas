include { BBDUK                                            } from '../modules/bbduk.nf'
include { EXTRACT_REFERENCE                                } from '../modules/extract_reference.nf'

workflow PREPARE {

  main:
  
  def input_csv = file("${params.input_csv}")
  //def input_ref = file("${params.input_ref}")


  
    // Create input channel from the contents of a CSV file
  read_ch = Channel.fromPath(input_csv)
  .splitCsv(header:true)
  .map { row -> [row.sample_id, file(row.fastq_1), file(row.fastq_2), file(row.input_ref)]}
  
  //TODO add check for reference file (genome size?)
  //ref_ch =Channel.fromPath(input_ref)
  //EXTRACT_REFERENCE(ref_ch)
  
  
  // Adapter trimming and post-trimming QC
  BBDUK(read_ch)
  


  
  emit:
  reads=BBDUK.out.trimmed_reads


}

