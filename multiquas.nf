#!/usr/bin/env nextflow

// Primary input
params.input_csv = "samples.csv"
params.input_ref = "data/referencia/referencia.fasta"
//params.input_ref2 = "data/referencia/referencia.fasta"

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { BBDUK } from './modules/bbduk.nf'
include { EXTRACT_REFERENCE } from './modules/extract_reference.nf'
include { MULTIPLE_ALIGNMENT } from './modules/multiple_aligment.nf'
include { SINGLE_ALIGNMENT } from './modules/single_alignment.nf'

include { SPLIT_ALIGNMENT } from './modules/split_alignment.nf'
include { MERGE_READS } from './modules/merge_reads.nf'
include { SPLIT_FILE } from './modules/split_file.nf'
include { COUNT } from './modules/count.nf'
include { RECONSTRUCT_QURE } from './modules/reconstruct_qure.nf'

include { MERGE_RECONSTRUCTED_QURE } from './modules/merge_reconstructed_qure.nf'
include { VARIANTS } from './modules/variants.nf'
include { HAPLOTYPE_EVALUATOR } from './modules/haplotype_evaluator.nf'

/*
 * Pipeline parameters
 */




workflow {

    // Create input channel from the contents of a CSV file
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, file(row.fastq_1), file(row.fastq_2)]}
        

    
    

    // Call processes
    //FASTQC(read_ch)
  // EXTRACT_REFERENCE(Channel.fromPath(params.input_ref))
   //EXTRACT_REFERENCE.out 
     //          .view()
    // Adapter trimming and post-trimming QC
   BBDUK(read_ch)

 //BBDUK.out.trimmed_reads.view()
 
   reads_trimmed_ch = BBDUK.out.combine(Channel.fromPath(params.input_ref))
   //reads_trimmed_ch_2 = BBDUK.out.combine(Channel.fromPath(params.input_ref))
  
   SINGLE_ALIGNMENT(reads_trimmed_ch)
   
   VARIANTS(SINGLE_ALIGNMENT.out.single_bam)
    // Run aligment of the reads over all references available
    MULTIPLE_ALIGNMENT(reads_trimmed_ch)
   // MULTIPLE_ALIGNMENT.out.multiple_bam.view()
    //MULTIPLE_ALIGNMENT.out.multiple_bam.view()
   split_reads_ch = MULTIPLE_ALIGNMENT.out.reference_list.first().splitCsv(strip:true)
	     .combine(MULTIPLE_ALIGNMENT.out.multiple_bam)

  
   SPLIT_ALIGNMENT(split_reads_ch)

  SPLIT_ALIGNMENT.out.splittted_reads
  MERGE_READS(SPLIT_ALIGNMENT.out.splittted_reads)
  //MERGE_READS1.out.view()
  
  //SPLIT_UNMAPPED(MULTIPLE_ALIGNMENT.out.multiple_bam)
  //SPLIT_UNMAPPED.out.view()
  //MERGE_READS2(SPLIT_UNMAPPED.out.merged_reads)
  //MERGE_READS2.out.view()
  RECONSTRUCT_QURE(MERGE_READS.out.reads)
  //RECONSTRUCT_QURE.out.haplotypes.groupTuple(by:2).view()
  //RECONSTRUCT_QURE2(MERGE_READS2.out.reads)

  MERGE_RECONSTRUCTED_QURE(RECONSTRUCT_QURE.out.haplotypes.groupTuple(by:2))
 // RECONSTRUCT_QURE.out.haplotypes.collectFile(name: "allhaplotypes.fasta").view()
  
  //concat(RECONSTRUCT_QURE2.out).view()
  
 haplo_ch=MERGE_RECONSTRUCTED_QURE.out
   .combine(VARIANTS.out,by:1)
  HAPLOTYPE_EVALUATOR(haplo_ch)
// new_ch.view()
}