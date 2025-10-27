#!/usr/bin/env nextflow

process FASTQC {
 tag "$sampleId"
    container "cacciabue/multiquas:developing"
    publishDir "results/${sample_id}/fastqc", mode: 'copy'

    input:
    tuple path(read1), path(read2), path(sample_id)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${read1} ${read2}
    """
}