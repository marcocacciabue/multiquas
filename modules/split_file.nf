#!/usr/bin/env nextflow

process SPLIT_FILE{
       tag "$sample_id"
    container "cacciabue/multiquas:developing"
  

    input:
    val(file)
    output:
    path("*contig.txt")

    cpus 4

    script:
    """


    """
}