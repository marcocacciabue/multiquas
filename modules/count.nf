#!/usr/bin/env nextflow

process COUNT{
       tag "$sample_id"
    container "cacciabue/multiquas:developing"
  

    input:
    path(file)
    output:
    val ("*")

    cpus 4

    script:
    """
    wc -c $file

    """
}