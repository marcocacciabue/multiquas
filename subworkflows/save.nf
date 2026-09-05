workflow GET_ALL_RECONSTRUCTIONS {
  
  take:
    reconstruction_ch
  main:
    save_all_ch =reconstruction_ch.map { sample_id , sample-> [sample_id, sample.reconstructer,
                                                               sample]
    }
  
  emit:
    save_all = save_all_ch
}

workflow GET_SUMMARY {
  
  take:
    reconstruction_ch_summary
  main:
    reconstruction_ch_summary
  .collectFile(name: 'multiquas_results_summary.txt', 
               newLine: false,
               keepHeader:true,
               storeDir:'results')
  .subscribe { file ->
      println "Summary is saved to file: $file"
    println "${file.text}"
  }
}
