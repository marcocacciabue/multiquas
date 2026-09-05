workflow GET_BEST_RECONSTRUCTION {
  
  take:
    reconstruction_ch
  main:
    data_all_ch= reconstruction_ch
  .groupTuple().map{tup ->
      [tup[0]] + tup[1].flatten()
  }.map{ 
    row ->
      def a = row[0]
      def candidates = row.subList(1, row.size()) 
      def elementMAX = candidates.max { it.r_sq }
      [elementMAX]
  }
  .flatten()
  
  emit:
    data_all = data_all_ch .map { sample -> [sample.sample_id, sample] }
}