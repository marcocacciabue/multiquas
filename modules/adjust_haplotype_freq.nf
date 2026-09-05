#!/usr/bin/env nextflow

process ADJUST_HAPLOTYPE_FREQ {
  label 'process_low'
  container "cacciabue/multiquas:vortex.v0.0.1"  
  containerOptions "--cpus=${task.cpus}"
  tag "${[sample_id, contig_name].join(':')}"

  
  input:
   tuple val(sample_id),
         path(first_ref),
         path(general_ref),
         val(reconstructer),
         val(method),
         val(contig_name),
         path(stats),
         path(reconstructedVariants),
         path(reads_per_haplotype)
  
  
  output:
  tuple val("$sample_id"),
        path("$first_ref"),
        path("$general_ref"),
        val("$reconstructer"),
        val("${method}"),
        val("${contig_name}"),
        path("${reconstructedVariants}"),     
        path("${sample_id}_${reconstructer}_${method}_${contig_name}_adjusted_proportions.txt"), emit: haplotypes

  script:
    """
    #!/usr/bin/Rscript
     library(ape)
   ## TODO maybe this process could be managed by the filtering step. adjusting
   ## proportions in bash instead of R.
      
   stats<-read.table('${stats}')
   names(stats)<-c("references","length","mapped_reads","unmapped_reads")
   merged_reads<- stats\$mapped_reads+ stats\$unmapped_reads
   merged_reads_proportion<-round(merged_reads/sum(merged_reads),4)

   haplotypes<-read.FASTA('${reconstructedVariants}')

if ('${reconstructer}'=='viquas_s' | '${reconstructer}'=='viquas_m'){
   
   temp<-vector()        
   proportions<-vector()
   for (i in 1:length(haplotypes)) {
  
   temp<-strsplit(names(haplotypes[i]),split="_")[[1]]
  
   proportions[i]<-round(as.numeric(temp[2]),2)*100
  }  
 #values must be readjusted to accomodate for missing haplotypes (filtered during reconstruction)
 proportions = 100*(proportions/sum(proportions))
}




if ('${reconstructer}'=='savage_s' | '${reconstructer}'=='savage_m'){
   
   proportions<-vector()
   stats_haplo<-read.table('${reads_per_haplotype}')
   names(stats_haplo)<-c("references","length","mapped_reads","unmapped_reads")
   merged_reads_haplo<- stats_haplo\$mapped_reads+ stats_haplo\$unmapped_reads
   proportions<-round(merged_reads_haplo/sum(merged_reads_haplo),4)*100
 
}






if ('${reconstructer}'=='qure_s' | '${reconstructer}'=='qure_m'){
   temp<-vector()        
   proportions<-vector()
   for (i in 1:length(haplotypes)) {
  
   temp<-strsplit(names(haplotypes[i]),split="_")[[1]]
  
   proportions[i]<-round(as.numeric(temp[2]),2)
  }  
}

if ('${reconstructer}'=='haploflow_s' | '${reconstructer}'=='haploflow_m'){
   temp<-vector()        
   proportions<-vector()
  for (i in 1:length(haplotypes)) {
  
  temp<-strsplit(names(haplotypes[i]),split="_")[[1]]
  
  proportions[i]<-round(as.numeric(temp[4]),2)
  
  
}      
#haploflow produces an overall abundance value. 
total_abundance=sum(proportions)

proportions=100*proportions/total_abundance 

} 

if ('${reconstructer}'=='clique_s' | '${reconstructer}'=='clique_m'){

temp<-vector()        
proportions<-vector()
for (i in 1:length(haplotypes)) {
  
  temp<-strsplit(names(haplotypes[i]),split="_")[[1]]
  #change this index value to get the correct info
  proportions[i]<-round(as.numeric(temp[3]),2)*100
  }      
}


if ('${contig_name}'=="unmapped"){
proportions<-round(proportions*merged_reads_proportion[stats\$references=="*"],2)
  }

if ('${method}'=="single"){
#no adjusting, just rounding up.
proportions <- round(proportions,2)
} 


if ('${contig_name}'!="unmapped" & '${method}'=="multiple"){
proportions<-round(proportions*merged_reads_proportion[stats\$references=='${contig_name}'],2)

}



proportions<-as.character(proportions)
out_name = paste0('${sample_id}',"_",'${reconstructer}',"_",'${method}',"_",'${contig_name}',"_adjusted_proportions.txt")
writeLines(proportions, out_name)



    """
}
