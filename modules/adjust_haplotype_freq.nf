#!/usr/bin/env nextflow

process ADJUST_HAPLOTYPE_FREQ {
  label 'process_low'
  tag "$sample_id"
  container "cacciabue/multiquas:developing"  
  containerOptions "--cpus=${task.cpus}"

  
  input:
   tuple path(reconstructedVariants),     
      val(sample_id),
      val(contig_name),
      path(stats)
  
  
  output:
  tuple    val("$sample_id"),
    path("${contig_name}_temporary_reconstructedVariants.txt"),     
              path("${contig_name}_adjusted_proportions.txt"),
              val("${contig_name}"), emit: haplotypes

  script:
    """
    #!/usr/bin/Rscript
     library(ape)

      
   stats<-read.table('${stats}')
   names(stats)<-c("references","length","mapped_reads","unmapped_reads")
   merged_reads<- stats\$mapped_reads+ stats\$unmapped_reads
   merged_reads_proportion<-round(merged_reads/sum(merged_reads),4)

   haplotypes<-read.FASTA('${reconstructedVariants}')

   temp<-vector()        
   proportions<-vector()
   for (i in 1:length(haplotypes)) {
  
   temp<-strsplit(names(haplotypes[i]),split="_")[[1]]
  
   proportions[i]<-round(as.numeric(temp[2]),2)
  
  
}      

if ('${contig_name}'=="unmapped"){
proportions<-round(proportions*merged_reads_proportion[stats\$references=="*"],2)
  }

if ('${contig_name}'=="single"){
#no adjusting
proportions <- proportions 
} 


if ('${contig_name}'!="unmapped" & '${contig_name}'!="single"){
proportions<-round(proportions*merged_reads_proportion[stats\$references=='${contig_name}'],2)

}



proportions<-as.character(proportions)
out_name = paste0('${contig_name}',"_adjusted_proportions.txt")
writeLines(proportions, out_name)


    """
}