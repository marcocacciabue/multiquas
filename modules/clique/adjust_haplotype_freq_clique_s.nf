#!/usr/bin/env nextflow

process ADJUST_HAPLOTYPE_FREQ_CLIQUE_S {
  label 'process_low'
  tag "$sample_id"
  container "cacciabue/multiquas:developing"  

  containerOptions "--cpus=${task.cpus}"

  
  input:
   tuple val(sample_id),
         val(reconstructer),
         path(reconstructedVariants), 
         val(contig_name),
         path(stats)
  
  
  output:
   tuple val("$sample_id"),
         val("$reconstructer"),
         path("${reconstructedVariants}"),     
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
  #change this index value to get the correct info
  proportions[i]<-round(as.numeric(temp[3]),2)
  
  
}      



  #no adjusting proportions

proportions<-as.character(proportions*100)
out_name = paste0('${contig_name}',"_adjusted_proportions.txt")
writeLines(proportions, out_name)

    """
}