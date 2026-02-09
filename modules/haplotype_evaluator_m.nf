#!/usr/bin/env nextflow

process HAPLOTYPE_EVALUATOR_M {
  label 'process_low'
  tag "${sample_id}"
  container "cacciabue/multiquas:developing"
  publishDir "results/${sample_id}/haplotypes_qure/evaluator_multiple", mode: 'copy'

  input:
  tuple val(sample_id), path(haplotypes), path(proportions), path(variants_vcf), val(contig_name)

  output:
  tuple val("${sample_id}"), path("${sample_id}_qure_haplotypes_aligned.fasta"), path("${proportions}"), path("${sample_id}_qure_graphs.png"), path("${sample_id}_R_sq.txt"), emit: reconstructed_data
  path ("${sample_id}_multiple_qure_results.txt"), emit: results

  script:
  """
     #!/usr/bin/Rscript

#####    Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))


#####    Load helper functions
#takes alignment and returns a dataframe with the SNPs 
generateSnpFromAlignment <- function(alignment){
  
  snps <- data.frame(ID = names(alignment))
  temp_snps <- vector()
  start <- 1
  end <- length(alignment[[1]])
  for (position in start:end) {
    for (entry in seq_along(alignment)) {
      temp_snps[entry]<-alignment[[entry]][position]
    }
    snps <- cbind(snps,temp_snps)
  }
  
  colnames(snps) <- c("ID",seq(start,end,by=1))
  
  #get positions with changes between reference and 
  #first sequence entry. Discard gaps.
  keepPosition <- (snps[2,]!=snps[1,]) & (snps[2,] != "-")
  
  #repeat process for the rest of the sequence entries (if more than one).
  if (length(snps[,1])>2){
    for (S_N in 3:length(snps[,1])) {
      keepTemp <- (snps[S_N,] != snps[1,]) & (snps[S_N,] != "-")
      keepPosition <- keepPosition | keepTemp
    }
  }  
  #if none of the positions meets the criteria
  #a dummy position is selected (first).
  #TODO: create a more elegant solution.
  if(sum(keepPosition)==1){
    
    snpsSelected <- snps[,1:2]
    
  }else{
    snpsSelected <- snps[,keepPosition]
  }
  
  keepPositionAll <- snpsSelected[2,] != snpsSelected[1,]
  
  #get from the dataframe only if difference from reference is met
  
  for (S_N in 3:length(snpsSelected[,1])) {
    temp <- (snpsSelected[S_N,] != snpsSelected[1,])
    keepPositionAll <- rbind(keepPositionAll,temp)
  }
  
  return(list(snps=snpsSelected,
              snpsLogical=keepPositionAll))
}

getComparisonData <- function(snpsSelected,
                              haplo_freq,
                              vcf){
  snp_haplotype <- snpsSelected\$snpsLogical*haplo_freq\$freq[row(snpsSelected\$snpsLogical)]
  snp_haplotype <- colSums(snp_haplotype)/100
  
  #remove reference
  snp_haplotype <- snp_haplotype[-1]
  
  haplotype_position <- as.numeric(names(snp_haplotype))
  
  haplotype_freq <- as.vector(snp_haplotype)*100
  Lofreq_freq <- (vcf@info\$AF)*100
  lofreq_position <- start(ranges(vcf))
  
  # combine position vectors, remove duplicates and sort
  Position <- (c(haplotype_position,lofreq_position))
  Position <- unique(Position)
  Position <- Position[order(Position)]
  
  #keep only position present in both dataset.
  Haplo <- match(Position,haplotype_position)
  Lofreq <- match(Position,lofreq_position)
  
  data <- data.frame(position = Position,
                     lofreq = Lofreq_freq[Lofreq],
                     haplo = haplotype_freq[Haplo])
  
  
  data[is.na(data)] <- 0
  return(data)
}
# FIRST sequence is consider as the reference!!

haplo_freq <- read.table('${proportions}')
names(haplo_freq) <- "freq"
alignment <- read.fasta('${sample_id}_qure_haplotypes_aligned.fasta',
                             seqtype = "DNA")
vcf <- readVcf('${variants_vcf}')
snpsSelected <- generateSnpFromAlignment(alignment)
data<-getComparisonData(snpsSelected,
                  haplo_freq,
                  vcf)

model <- lm(haplo ~ lofreq,
            data = data)

R_sq <- round(summary(model)\$r.squared,3)


eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(        a = format(unname(coef(model))[1], digits = 4),
                              b = format(unname(coef(model))[2], digits = 4),
                              r2 = format(summary(model)\$r.squared, digits = 3)))

dftext <- data.frame(x = 20, y = 60, eq = as.character(as.expression(eq)))

# plots graph
p<-ggplot(data,aes(x=lofreq,
                   y=haplo,
                   label=position))+
  stat_smooth(method="lm")+
  geom_point()+
  geom_label_repel(max.overlaps = 20)+
  ylab("Reconstructed SNVs (%)")+
  xlab("Lofreq SNVs (%)")+
  theme(axis.title=element_text(size=14))+
  geom_text(mapping = aes(x=x, y=y, label = eq), size=7, data = dftext, parse = TRUE)

ggsave('${sample_id}_qure_graphs.png',
       p,
       units="in", 
       width=10, 
       height=8,
       dpi=600)

out_table <- data.frame(sample_id = '${sample_id}',
                        rsquared = R_sq,
                        method = "multiple",
                        reconstructer = "QuRe")

write.table(out_table,
            '${sample_id}_multiple_qure_results.txt', 
            sep=",",
            row.names=FALSE)

write.table(R_sq,
            '${sample_id}_R_sq.txt', 
            sep=",",
            col.names=FALSE,
            row.names=FALSE)
"""
}
