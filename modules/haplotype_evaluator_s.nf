#!/usr/bin/env nextflow

process HAPLOTYPE_EVALUATOR_S {
   label 'process_low'
  tag "$sample_id"
  container "cacciabue/multiquas:developing"
  publishDir "results/${sample_id}/haplotypes_qure/evaluator_single", mode: 'copy'

  
  input:
    tuple val(sample_id), 
    path(haplotypes),
    path(proportions), 
    path(variants_vcf),
    val(contig_name)
    
  
  output:
     tuple val("${sample_id}"),
     path("${sample_id}_qure_haplotypes_aligned.fasta"),
     path("${proportions}"),
     path("${sample_id}_qure_graphs.png"),
     path("${sample_id}_R_sq.txt"),emit:reconstructed_data
     path("${sample_id}_single_qure_results.txt"),emit:results

 
  script:
    """
     #!/usr/bin/Rscript


suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
# FIRST sequence is consider as the reference!!

haplo_freq <- read.table('${proportions}')
names(haplo_freq) <- "freq"
query_alineado <- read.fasta('${sample_id}_qure_haplotypes_aligned.fasta',
                             seqtype = "DNA")
vcf <- readVcf('${variants_vcf}')

results_table <- data.frame(ID = names(query_alineado))
temp_aa <- vector()


start <- 1

end <- length(query_alineado[[1]])
for (i in start:end) {
  for (h in seq_along(query_alineado)) {
    temp_aa[h]<-query_alineado[[h]][i]
  }
  results_table <- cbind(results_table,temp_aa)
}
colnames(results_table) <- c("ID",seq(start,end,by=1))

temp2 <- results_table[2,]!=results_table[1,]
temp2 <- temp2 & (results_table[2,] != "-")

if (length(results_table[,1])>2){
for (S_N in 3:length(results_table[,1])) {
  temp <- (results_table[S_N,] != results_table[1,])
  temp <- temp & (results_table[S_N,] != "-")
  temp2 <- temp2 | temp
}

results_table_final <- results_table[,temp2]

temp2 <- results_table_final[2,] != results_table_final[1,]
for (S_N in 3:length(results_table_final[,1])) {
  temp <- (results_table_final[S_N,] != results_table_final[1,])
  temp2 <- rbind(temp2,temp)
  
}
}else{
  results_table_final <- results_table[,temp2]
  
  temp2 <- results_table_final[2,] != results_table_final[1,]
}

snp_haplotype <- temp2*haplo_freq\$freq[row(temp2)]

snp_haplotype <- colSums(snp_haplotype)/100

snp_haplotype <- snp_haplotype[-1]

haplotype_position <- as.numeric(names(snp_haplotype))

haplotype_freq <- as.vector(snp_haplotype)*100
Lofreq_freq <- (vcf@info\$AF)*100
lofreq_position <- start(ranges(vcf))

# combino vectores de posicion
Position <- (c(haplotype_position,lofreq_position))

#  saco duplicados que exceden

Position <- unique(Position)
# ordeno de menor a mayor
Position <- Position[order(Position)]

Haplo <- match(Position,haplotype_position)
Lofreq <- match(Position,lofreq_position)

data <- data.frame(position = Position,
lofreq = Lofreq_freq[Lofreq],
haplo = haplotype_freq[Haplo])


data[is.na(data)] <- 0
model <- lm(haplo ~ lofreq,
            data = data)
summary(model)

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
                    method = "single",
                    reconstructer = "QuRe")
                    


write.table(out_table,
'${sample_id}_single_qure_results.txt', 
sep=",",
row.names=FALSE)

write.table(R_sq,
'${sample_id}_R_sq.txt', 
sep=",",
col.names=FALSE,
row.names=FALSE)

"""
}
