#!/usr/bin/env Rscript
library(foreach)
print("Iniciando")
freq <- read.delim2("SABE1171.hg38.VQSR-AS.biallelic.recodecleanFBFD.Anno3.afreq", header = TRUE,  sep = "\t")
print("Archivo de frecuencias cargado")
length(freq$ALT_FREQS)
freq <- freq[as.numeric(freq$ALT_FREQS) >= 0.01,]
print("Filtro por MAF terminado")
length(freq$ALT_FREQS)
#Map <- read.delim2("../SABE1171.hg38.VQSR-AS.biallelic.recodecleanFBFD.Anno3.bim", header = F,  sep = "\t")
Map <- read.delim2("SABE1171.hg38.VQSR-AS.biallelic.recodecleanFBFD.Anno3.bim", header = F,  sep = "\t")
print("Archivo de coordenadas cargado")
Map <- Map[c(2,4)]
variant_info <- merge(freq, Map, by.y = "V2", by.x = "ID")
names(variant_info)[2] <- "CHROM"
names(variant_info)
variant_info$CHROM <- paste("chr",variant_info$CHROM, sep = "")
remove(Map)
remove(freq)
#Ref <- read.delim2("Human_lncRNAs.bed12", header = F,  sep = "\t")
Ref <- read.delim2("../Tesis/Human_lncRNAs.bed12", header = F,  sep = "\t")
print("Referencia de lncRNAs cargada")
summary_df <- data.frame(lncRNA_ID=character(),
                Num_Var=character(),
                Variants=character(), 
                stringsAsFactors=FALSE)
full_list <- data.frame(lncRNA_ID=character(),                
                Variants=character(), 
                stringsAsFactors=FALSE)
for (lncRNA in 1:length(Ref$V4)) {
    lncRNA_variants <- variant_info[variant_info$CHROM == Ref$V1[lncRNA] & variant_info$V4 >= Ref$V2[lncRNA] &variant_info$V4 <= Ref$V3[lncRNA],]
    if (length(lncRNA_variants$ID != 0)) {
       masisi <- paste(lncRNA_variants$ID,"(",lncRNA_variants$ALT_FREQS,")", sep = "")
       variants_number <- length(masisi)       
       masisi <- paste(masisi, collapse = " ")               
       Variants2 <- lncRNA_variants[c("ID","ALT_FREQS")]
       summary_df <- rbind(summary_df, data.frame(Ref$V4[lncRNA], variants_number, masisi, fix.empty.names = F))
       full_list <- rbind(full_list, data.frame(Ref$V4[lncRNA],Variants2, fix.empty.names = F)) }
    }

colnames(full_list) <- c("lncRNA_ID","Variant", "Freq")
colnames(summary_df) <- c("lncRNA_ID","Number_Variants", "Variants(Freq)")
print("Escribiendo")

#write.table(summary_df, sep = "\t", file = "lncRNA_var_SABE1171_summary.tab", row.names = F, quote = F, col.names = T)
View(summary_df)
remove(summary_df)
#write.table(full_list, sep = "\t", file = "lncRNA_var_SABE1171_full_list.tab", row.names = F, quote = F, col.names = T)
View(full_list)
remove(full_list)

#--------------------------------------------------------------------------------------------------------------------------------#

out_file <- read.delim2("../Trabajillos/lncRNAs/no-nan_lncRNA_var_SABE1171_full_list2.tab", sep = "\t", header = T)
out_file$lncRNA_ID <- as.factor(out_file$lncRNA_ID)
total <- length(levels(out_file$lncRNA_ID))

summary_df2 <- data.frame(lncRNA_ID=character(),
                Num_Var=character(),
                Variants=character(), 
                stringsAsFactors=FALSE)

for (lncRNA_num in 1:length(levels(out_file$lncRNA_ID))) {
    print(paste(lncRNA_num, " ",lncRNA_num/total*100, "% completado", sep = ""))
    lncRNA <- levels(out_file$lncRNA_ID)[lncRNA_num]
    lncRNA_variants <- out_file[out_file$lncRNA_ID == lncRNA ,]
    if (length(lncRNA_variants$Variant != 0)) {
       masisi <- paste(lncRNA_variants$Variant,"(",lncRNA_variants$Freq,")", sep = "")
       variants_number <- length(masisi)       
       masisi <- paste(masisi, collapse = " ")      
       summary_df2 <- rbind(summary_df2, data.frame(lncRNA, variants_number, masisi, fix.empty.names = F))}
    }
colnames(summary_df2) <- c("lncRNA_ID","Number_Variants", "Variants(Freq)")
write.table(summary_df2, sep = "\t", file = "../Trabajillos/lncRNAs/no-nan_lncRNA_var_SABE1171_summary.tab", row.names = F, quote = F, col.names = T)