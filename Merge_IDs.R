#!/usr/bin/env Rscript
setwd("../Tesis")
file1 <- "IDs_matrix.txt"
file2 <- "Human_lncRNAs.bed12"

IDs1 <- read.table(file1, header=F, sep = "|")
IDs2 <- read.table(file2, header=F, sep = "\t")

IDs1 <- IDs1[c(1,2,5,6,8)]
IDs2 <- IDs2[c(1:6)]
colnames(IDs1) <- c("Transcrip_ID","Gene_ID","Transcript","Gene","Type")
colnames(IDs2) <- c("chrom","chromStart","chromEnd","name","score","strand")

IDs1$Transcrip_ID <- strsplit(as.character(IDs1$Transcrip_ID), "\\.\\w*")
IDs2$name <- strsplit(as.character(IDs2$name), "\\.\\w*")

IDs1 <- transform(IDs1,Transcrip_ID=unlist(Transcrip_ID))
IDs2 <- transform(IDs2,name=unlist(name))

Merged_file <- as.data.frame(merge(IDs1, IDs2, by.x = "Transcrip_ID", by.y = "name"))

Merged_file_protein <- Merged_file[ which(Merged_file$Type=='protein_coding'), ]

#write.table(Merged_file_protein, sep = "\t",
#            file = "Merge_protein_lncRNAs.tsv", 
#            row.names = F, quote = F, col.names = T)
write.table(Merged_file, sep = "\t",
            file = "Merge_lncRNAs.tsv", 
            row.names = F, quote = F, col.names = T)
Merged_file <- read.table("Merge_lncRNAs.tsv", header=T, sep = "\t")


Replacing_file <- read.table(file1, header=F, sep = "|")


Replacing_file$V8 <- as.character(Replacing_file$V8) 
contador <- 0
contador2 <- 0
for (Ensembl_Transcript_ID in Replacing_file$V1){
  if (length(grep(substr(Ensembl_Transcript_ID, start = 1, stop = 15), Merged_file$Transcrip_ID))!=0){
    Replacing_file$V8[grep(substr(Ensembl_Transcript_ID, start = 1, stop = 15), Replacing_file$V1)] <- "lncRNA"
    contador2 <- contador2+1
  }
  contador <- contador +1
  porcentaje <- round(contador/length(Replacing_file$V1)*100, digits = 1)
  print(paste(porcentaje, " % completado, ",contador2, " elementos cambiados", sep = ""))
  }
Replacing_file$V8 <- as.factor(Replacing_file$V8)
write.table(Replacing_file, sep = "|",
            file = "lncRNAs_añadidos.txt", 
            row.names = F, quote = F, col.names = F)
