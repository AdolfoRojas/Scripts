setwd("../Trabajillos")

library(mogene20sttranscriptcluster.db)
Annot <- data.frame(ID = rep("---", times = 41345),
                    SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=","), 
                    DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=","),
                    ENTREZID=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=","),
                    ENSEMBLID=sapply(contents(mogene20sttranscriptclusterENSEMBL), paste, collapse=","))

Annot$ID <- rownames(Annot)
Annot$ID <- as.factor(Annot$ID)
affi2 <- Annot[c(1,2)]
levels(affi2$SYMBOL)[levels(affi2$SYMBOL) == "NA"] <- "---"

Files <- c("GSE119257 M0 vs M1", "GSE119257 M2 vs M1") 
for (File in Files) {
  Original <- read.delim(paste(File, ".txt", sep = ""), header=T, sep = "\t")
  Original <- Original[c(1:8)]
  print(length(Original$ID))
  colnames(affi2) <- c("ID", "Gene_ID")
  archivo_final <- as.data.frame(merge(Original, affi2, by = "ID"))
  print(length(archivo_final$ID))
  write.table(archivo_final, sep = "\t", file = paste(File, "_Gene_ID.tsv", sep = ""), row.names = F, quote = F, col.names = T)
  
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------#


original2 <- read.delim("GPL16570-1802.txt", header=T, sep = "\t")
original2_recortado <- original2[c(1,9:10)]
Original <- read.delim("GSE119257 M0 vs M1.txt", header=T, sep = "\t")
original_recortado <- Original[1]
merged <- as.data.frame(merge(original_recortado, original2_recortado, by = "ID"))

merged$gene_assignment <- as.character(merged$gene_assignment)
merged$mrna_assignment <- as.character(merged$mrna_assignment)
test <- merged[18210:18230,]

c <- data.frame(ID=character(),
                Gene_Symbol=character(), 
                stringsAsFactors=FALSE) 

count = 0
for (Row in 1:length(merged$ID)) {
  Symbol <- strsplit(merged$gene_assignment[Row], " // ")[[1]][2]
  ID <- merged$ID[Row]
  if (is.na(Symbol)==T){
    Symbol <- strsplit(merged$mrna_assignment[Row], " // ")[[1]][3]
  }
  c = rbind(c, data.frame(ID,Symbol, fix.empty.names = F))
  count = count+1
  print(count)
}

ya_reali <- read.delim("GSE119257_M0_vs_M1_Gene_ID.tsv", sep = "\t", header = T)