#!/usr/bin/env Rscript
library(miRBaseConverter)
matrix <- read.delim("matrix_miRNAs_isoforms.tab", sep ='\t')
Accessions=matrix$X
result2=miRNA_AccessionToName(Accessions,targetVersion="v22")
result2 = result2[complete.cases(result2), ]
names(matrix)[1] <- "Accession"
matrix2 <- merge(result2, matrix, by = "Accession")
length(matrix$Accession)
length(matrix2$TargetName)
length(result2$TargetName)
names(matrix2) <- gsub("\\.", "-", names(matrix2))
write.table(matrix2, file = "matrix_miRNAs_isoforms_with_names.tab", sep ='\t', row.names = F, quote = F)