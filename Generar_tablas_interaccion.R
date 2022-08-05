#!/usr/bin/env Rscript

genesymbol_to_geneid <- read.delim("gencode_genes.v38.annotation.tab")
genesymbol_to_geneid <- genesymbol_to_geneid[c("gene_name","gene_id")]

## PPI
protein_int <- read.delim("GTEx_TFs_PPI/2_interacciones_sobre_700_gene_symbol_in_GTEx.tab")
names(protein_int) <- c("element1","element2")
protein_int <- protein_int[!duplicated(protein_int),]

protein_int2 <- protein_int
names(protein_int2) <- c("element2","element1")
protein_int <- rbind(protein_int,protein_int2)
protein_int$int_type <- "PPI"
protein_int <- protein_int[!duplicated(protein_int),]

## TFs
tfs_int <- read.delim("GTEx_TFs_PPI/5_Interacciones_TF-Target_tejido_mamario.tsv")
tfs_int <- tfs_int[c("tf","target")]
names(tfs_int) <- c("element1","element2")
tfs_int <- tfs_int[!duplicated(tfs_int),] # 21666     2
tfs_int <- tfs_int[!duplicated(tfs_int),]
tfs_int$int_type <- "TF-Target"

## NPInter

NPInter_int <- read.delim("lncRNAfunc_data/lncRNAfunc_data_expresado.tab", sep = "\t",header=F)
names(NPInter_int)[1:2] <- c("element1","element2")
NPInter_int$V3 <- NULL
NPInter_int$V4 <- NULL
NPInter_int <- NPInter_int[!duplicated(NPInter_int),] 

NPInter_int2 <- NPInter_int
names(NPInter_int2) <- c("element2","element1")
NPInter_int <- rbind(NPInter_int,NPInter_int2)
NPInter_int$int_type <- "lncRNAfunc"
NPInter_int <- NPInter_int[!duplicated(NPInter_int),]


## ceRNA
load("Sponge/1_entorno_completo_sponge2.RData")
ceRNA_int <- ceRNA_interactions_fdr[c("geneA","geneB")]
names(ceRNA_int) <- c("element2","element1")
ceRNA_int <- ceRNA_int[!duplicated(ceRNA_int),]

ceRNA_int2 <- ceRNA_int
names(ceRNA_int2) <- c("element1","element2")
ceRNA_int <- rbind(ceRNA_int2,ceRNA_int)
ceRNA_int$int_type <- "ceRNA"
ceRNA_int <- ceRNA_int[!duplicated(ceRNA_int),]

genes_miRNA_candidates2 <- NULL
head(genes_miRNA_candidates)
for (interactor in names(genes_miRNA_candidates)){
   if (length(genes_miRNA_candidates[[interactor]][,1]) != 0){
      for(i in 1:length(genes_miRNA_candidates[[interactor]][,1])){
         genes_miRNA_candidates2 <- rbind(genes_miRNA_candidates2,as.data.frame(cbind(interactor, genes_miRNA_candidates[[interactor]][i,])))}}}
genes_miRNA_candidates2$coefficient <- NULL
names(genes_miRNA_candidates2) <- c("element2","element1")
genes_miRNA_candidates2 <- genes_miRNA_candidates2[!duplicated(genes_miRNA_candidates2),]
genes_miRNA_candidates2$int_type <- "miRNA-mRNA"
genes_miRNA_candidates2 <- genes_miRNA_candidates2[!duplicated(genes_miRNA_candidates2),]

final_int <- rbind(NPInter_int,ceRNA_int, tfs_int, protein_int, genes_miRNA_candidates2)
write.table(final_int, sep = "\t", file = "final_interaction.tab", row.names = F, quote = F, col.names = T)