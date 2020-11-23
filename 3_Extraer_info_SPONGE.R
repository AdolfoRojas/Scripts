#!/usr/bin/env Rscript
library(pdftools)
library(limma)
library(edgeR)
load("entorno_completo_sponge2.RData")

pdf("df_plot1.pdf", height = 9, width = 16)
boxplot(norm_RNA$counts[,1:100], las = 2)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 500)
png::writePNG(bitmap, "norm_rna_boxplot.png")
pdf("df_plot1.pdf", height = 9, width = 16)
boxplot(norm_mir$counts[,1:100], las = 2)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 500)
png::writePNG(bitmap, "norm_mir_boxplot.png")

gencode <- read.delim("gencodeV35_fixed.txt")
ids <- as.data.frame(unique(append(ceRNA_interactions_fdr$geneA,ceRNA_interactions_fdr$geneB)))
colnames(ids) <- "gene_id"
gencode <- merge(gencode, ids, by = "gene_id")
levels(as.factor(gencode$gene_type))

gene_type_protein_coding <- gencode[gencode$gene_type=="protein_coding",]
levels(as.factor(gene_type_protein_coding$transcript_type))

gencode <- unique(gencode[c("gene_id", "gene_type", "gene_name")])

gene_type_lncRNA <- gencode[gencode$gene_type=="lncRNA",]
gene_type_lncRNA <- gene_type_lncRNA["gene_id"]
DE_info <- read.delim("../Resultados_expresion_diferencial_mRNAs-gene_level/all_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
DE_info <- DE_info[DE_info$pvalue < 0.05,]
DE_info <- DE_info[DE_info$log2FoldChange < -1 | DE_info$log2FoldChange > 1,]

colnames(DE_info)[1] <- "gene_id"

DE_lncRNA <- merge(gene_type_lncRNA,DE_info, by = "gene_id")
DE_lncRNA_ids <- as.data.frame(DE_lncRNA$gene_id)
DE_lncRNA_up <- DE_lncRNA[DE_lncRNA$log2FoldChange < 1,]
DE_lncRNA_down <- DE_lncRNA[DE_lncRNA$log2FoldChange > -1,]
colnames(DE_lncRNA_ids) <- "gene_id"
new_table <- merge(gencode, ceRNA_interactions_fdr, by.x = "gene_id", by.y = "geneB")
colnames(new_table)[1:3] <- c("GeneB", "gene_typeB", "gene_nameB")
new_table <- merge(gencode, new_table, by.x = "gene_id", by.y = "geneA")
colnames(new_table)[1:3] <- c("GeneA", "gene_typeA", "gene_nameA")
head(new_table)
lncRNA_ceRNA_interactions <- new_table[grepl("lncRNA", new_table$gene_typeA) | grepl("lncRNA", new_table$gene_typeB),]
A <- merge(DE_lncRNA_ids, lncRNA_ceRNA_interactions, by.x = "gene_id", by.y = "GeneA")
colnames(A)[1] <- "GeneA"
B <- merge(DE_lncRNA_ids, lncRNA_ceRNA_interactions, by.x = "gene_id", by.y = "GeneB")
colnames(B)[1] <- "GeneB"
DE_lncRNA_ceRNA_interactions <- rbind(A,B) #####
save(DE_lncRNA_ceRNA_interactions, file = "DE_lncRNA_ceRNA_interactions.Rdata")
save(DE_lncRNA, file = "DE_lncRNA.Rdata")
save(genes_miRNA_candidates, file = "genes_miRNA_candidates.Rdata")
mirRNAs_unicos <- sapply(genes_miRNA_candidates, function(x){x[1]})
mirRNAs_unicos <- unlist(mirRNAs_unicos, use.names = F)
mirRNAs_unicos <- unique(mirRNAs_unicos)
length(mirRNAs_unicos)