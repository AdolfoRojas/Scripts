#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)
library(ggrepel)
my_comparisons <- list(c("Tumoral", "Normal"))
BioTypes <- c("lncRNA", "protein_coding") # "lncRNAs", 
subtype <- "all"
setwd("/media/run-projects/Adolfo/Datos_tesis/2do_Objetivo/analisis_muestras_y_DE")
gencode <- read.delim("../gencode_genes.v38.annotation.tab", sep = "\t")
gencode <- gencode[c("gene_id","gene_name", "gene_type")]
gencode <- unique(gencode)
BioType_ids <- gencode
biotype <- ""

de_inicial <- ""
de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")

colnames(de_inicial)[1] <- "Transcrito"
de <- merge(de_inicial, BioType_ids, by.x = "Transcrito", by.y = "gene_id")
de$diffexpressed <- "NO"
de$important <- "NO"
de$diffexpressed[de$log2FoldChange > 2 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$gene_type == "lncRNA" & de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -2 & de$padj < 0.05] <- "DOWN"
de$diffexpressed[de$gene_type == "lncRNA" & de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
de <- de[order(de$log2FoldChange, decreasing = T),]
de[1:3,]$important <- "UP"
de[(nrow(de)-2):nrow(de),]$important <- "DOWN"  
de$delabel <- NA
de$delabel[de$important != "NO"] <- de$gene_name[de$important != "NO"]
  
# plot adding up all layers we have seen so far
p1 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  #ggtitle(paste("Breast cancer diferential expression analysis\n","Gene level\n", subtype, " samples", sep = ""))+
  geom_point(size = 1.4) + 
  theme_minimal() +
  geom_text_repel(segment.size = 0.2, force = 8) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  #geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  theme(plot.title = element_text(hjust = 0, size=20, face="bold.italic", vjust= 0))  
#############################################################################################################################################
#                   miRNAs
#############################################################################################################################################
de <- ""
de <- read.delim(paste("Corregido_resultados_expresion_diferencial_miRNAs-gene_level/", subtype,"_samples_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  
colnames(de)[1] <- "Transcrito"
de$diffexpressed <- "NO"
de$important <- "NO"
de$diffexpressed[de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
de <- de[order(de$log2FoldChange, decreasing = T),]
de[1:3,]$important <- "UP"
de[(nrow(de)-2):nrow(de),]$important <- "DOWN"  
de$delabel <- NA
de$delabel[de$important != "NO"] <- de$Transcrito[de$important != "NO"]
de$delabel <- gsub("hsa-", "", de$delabel)
# plot adding up all layers we have seen so far
p2 <- ggplot(data = de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  #ggtitle(paste("Breast cancer diferential expression analysis\nmiRNAs\n", subtype, " samples", sep = ""))+
  geom_point(size = 1.4) + theme_minimal() +
  geom_text_repel(segment.size = 0.2, force = 8) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  #geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  theme(plot.title = element_text(hjust = 0, size=20, face="bold.italic", vjust= 0))
  
ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("Figura_volcano_RNA-miRNA.tiff", width = 16, height = 8, dpi = 300, units = "in")
system("ls -lh Figura_volcano_RNA-miRNA.tiff")