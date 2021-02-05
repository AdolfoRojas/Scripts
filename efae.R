#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)
library(ggrepel)
setwd("../Tesis/Expresion Diferencial")
my_comparisons <- list(c("Tumoral", "Normal"))
BioTypes <- c("lncRNA", "protein_coding") # "lncRNAs", 
subtipos <- c("all", "Her2", "LumB", "LumA", "Basal")
gencode <- read.delim("gencodeV35_fixed.txt", sep = "\t")
#gencode <- gencode[c("transcript_id","gene_name","transcript_name")]
gencode <- gencode[c("gene_id","gene_name", "gene_type")]
gencode <- unique(gencode)
for (biotype in BioTypes){
  #BioType_ids <- read.delim(paste(biotype, "_transcript_type_IDs_BED_Gencode.txt", sep = ""), sep = "\t")
  #BioType_ids <- merge(BioType_ids, gencode, by ="transcript_id")
  BioType_ids <- gencode[gencode$gene_type == biotype,]
  for (subtype in subtipos){
  de_inicial <- read.delim(paste("Resultados_expresion_diferencial_mRNAs/", subtype,"_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  colnames(de_inicial)[1] <- "Transcrito"
  de <- merge(de_inicial, BioType_ids, by.x = "Transcrito", by.y = "gene_id")
  de$diffexpressed <- "NO"
  de$important <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 2 & de$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2FoldChange < -2 & de$padj < 0.05] <- "DOWN"

  if(biotype == "lncRNAs") {
    de$important[de$log2FoldChange > 5 & de$padj < 0.05] <- "UP"
    de$important[de$log2FoldChange < -5 & de$padj < 0.05] <- "DOWN"
  }else {
    de$important[de$log2FoldChange > 7 & de$padj < 0.05] <- "UP"
    de$important[de$log2FoldChange < -6.75 & de$padj < 0.05] <- "DOWN"
  }
  
  de$delabel <- NA
  de$delabel[de$important != "NO"] <- de$gene_name[de$important != "NO"]
  
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    ggtitle(paste("Breast cancer diferential expression analysis\n", biotype, " - Gene level\n", subtype, " samples", sep = ""))+
    geom_point(size = 0.4) + 
    theme_minimal() +
    geom_text_repel(segment.size = 0.2, force = 8) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-2, 2), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    theme(plot.title = element_text(hjust = 0, size=20, face="bold.italic", vjust= 0))
  
  ggsave(paste("Volcano_plot_", subtype,"_samples-", biotype, ".png"), width = 16, height = 9, dpi = 500, units = "in")

  # Create test data.
  data <- data.frame(
    category=c("Upregulated", "Downregulated"),
    count=c(length(de[de$diffexpressed == "UP",]$diffexpressed), length(de[de$diffexpressed == "DOWN",]$diffexpressed)))
  # Compute percentages
  data$fraction <- data$count / sum(data$count)
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  # Compute a good label
  data$label <- paste0(data$category, "\n", data$count)
  # Make the plot
  ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    ggtitle(paste("Breast cancer diferential expression analysis\n", biotype, " - Gene level\n", subtype, " samples", sep = "")) +
    geom_rect(colour = "white") +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=6, colour = "white") +
    scale_fill_manual(values = c("blue", "red")) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size=20, face="bold.italic", vjust=-5))  
  ggsave(paste("Donut_plot_", subtype, "_gene_level_samples-", biotype, ".png", sep = ""), width = 16, height = 9, dpi = 500, units = "in")
  }
}