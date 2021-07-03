#!/usr/bin/env Rscript
#############################################################################################################################################
#                   miRNAs
#############################################################################################################################################
library(ggplot2)
library(ggpubr)
library(ggrepel)
setwd("../Tesis/Expresion Diferencial")
my_comparisons <- list(c("Tumoral", "Normal"))
subtipos <- c("all", "Her2", "LumB", "LumA", "Basal")
for (subtype in subtipos){
  de <- read.delim(paste("Resultados_expresion_diferencial_miRNAs/", subtype,"_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  colnames(de)[1] <- "Transcrito"
  de$diffexpressed <- "NO"
  de$important <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
  
  de$important[de$log2FoldChange > 2.5 & de$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$important[de$log2FoldChange < -2.5 & de$padj < 0.05] <- "DOWN"
  de$delabel <- NA
  de$delabel[de$important != "NO"] <- de$Transcrito[de$important != "NO"]
  de$delabel <- gsub("hsa-", "", de$delabel)
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    ggtitle(paste("Breast cancer diferential expression analysis\nmiRNAs\n", subtype, " samples", sep = ""))+
    geom_point(size = 0.4) + 
    theme_minimal() +
    geom_text_repel(
      segment.size = 0.2,
      force = 8
    ) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    theme(plot.title = element_text(hjust = 0, size=20, face="bold.italic", vjust= 0))
  
  ggsave(paste("Volcano_plot_", subtype,"_samples-miRNAs.png"), width = 16, height = 9, dpi = 500, units = "in")
  
  de <- read.delim(paste("Resultados_expresion_diferencial_miRNAs/", subtype,"_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  colnames(de)[1] <- "Transcrito"
  de$diffexpressed <- "NO"
  de$important <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 1 & de$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2FoldChange < -1 & de$pvalue < 0.05] <- "DOWN"
  # Create test data.
  data <- data.frame(
    category=c("Upregulated", "Downregulated"),
    count=c(length(de[de$diffexpressed == "UP",]$diffexpressed), length(de[de$diffexpressed == "DOWN",]$diffexpressed))
  )
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
    ggtitle(paste("Breast cancer diferential expression analysis\nmiRNAs - Gene level\n", subtype, " samples", sep = "")) +
    geom_rect(colour = "white") +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=6, colour = "white") +
    scale_fill_manual(values = c("blue", "red")) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size=20, face="bold.italic", vjust=-5))
  ggsave(paste("Donut_plot_", subtype, "_samples-miRNAs.png", sep = ""), width = 16, height = 9, dpi = 500, units = "in")
  
  mdp <- read.delim(paste("Resultados_expresion_diferencial_miRNAs/", subtype, "_all_controls_sample_scores.tsv", sep = ""), sep = " ")
  
  ggdensity(mdp, x = "allgenes.Score",
            add = "mean", rug = TRUE,
            color = "allgenes.Class", fill = "allgenes.Class",
            palette = c("#00AFBB", "#E7B800"))+
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples - miRNAs", sep = ""))+ 
    theme(legend.position = c(0.75,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), 
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.direction = "horizontal")+
    ylab("Density")+ xlab("Sample score - all") + labs(fill = "Sample type", color ="Sample type")
  ggsave(paste("Density_plot_", subtype,"_samples_distribtion-miRNAs.png"), width = 16, height = 9, dpi = 500, units = "in")
  
  
  ggbarplot(mdp, x = "allgenes.Sample", y = "allgenes.Score",
            fill = "allgenes.Class",               # change fill color by cyl
            color = "allgenes.Class",            # Set bar border colors to white
            palette = "jco",            # jco journal color palett. see ?ggpar
            sort.val = "asc",           # Sort the value in dscending order
            sort.by.groups = TRUE,      # Sort inside each group
            x.text.angle = 90           # Rotate vertically x axis texts
  )+ theme(axis.text.x=element_blank(), legend.position = c(0.15,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), 
           legend.title = element_text(size = 15),
           legend.text = element_text(size = 12),
           legend.direction = "horizontal")+
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples - miRNAs", sep = ""))+
    ylab("All genes score")+ xlab("Samples")+ labs(fill = "Sample type", color = "Sample type")
  ggsave(paste("Bar_plot_", subtype,"_samples_distribtion-miRNAs.png"), width = 16, height = 9, dpi = 500, units = "in")
  
  
  stat.test <- compare_means(allgenes.Score ~ allgenes.Class, mdp, comparisons = my_comparisons)
  stat.test$p<- formatC(stat.test$p, format = "e", digits = 2)
  ggboxplot(mdp, x = "allgenes.Class", y = "allgenes.Score",
            color = "allgenes.Class", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
            add = "jitter", shape = "allgenes.Class") +
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples - miRNAs", sep = ""))+
    theme(legend.position = "none")+
    stat_pvalue_manual(stat.test, label = "{p.signif}\np ={p} {method} test", y.position = max(mdp$allgenes.Score)+0.5, vjust = 0.8)+
    ylab("All genes score")+ xlab("Sample type")
  ggsave(paste("Scatter_plot_", subtype,"_samples-miRNAs.png"), width = 16, height = 9, dpi = 500, units = "in")
}