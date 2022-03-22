#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)
library(ggrepel)
my_comparisons <- list(c("Tumoral", "Normal"))
BioTypes <- c("lncRNA", "protein_coding") # "lncRNAs", 
subtipos <- c("all", "Her2", "LumB", "LumA", "Basal")
gencode <- read.delim("../gencode_genes.v38.annotation.tab", sep = "\t")
gencode <- gencode[c("gene_id","gene_name", "gene_type")]
gencode <- unique(gencode)
mdp <- ""
for (biotype in BioTypes){  
  BioType_ids <- gencode[gencode$gene_type == biotype,]
  for (subtype in subtipos){
  de_inicial <- ""
  if (subtype != "all") {
     de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  } else {
     de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  }
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
    de$important[de$log2FoldChange < -6.75 & de$padj < 0.05] <- "DOWN"}
  
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
  
  ggsave(paste("Volcano_plot_", subtype,"_samples-", biotype, ".png", sep =""), width = 16, height = 9, dpi = 500, units = "in")}}
BioType_ids <- gencode
biotype <- ""

for (subtype in subtipos){
  de_inicial <- ""
  if (subtype != "all") {
     de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  } else {
     de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  }
  colnames(de_inicial)[1] <- "Transcrito"
  de <- merge(de_inicial, BioType_ids, by.x = "Transcrito", by.y = "gene_id")
  de <- de[de$padj < 0.05,]
  de <- de[de$log2FoldChange < -1 | de$log2FoldChange > 1,]
  de <- de[c("Transcrito","log2FoldChange","padj")]
  de <- na.omit(de)
  print(dim(de))
  write.table(de, file = paste("../Sponge/", subtype, "_DE_Sponge.tab", sep = ""), row.names = F, quote = F, col.names = T, sep = "\t")
  }



for (subtype in subtipos){
  de_inicial <- ""
  if (subtype != "all") {
     de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  } else {
     de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  }
  colnames(de_inicial)[1] <- "Transcrito"
  de <- merge(de_inicial, BioType_ids, by.x = "Transcrito", by.y = "gene_id")
  de$diffexpressed <- "NO"
  de$important <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 2 & de$padj < 0.05] <- "UP"
  de$diffexpressed[de$gene_type == "lncRNA" & de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
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
    geom_point(size = 0.4) + 
    theme_minimal() +
    geom_text_repel(segment.size = 0.2, force = 8) +
    scale_color_manual(values=c("blue", "black", "red")) +
    #geom_vline(xintercept=c(-2, 2), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    theme(plot.title = element_text(hjust = 0, size=20, face="bold.italic", vjust= 0))  
  ggsave(paste("Volcano_plot_", subtype,"_samples.png", sep = ""), width = 12, height = 6, dpi = 500, units = "in")
  system("scp -P 1313 Volcano_plot_all_samples.png adolfo@200.89.65.156:/media/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/DE")

    if (subtype == "all") {
     mdp <- read.delim("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/all_samples_sample_scores.tsv", sep = " ")
  } else {
     mdp <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype, "_all_controls_sample_scores.tsv", sep = ""), sep = " ")
  }
  
  ggdensity(mdp, x = "allgenes.Score",
    add = "mean", rug = TRUE,
    color = "allgenes.Class", fill = "allgenes.Class",
    palette = c("#00AFBB", "#E7B800"))+
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples", sep = ""))+ 
    theme(legend.position = c(0.75,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), 
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.direction = "horizontal")+
    ylab("Density")+ xlab("Sample score - all") + labs(fill = "Sample type", color ="Sample type")
  ggsave(paste("Density_plot_", subtype,"_samples_distribtion.png", sep = ""), width = 16, height = 9, dpi = 500, units = "in")
  
  
  ggbarplot(mdp, x = "allgenes.Sample", y = "allgenes.Score",
    fill = "allgenes.Class",               # change fill color by cyl
    color = "allgenes.Class",            # Set bar border colors to white
    palette = "jco",            # jco journal color palett. see ?ggpar
    sort.val = "asc",           # Sort the value in dscending order
    sort.by.groups = TRUE,      # Sort inside each group
    x.text.angle = 90)+         # Rotate vertically x axis texts
    theme(axis.text.x=element_blank(), legend.position = c(0.15,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), 
    legend.title = element_text(size = 15), legend.text = element_text(size = 12), legend.direction = "horizontal")+
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples", sep = ""))+
    ylab("All genes score")+ xlab("Samples")+ labs(fill = "Sample type", color = "Sample type")
  ggsave(paste("Bar_plot_", subtype,"_samples_distribtion.png", sep =""), width = 16, height = 9, dpi = 500, units = "in")
  
  
  stat.test <- compare_means(allgenes.Score ~ allgenes.Class, mdp, comparisons = my_comparisons)
  stat.test$p<- formatC(stat.test$p, format = "e", digits = 2)
  ggboxplot(mdp, x = "allgenes.Class", y = "allgenes.Score",
            color = "allgenes.Class", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
            add = "jitter", shape = "allgenes.Class") +
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples", sep = ""))+
    theme(legend.position = "none")+
    stat_pvalue_manual(stat.test, label = "{p.signif}\np ={p} {method} test", y.position = max(mdp$allgenes.Score)+0.5, vjust = 0.8)+
    ylab("All genes score")+ xlab("Sample type")
  ggsave(paste("Scatter_plot_", subtype,"_samples.png", sep = ""), width = 16, height = 9, dpi = 500, units = "in")  

#############################################################################################################################################
#                   miRNAs
#############################################################################################################################################

  de <- ""
  if (subtype != "all") {
     de <- read.delim(paste("Corregido_resultados_expresion_diferencial_miRNAs-gene_level/", subtype,"_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  } else {
     de <- read.delim(paste("Corregido_resultados_expresion_diferencial_miRNAs-gene_level/", subtype,"_samples_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  }
  
  colnames(de)[1] <- "Transcrito"
  de$diffexpressed <- "NO"
  de$important <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 1 & de$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2FoldChange < -1 & de$padj < 0.05] <- "DOWN"
  de <- de[order(de$log2FoldChange, decreasing = T),]
  de[1:3,]$important <- "UP"
  de[(nrow(de)-2):nrow(de),]$important <- "DOWN"  
  de$delabel <- NA
  de$delabel[de$important != "NO"] <- de$Transcrito[de$important != "NO"]
  de$delabel <- gsub("hsa-", "", de$delabel)
  # plot adding up all layers we have seen so far
  p2 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    #ggtitle(paste("Breast cancer diferential expression analysis\nmiRNAs\n", subtype, " samples", sep = ""))+
    geom_point(size = 0.4) + theme_minimal() +
    geom_text_repel(segment.size = 0.2, force = 8) +
    scale_color_manual(values=c("blue", "black", "red")) +
    #geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    theme(plot.title = element_text(hjust = 0, size=20, face="bold.italic", vjust= 0))
  
  ggsave(paste("Volcano_plot_", subtype,"_samples-miRNAs.png", sep =""), width = 12, height = 6, dpi = 500, units = "in")
    
  if (subtype == "all") {
     mdp <- read.delim("Corregido_resultados_expresion_diferencial_miRNAs-gene_level/all_samples_sample_scores.tsv", sep = " ")
  } else {
     mdp <- read.delim(paste("Corregido_resultados_expresion_diferencial_miRNAs-gene_level/", subtype, "_all_controls_sample_scores.tsv", sep = ""), sep = " ")
  }
  
p1 <- ggdensity(mdp, x = "allgenes.Score",
            add = "mean", rug = TRUE,
            color = "allgenes.Class", fill = "allgenes.Class",
            palette = c("#00AFBB", "#E7B800")) +
            ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples - miRNAs", sep = "")) +
            theme(legend.position = c(0.75,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0),          
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            legend.direction = "horizontal")+
            ylab("Density")+ xlab("Sample score - all") + labs(fill = "Sample type", color ="Sample type")
  ggsave(paste("Density_plot_", subtype,"_samples_distribtion-miRNAs.png", sep = ""), width = 16, height = 9, dpi = 500, units = "in")
  
  
p2 <- ggbarplot(mdp, x = "allgenes.Sample", y = "allgenes.Score",
            fill = "allgenes.Class",               # change fill color by cyl
            color = "allgenes.Class",            # Set bar border colors to white
            palette = "jco",            # jco journal color palett. see ?ggpar
            sort.val = "asc",           # Sort the value in dscending order
            sort.by.groups = TRUE,      # Sort inside each group
            x.text.angle = 90) +        # Rotate vertically x axis texts
            theme(axis.text.x=element_blank(), legend.position = c(0.15,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), 
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            legend.direction = "horizontal")+
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples - miRNAs", sep = ""))+
    ylab("All genes score")+ xlab("Samples")+ labs(fill = "Sample type", color = "Sample type")
  ggsave(paste("Bar_plot_", subtype,"_samples_distribtion-miRNAs.png", sep = ""), width = 16, height = 9, dpi = 500, units = "in")
  
  stat.test <- compare_means(allgenes.Score ~ allgenes.Class, mdp, comparisons = my_comparisons)
  stat.test$p<- formatC(stat.test$p, format = "e", digits = 2)
p3 <-  ggboxplot(mdp, x = "allgenes.Class", y = "allgenes.Score",
            color = "allgenes.Class", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
            add = "jitter", shape = "allgenes.Class") +
    ggtitle(paste("Molecular degree of perturbation, ", subtype, " samples - miRNAs", sep = ""))+
    theme(legend.position = "none")+
    stat_pvalue_manual(stat.test, label = "{p.signif}\np ={p} {method} test", y.position = max(mdp$allgenes.Score)+0.5, vjust = 0.8)+
    ylab("All genes score")+ xlab("Sample type")
  ggsave(paste("Scatter_plot_", subtype,"_samples-miRNAs.png", sep =""), width = 16, height = 9, dpi = 500, units = "in")

 p1 <- ggarrange(p1, p3, labels = c("A", "B"), ncol = 2, nrow = 1)
 ggarrange(p1, p2, labels = c("", "C"), ncol = 1, nrow = 2)

 ggsave("filename.png", width = 16, height = 9, dpi = 600, units = "in")
 ggsave("Figura_volcano_RNA-miRNA.png", width = 12, height = 12, dpi = 600, units = "in")
 system("scp -P 1313 Figura_volcano_RNA-miRNA.png adolfo@200.89.65.156:/media/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/DE")
  
  }
  system("scp -P 1313 Volcano_plot_all_samples-miRNAs.png adolfo@200.89.65.156:/media/run-projects/Adolfo/Volcano_plot_all_samples-miRNAs.png")
  system("scp -P 1313 Volcano_plot_all_samples.png adolfo@200.89.65.156:/media/run-projects/Adolfo/Volcano_plot_all_samples.png")
  #system("scp -P 1313 filename.png adolfo@200.89.65.156:/media/run-projects/Adolfo/filename.png")