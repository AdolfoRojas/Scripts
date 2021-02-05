setwd("../Tesis/Expresion Diferencial/")

Basal <- read.delim("Resultados_expresion_diferencial_mRNAs/Basal_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
LumA <- read.delim("Resultados_expresion_diferencial_mRNAs/LumA_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
LumB <- read.delim("Resultados_expresion_diferencial_mRNAs/LumB_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
Her2 <- read.delim("Resultados_expresion_diferencial_mRNAs/Her2_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")

Basal_all <- Basal[Basal$padj < 0.05 & (Basal$log2FoldChange < -2| Basal$log2FoldChange > 2),]
LumA_all <- LumA[LumA$padj < 0.05 & (LumA$log2FoldChange < -2| LumA$log2FoldChange > 2),]
LumB_all <- LumB[LumB$padj < 0.05 & (LumB$log2FoldChange < -2| LumB$log2FoldChange > 2),]
Her2_all <- Her2[Her2$padj < 0.05 & (Her2$log2FoldChange < -2| Her2$log2FoldChange > 2),]

Basal_up <- Basal[Basal$padj < 0.05 & Basal$log2FoldChange > 2,]
LumA_up <- LumA[LumA$padj < 0.05 & LumA$log2FoldChange > 2,]
LumB_up <- LumB[LumB$padj < 0.05 & LumB$log2FoldChange > 2,]
Her2_up <- Her2[Her2$padj < 0.05 & Her2$log2FoldChange > 2,]

Basal_down <- Basal[Basal$padj < 0.05 & Basal$log2FoldChange < -2,]
LumA_down <- LumA[LumA$padj < 0.05 & LumA$log2FoldChange < -2,]
LumB_down <- LumB[LumB$padj < 0.05 & LumB$log2FoldChange < -2,]
Her2_down <- Her2[Her2$padj < 0.05 & Her2$log2FoldChange < -2,]


library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(Basal_all$X, Her2_all$X, LumA_all$X, LumB_all$X),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_RNA_all.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")

venn.diagram(
  x = list(Basal_up$X, Her2_up$X, LumA_up$X, LumB_up$X),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_RNA_up.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")

venn.diagram(
  x = list(Basal_down$X, Her2_down$X, LumA_down$X, LumB_down$X),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_RNA_down.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")

Basal <- read.delim("Resultados_expresion_diferencial_miRNAs/Basal_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
LumA <- read.delim("Resultados_expresion_diferencial_miRNAs/LumA_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
LumB <- read.delim("Resultados_expresion_diferencial_miRNAs/LumB_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
Her2 <- read.delim("Resultados_expresion_diferencial_miRNAs/Her2_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")

Basal_all <- Basal[Basal$padj < 0.05 & (Basal$log2FoldChange < -1| Basal$log2FoldChange > 1),]
LumA_all <- LumA[LumA$padj < 0.05 & (LumA$log2FoldChange < -1| LumA$log2FoldChange > 1),]
LumB_all <- LumB[LumB$padj < 0.05 & (LumB$log2FoldChange < -1| LumB$log2FoldChange > 1),]
Her2_all <- Her2[Her2$padj < 0.05 & (Her2$log2FoldChange < -1| Her2$log2FoldChange > 1),]

Basal_up <- Basal[Basal$padj < 0.05 & Basal$log2FoldChange > 1,]
LumA_up <- LumA[LumA$padj < 0.05 & LumA$log2FoldChange > 1,]
LumB_up <- LumB[LumB$padj < 0.05 & LumB$log2FoldChange > 1,]
Her2_up <- Her2[Her2$padj < 0.05 & Her2$log2FoldChange > 1,]

Basal_down <- Basal[Basal$padj < 0.05 & Basal$log2FoldChange < -1,]
LumA_down <- LumA[LumA$padj < 0.05 & LumA$log2FoldChange < -1,]
LumB_down <- LumB[LumB$padj < 0.05 & LumB$log2FoldChange < -1,]
Her2_down <- Her2[Her2$padj < 0.05 & Her2$log2FoldChange < -1,]

venn.diagram(
  x = list(Basal_all$X, Her2_all$X, LumA_all$X, LumB_all$X),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_miRNA_all.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")

venn.diagram(
  x = list(Basal_up$X, Her2_up$X, LumA_up$X, LumB_up$X),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_miRNA_up.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")

venn.diagram(
  x = list(Basal_down$X, Her2_down$X, LumA_down$X, LumB_down$X),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_miRNA_down.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")
