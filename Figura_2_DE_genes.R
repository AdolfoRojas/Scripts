#!/usr/bin/env Rscript
subtipos <- c("Her2", "LumB", "LumA", "Basal")

df <- data.frame()
for (subtype in subtipos){
  de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  colnames(de_inicial)[1] <- "Transcrito"
  de_inicial$subtype <- subtype
  de_inicial <- de_inicial[de_inicial$padj < 0.05,]
  de_inicial <- de_inicial[de_inicial$log2FoldChange < -1 | de_inicial$log2FoldChange > 1,]
  df <- rbind(df,de_inicial)
  df <- na.omit(df)
}
gencode <- read.delim("gencodeV35_fixed.txt", sep = "\t")
gencode <- unique(gencode[c("gene_id", "gene_type")])
names(gencode)[2] <- "biotype"
df <- merge(df, gencode, by.x = "Transcrito", by.y = "gene_id")
df2 <- data.frame()
for (subtype in subtipos){
  de_inicial <- read.delim(paste("Corregido_resultados_expresion_diferencial_miRNAs-gene_level/", subtype,"_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
  colnames(de_inicial)[1] <- "Transcrito"
  de_inicial$subtype <- subtype
  de_inicial$biotype <- "miRNA"
  de_inicial <- de_inicial[de_inicial$padj < 0.05,]
  de_inicial <- de_inicial[de_inicial$log2FoldChange < -1 | de_inicial$log2FoldChange > 1,]
  df2 <- rbind(df2,de_inicial)
  df2 <- na.omit(df2)
}


df[df$subtype == "LumA",]$subtype <- "Luminal A"
df[df$subtype == "LumB",]$subtype <- "Luminal B"
df2[df2$subtype == "LumA",]$subtype <- "Luminal A"
df2[df2$subtype == "LumB",]$subtype <- "Luminal B"


df$regulation <- "Up"
df[df$log2FoldChange < 0,]$regulation <- "Down"

df1 <- df[df$biotype == "protein_coding",]
df1 <- df1[df1$log2FoldChange < -2 | df1$log2FoldChange > 2,]
df3 <- df[df$biotype == "lncRNA",]
df4 <- df[df$biotype != "protein_coding" & df$biotype != "lncRNA",]

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(df1[df1$subtype == "Basal" & df1$log2FoldChange > 0,]$Transcrito,
  df1[df1$subtype == "Her2" & df1$log2FoldChange > 0,]$Transcrito,
  df1[df1$subtype == "Luminal A" & df1$log2FoldChange > 0,]$Transcrito, 
  df1[df1$subtype == "Luminal B" & df1$log2FoldChange > 0,]$Transcrito),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_Protein_coding_up.png',
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
  x = list(df1[df1$subtype == "Basal" & df1$log2FoldChange < 0,]$Transcrito,
  df1[df1$subtype == "Her2" & df1$log2FoldChange < 0,]$Transcrito,
  df1[df1$subtype == "Luminal A" & df1$log2FoldChange < 0,]$Transcrito, 
  df1[df1$subtype == "Luminal B" & df1$log2FoldChange < 0,]$Transcrito),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_Protein_coding_down.png',
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
  x = list(df3[df3$subtype == "Basal" & df3$log2FoldChange < 0,]$Transcrito,
  df3[df3$subtype == "Her2" & df3$log2FoldChange < 0,]$Transcrito,
  df3[df3$subtype == "Luminal A" & df3$log2FoldChange < 0,]$Transcrito, 
  df3[df3$subtype == "Luminal B" & df3$log2FoldChange < 0,]$Transcrito),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_lncRNA_down.png',
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
  x = list(df3[df3$subtype == "Basal" & df3$log2FoldChange > 0,]$Transcrito,
  df3[df3$subtype == "Her2" & df3$log2FoldChange > 0,]$Transcrito,
  df3[df3$subtype == "Luminal A" & df3$log2FoldChange > 0,]$Transcrito, 
  df3[df3$subtype == "Luminal B" & df3$log2FoldChange > 0,]$Transcrito),
  category.names = c("Basal","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_lncRNA_up.png',
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
  x = list(df2[df2$subtype == "Basal" & df2$log2FoldChange > 0,]$Transcrito,
  df2[df2$subtype == "Her2" & df2$log2FoldChange > 0,]$Transcrito,
  df2[df2$subtype == "Luminal A" & df2$log2FoldChange > 0,]$Transcrito, 
  df2[df2$subtype == "Luminal B" & df2$log2FoldChange > 0,]$Transcrito),
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
  x = list(df2[df2$subtype == "Basal" & df2$log2FoldChange < 0,]$Transcrito,
  df2[df2$subtype == "Her2" & df2$log2FoldChange < 0,]$Transcrito,
  df2[df2$subtype == "Luminal A" & df2$log2FoldChange < 0,]$Transcrito, 
  df2[df2$subtype == "Luminal B" & df2$log2FoldChange < 0,]$Transcrito),
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
   
system("scp Venn_DE_*.png adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/")



############################################################################################
library(plyr)
library(ggplot2)
counts <- ddply(df1, .(df1$subtype, df1$regulation), nrow) # protein_coding
names(counts) <- c("Subtype", "Regulation", "Freq")

counts2 <- ddply(df3, .(df3$subtype, df3$regulation), nrow) # lncRNA
names(counts2) <- c("Subtype", "Regulation", "Freq")

counts4 <- ddply(df4, .(df4$subtype, df4$regulation), nrow) # Others
names(counts4) <- c("Subtype", "Regulation", "Freq")


df2$regulation <- "Up"
df2[df2$log2FoldChange < 0,]$regulation <- "Down"
counts3 <- ddply(df2, .(df2$subtype, df2$regulation), nrow) # miRNA
names(counts3) <- c("Subtype", "Regulation", "Freq")

protein_coding <- ggplot(counts, aes(x = Subtype, y = Regulation)) + 
  ggtitle("Protein coding")+
  geom_point(aes(color = Regulation, size = Freq), alpha = 0.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_size(range = c(20, 35)) + # Adjust the range of points size
  geom_text(label= counts$Freq, size=6, color = "white")+ 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), plot.title = element_text(face = "bold", hjust = 0.5), axis.text.y = element_blank())

lncRNA <- ggplot(counts2, aes(x = Subtype, y = Regulation)) + 
  ggtitle("lncRNAs")+
  geom_point(aes(color = Regulation, size = Freq), alpha = 0.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_size(range = c(20, 35)) + # Adjust the range of points size
  geom_text(label= counts2$Freq, size=6, color = "white")+ 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), plot.title = element_text(face = "bold", hjust = 0.5), axis.text.y = element_blank())
miRNA <- ggplot(counts3, aes(x = Subtype, y = Regulation)) + 
  ggtitle("miRNAs")+
  geom_point(aes(color = Regulation, size = Freq), alpha = 0.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_size(range = c(20, 35)) + # Adjust the range of points size
  geom_text(label= counts3$Freq, size=6, color = "white")+ 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), plot.title = element_text(face = "bold", hjust = 0.5), axis.text.y = element_blank())
Others <- ggplot(counts4, aes(x = Subtype, y = Regulation)) + 
  ggtitle("Others")+
  geom_point(aes(color = Regulation, size = Freq), alpha = 0.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_size(range = c(20, 35)) + # Adjust the range of points size
  geom_text(label= counts4$Freq, size=6, color = "white")+ 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), plot.title = element_text(face = "bold", hjust = 0.5), axis.text.y = element_blank())
require(gridExtra)
plot_all <- grid.arrange(protein_coding, lncRNA, miRNA, Others, ncol=2, nrow = 2)
ggsave(plot = plot_all,filename = "Figura_DE_genes.png", height = 9, width = 16, dpi = 600)
system("scp Figura_DE_genes.png adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Figura_DE_genes.png")
