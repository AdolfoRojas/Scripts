library(ggplot2)
library(stringr)
library(tidyr)
library(pdftools)
require(gridExtra)
setwd("../Tesis/co-expresion")
file <- read.delim("normal_vs_tumoral/Tables/enrichment_nes.tsv")
colnames(file)[1] <- "Modules"
df <- pivot_longer(file, cols=2:3, names_to = "Class", values_to = "NES")
df <- df[str_order(df$Modules, numeric = TRUE),]
df$Modules <- factor(df$Modules,levels=rev(unique(df$Modules)))



file <- read.delim("subtipos_vs_normal/Tables/enrichment_nes.tsv")
colnames(file)[1] <- "Modules"
df2 <- pivot_longer(file, cols=2:6, names_to = "Class", values_to = "NES")
df2 <- df2[str_order(df2$Modules, numeric = TRUE),]
df2$Modules <- factor(df2$Modules,levels=rev(unique(df2$Modules)))



p <- ggplot(df, aes(Class, Modules)) +
  ggtitle("A")+
  geom_point(aes(size = abs(NES), colour=NES)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() + theme(plot.title = element_text(face = "bold.italic", hjust= -0.1, size = 24)) + 
  labs(size = "NES")

q <- ggplot(df2, aes(Class, Modules)) +
  ggtitle("B")+
  geom_point(aes(size = abs(NES), colour=NES)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal()+ theme(plot.title = element_text(face = "bold.italic", hjust= -0.1, size = 24)) + 
  labs(size = "NES") + theme(legend.position = "none")


pdf("GSEA_co-expresion.pdf", height = 9, width =16)
grid.arrange(p, q, ncol=2)
dev.off()
bitmap <- pdf_render_page("GSEA_co-expresion.pdf", page = 1, dpi = 600)
png::writePNG(bitmap, "GSEA_co-expresion.png")
unlink("GSEA_co-expresion.pdf")
