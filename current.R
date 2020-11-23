#!/usr/bin/env Rscript
library(pdftools)
library(DESeq2)
library(ggplot2)
library(BiocParallel)
register(MulticoreParam(5))
GTEx_controls <- read.csv("Female_Matrix_replaced_breast_tissue_counts_+1_less_PAR_transcripts.tab",sep="\t")
TCGA_files <- c("Her2", "Basal", "LumA", "LumB")
for (File in TCGA_files){
  TCGA_Data <- read.csv(paste(File,"_counts_all_+1.tab", sep=""),sep="\t")
  de <- merge(GTEx_controls, TCGA_Data, by=1)
  row.names(de) <- de[,1]
  print("row names asignados")
  de[,1] <- NULL
  print("columna eliminada")
  de <- round(de, 0)
  print("Redondeo")
  cts <- as.matrix(de)
  print("As matrix")
  coldata <- data.frame(column=character(), Condition=character(), value=character(), stringsAsFactors=FALSE)
  for (column in colnames(cts)){
    if (substr(column, 1, 2) == "NT"){
      Condition_loop= "Normal"
      Type_loop = "TCGA"
      } else if (substr(column, 1, 2) == "TP"){
        Condition_loop= "Tumoral"
        Type_loop = "TCGA"
      } else if (substr(column, 1, 2) == "GT"){
        Condition_loop= "Normal"
        Type_loop = "GTEx"
      }
    coldata <- rbind(coldata, data.frame(column, Condition_loop, Type_loop, fix.empty.names = F))}
  colnames(coldata) <- c("column", "Condition", "Type")
  row.names(coldata) <- coldata[,1]
  coldata[,1] <- NULL
  ######################################################################################################
  Part_1 = "all_controls"
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

  keep <- rowSums(counts(dds)) >= 10 ## Pre-Filtering
  dds <- dds[keep,]
  dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
  dds$Condition <- droplevels(dds$Condition)
  dds <- DESeq(dds, parallel = TRUE)

  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("Condition", "Type"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Type)) +
    ggtitle(paste("Samples for", File, "analysis\n", Part_1, sep =" ")) + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12, face="bold.italic"), axis.title.y = element_text(size=12, face="bold.italic")) + 
    geom_point(size=1) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + 
    scale_color_manual(values=c("#56B4E9", "red"))
  #ggsave("mtcars.pdf", width = 20, height = 20, units = "cm")
  print("grafico Generado")
  pdf(paste(File, "_PCA_plotx.pdf", sep=""), height = 8.5, width = 8.5)
  print(p)
  dev.off()
  bitmap <- pdf_render_page(paste(File, "_PCA_plotx.pdf", sep=""), page = 1, dpi = 300)
  png::writePNG(bitmap, paste(File, "_PCA_plot_",Part_1,".png",sep = ""))

  res <- results(dds, parallel = TRUE)
  pdf(paste(File, "_MA_plotx.pdf", sep=""), height = 6.5, width = 10)
  plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "analysis\n", Part_1, sep =" "), ylim=c(-10,10))
  dev.off()
  bitmap <- pdf_render_page(paste(File, "_MA_plotx.pdf", sep=""), page = 1, dpi = 300)
  png::writePNG(bitmap, paste(File, "_MA_plot_", Part_1, ".png",sep = ""))
  write.csv( res, file=paste(File, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE_", Part_1, ".tab", sep= ""))
  ####################################################################################################################
  Part_2 = "without_TCGA_controls"
  coldata2 <- coldata[coldata$Condition == "Tumoral" | coldata$Type == "GTEx",]
  cts2 <- cts[, row.names(coldata2)]
  dds <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata2, design = ~ Condition)

  keep <- rowSums(counts(dds)) >= 10 ## Pre-Filtering
  dds <- dds[keep,]
  dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
  dds$Condition <- droplevels(dds$Condition)
  dds <- DESeq(dds, parallel = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("Condition", "Type"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Type)) + 
    ggtitle(paste("Samples for", File, "analysis\n", Part_2, sep =" ")) + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12, face="bold.italic"), axis.title.y = element_text(size=12, face="bold.italic")) + 
    geom_point(size=1) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + 
    scale_color_manual(values=c("#56B4E9", "red"))
  pdf(paste(File, "_PCA_plotx.pdf", sep=""), height = 8.5, width = 8.5)
  print(p)
  dev.off()
  bitmap <- pdf_render_page(paste(File, "_PCA_plotx.pdf", sep=""), page = 1, dpi = 300)
  png::writePNG(bitmap, paste(File, "_PCA_plot_",Part_2,".png",sep = ""))
  res <- results(dds, parallel = TRUE)
  pdf(paste(File, "_MA_plotx.pdf", sep=""), height = 6.5, width = 10)
  plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "analysis\n", Part_2, sep =" "), ylim=c(-10,10))
  dev.off()
  bitmap <- pdf_render_page(paste(File, "_MA_plotx.pdf", sep=""), page = 1, dpi = 300)
  png::writePNG(bitmap, paste(File, "_MA_plot_", Part_2, ".png",sep = ""))
  write.csv( res, file=paste(File, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE_", Part_2, ".tab", sep= ""))
  ####################################################################################################################
  Part_3 = "without_GTEx_controls"  
  coldata3 <- coldata[coldata$Type != "GTEx",]
  cts3 <- cts[, row.names(coldata3)]

  dds <- DESeqDataSetFromMatrix(countData = cts3, colData = coldata3, design = ~ Condition)

  keep <- rowSums(counts(dds)) >= 10 ## Pre-Filtering
  dds <- dds[keep,]
  dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
  dds$Condition <- droplevels(dds$Condition)
  dds <- DESeq(dds, parallel = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("Condition", "Type"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Type)) + 
    ggtitle(paste("Samples for", File, "analysis\n", Part_3, sep =" ")) + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12, face="bold.italic"), axis.title.y = element_text(size=12, face="bold.italic")) + 
    geom_point(size=1) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + 
    scale_color_manual(values=c("#56B4E9", "red"))
  pdf(paste(File, "_PCA_plotx.pdf", sep=""), height = 8.5, width = 8.5)
  print(p)
  dev.off()
  bitmap <- pdf_render_page(paste(File, "_PCA_plotx.pdf", sep=""), page = 1, dpi = 300)
  png::writePNG(bitmap, paste(File, "_PCA_plot_",Part_3,".png",sep = ""))
  res <- results(dds, parallel = TRUE)
  pdf(paste(File, "_MA_plotx.pdf", sep=""), height = 6.5, width = 10)
  plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "analysis\n", Part_3, sep =" "), ylim=c(-10,10))
  dev.off()
  bitmap <- pdf_render_page(paste(File, "_MA_plotx.pdf", sep=""), page = 1, dpi = 300)
  png::writePNG(bitmap, paste(File, "_MA_plot_", Part_3, ".png",sep = ""))
  write.csv( res, file=paste(File, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE_", Part_3,".tab", sep= ""))
  #unlink(paste(File, "_PCA_plotx.pdf", sep=""))
  unlink(paste(File, "_MA_plotx.pdf", sep=""))

  }





















  #!/usr/bin/env python3
import pandas as pd
print("Iniciando Carga")        
raw_matrix = pd.read_csv("expressed_transcripts.txt", sep="\t",header=0)  # Requerido
annotated_matrix = pd.read_csv("../ID_matrix_New_lncRNAs_anotated.anotate", sep="|",header=0).drop(columns="8") # Requerido
raw_matrix["transcript_id"] = raw_matrix["transcript_id"].str.split(".", n=1, expand = True)[0]
    ### ID_matrix_New_lncRNAs_anotated.anotate expressed_transcripts.txt

merge_data = raw_matrix.merge(annotated_matrix, left_on= "transcript_id", right_on= "Transcripts").drop(columns=["Transcripts","1","2","3","6"])
without_anotarion = merge_data.merge(raw_matrix, on= "transcript_id", how = "outer")
without_anotarion = without_anotarion[without_anotarion['7'].isnull()].drop(columns=["4","5","7"])
print(merge_data["7"].value_counts())
    
merge_data.to_csv("Anotacion_incompleta.txt", sep="\t",header=True, index=False)
without_anotarion.to_csv("Transcritos_no_anotados.txt", sep="\t",header=True, index=False)

