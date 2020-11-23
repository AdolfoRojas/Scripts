#!/usr/bin/env Rscript
require(TCGAbiolinks)
library(pdftools)
library(DESeq2)
library(ggplot2)
library(BiocParallel)
register(MulticoreParam(20))
system("mkdir -p Resultados_expresion_diferencial_mRNAs")
Result_directory <- "Resultados_expresion_diferencial_mRNAs/"

matrix <- read.delim("TCGA_BRCA_counts.tsv")
sample_data <- read.delim("Samples_Subtype_BRCA.tsv")
sample_data <- sample_data[sample_data$BRCA_Subtype_PAM50 != "Normal",]
sample_data <- sample_data[sample_data$BRCA_Subtype_PAM50 != "NA",]

sample_data$CGHubAnalysisID <- gsub("-", "\\.", sample_data$CGHubAnalysisID)

rownames(matrix) <- sapply(strsplit(as.character(matrix$X),'\\.'), "[", 1)
matrix$X <- NULL
matrix <- round(matrix, 0)

for (AnalysisID in sample_data$CGHubAnalysisID) {
    names(matrix)[grep(substr(AnalysisID, start = 4 , stop = 40 ), names(matrix))] <- AnalysisID
	}
matrix <- matrix[sample_data$CGHubAnalysisID]
sample_data$sample_type<- gsub("TP", "Tumoral", sample_data$sample_type)
sample_data$sample_type<- gsub("NT", "Normal", sample_data$sample_type)
sample_data <- sample_data[c("CGHubAnalysisID", "sample_type", "BRCA_Subtype_PAM50")]
colnames(sample_data) <- c("column", "Condition", "Subtype")
sample_data <- sample_data[is.na(sample_data$Subtype) != T,]
sample_data <- as.data.frame(sample_data)

######################################################################################################################################################
coldata <- sample_data
cts <- matrix[coldata$column] + 1

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
keep <- rowSums(counts(dds)) > length(cts) ## Pre-Filtering
dds <- dds[keep,]
dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
dds$Condition <- droplevels(dds$Condition)
dds <- DESeq(dds, fitType="local", parallel = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Subtype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Subtype)) +
  ggtitle("All samples for mRNAs\n diferential expression analysis") +
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values=c("#56B4E9", "red"))
ggsave(paste(Result_directory, "all_samples_mRNAs_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in")
res <- results(dds, parallel = TRUE)
pdf("all_samples_mRNAs_MA_plotx.pdf", height = 6.5, width = 10)
plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in all samples mRNAs analysis\n", sep = " "), ylim=c(-10,10))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()
bitmap <- pdf_render_page("all_samples_mRNAs_MA_plotx.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, paste(Result_directory, "all_samples_mRNAs_MA_plot.png", sep = ""))
unlink("all_samples_mRNAs_MA_plotx.pdf")
write.csv( res, file=paste(Result_directory, "all_samples_mRNAs_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))


TCGA_files <- c("LumA", "LumB", "Her2", "Basal")
for (File in TCGA_files){
    coldata <- sample_data[sample_data$Subtype == File | sample_data$Condition == "Normal",]
    cts <- matrix[coldata$column] + 1
    
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
    keep <- rowSums(counts(dds)) > length(cts) ## Pre-Filtering
    dds <- dds[keep,]
    dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
    dds$Condition <- droplevels(dds$Condition)
    dds <- DESeq(dds, fitType="local", parallel = TRUE)
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup=c("Condition", "Subtype"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Subtype)) +
      ggtitle(paste(File, "samples for mRNAs\n diferential expression analysis", sep = " ")) +
      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
      geom_point(size=1) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      scale_color_manual(values=c("#56B4E9", "red"))
    ggsave(paste(Result_directory, File, "_samples_all_controls_mRNAs_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in") #####################
    res <- results(dds, parallel = TRUE)
    pdf("all_samples_mRNAs_MA_plotx.pdf", height = 6.5, width = 10)
    plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "samples mRNAs analysis\n", sep = " "), ylim=c(-10,10))
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    bitmap <- pdf_render_page("all_samples_mRNAs_MA_plotx.pdf", page = 1, dpi = 300)
    png::writePNG(bitmap, paste(Result_directory, File, "_samples_all_controls_mRNAs_MA_plot.png", sep = "")) #####################
    unlink("all_samples_mRNAs_MA_plotx.pdf")
    write.csv( res, file=paste(Result_directory, File, "_samples_all_controls_mRNAs_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))#######################

    ####################################################################################################################################################################

    coldata <- sample_data[sample_data$Subtype == File,]
    cts <- matrix[coldata$column] + 1
    
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
    keep <- rowSums(counts(dds)) > length(cts) ## Pre-Filtering
    dds <- dds[keep,]
    dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
    dds$Condition <- droplevels(dds$Condition)
    dds <- DESeq(dds, fitType="local", parallel = TRUE)
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup=c("Condition", "Subtype"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Subtype)) +
      ggtitle(paste(File, "samples for mRNAs\n diferential expression analysis", sep = " ")) +
      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
      geom_point(size=1) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      scale_color_manual(values=c("#56B4E9", "red"))
    ggsave(paste(Result_directory, File, "_samples_subtype_controls_mRNAs_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in") #####################
    res <- results(dds, parallel = TRUE)
    pdf("all_samples_mRNAs_MA_plotx.pdf", height = 6.5, width = 10)
    plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "samples mRNAs analysis\n", sep = " "), ylim=c(-10,10))
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    bitmap <- pdf_render_page("all_samples_mRNAs_MA_plotx.pdf", page = 1, dpi = 300)
    png::writePNG(bitmap, paste(Result_directory, File, "_samples_subtype_controls_mRNAs_MA_plot.png", sep = "")) #####################
    unlink("all_samples_mRNAs_MA_plotx.pdf")
    write.csv( res, file=paste(Result_directory, File, "_samples_subtype_controls_mRNAs_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))#######################    
}
system(paste("tar -czvf resultados_mRNAs.tar.gz", Result_directory, sep = " "))