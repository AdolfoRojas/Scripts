#!/usr/bin/env Rscript
require(TCGAbiolinks)
library(pdftools)
library(DESeq2)
library(ggplot2)
library(BiocParallel)
register(MulticoreParam(20))
CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("GDC/",gsub("-","_",CancerProject))
system("mkdir -p Resultados_expresion_diferencial_miRNAs")
Result_directory <- "Resultados_expresion_diferencial_miRNAs/"
FileNameData <- paste0(DataDirectory, "_","miRNA_gene_quantification",".rda")
query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE)
samplesDown.miR <- getResults(query.miR,cols=c("cases"))
dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
GDCdownload(query = queryDown.miR, directory = DataDirectory)
dataAssy.miR <- GDCprepare(query = queryDown.miR, save = TRUE, save.filename = FileNameData, summarizedExperiment = TRUE, directory =DataDirectory)
rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID
read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))] # using read_count's data
matrix <- dataAssy.miR[, read_countData]
colnames(matrix) <- gsub("read_count_","", colnames(matrix))
names <- rownames(dataAssy.miR)
matrix<- cbind(names,matrix)
write.table(matrix, sep = "\t", file = "miRNAs_counts.tab", row.names = F, quote = F, col.names = T)
write.table(dataSmTP.miR, sep = "\t", file = "miRNAs_TP_IDs.tab", row.names = F, quote = F, col.names = T)
write.table(dataSmNT.miR, sep = "\t", file = "miRNAs_NT_IDs.tab", row.names = F, quote = F, col.names = T)

clinical_data <- GDCquery_clinic(CancerProject, type = "clinical", save.csv = F)
clinical_data <- clinical_data[c("submitter_id", "primary_diagnosis", "ethnicity", "race", "gender","age_at_diagnosis")] # Caracteristicas Clinicaas de Interes
Subtype_BCa <- TCGAquery_subtype(tumor = "brca")
Subtype_BCa <- Subtype_BCa[c(1,9,12)]
names(clinical_data)[names(clinical_data) == "submitter_id"] <- "patient"
df <- as.data.frame(merge(clinical_data, Subtype_BCa, by = "patient"))
df <- df[is.na(df$gender) != T,]
df <- df[df$gender == "female",]
df <- df[df$BRCA_Subtype_PAM50 != "Normal",]
df <- df[is.na(df$BRCA_Subtype_PAM50) != T,]
df <- df[df$BRCA_Subtype_PAM50 != "NA",]
dataSmTP.miR2 <- data.frame(sample_id= dataSmTP.miR, sample_type = rep("Tumoral"))
dataSmNT.miR2 <- data.frame(sample_id= dataSmNT.miR, sample_type = rep("Normal"))
NT_TP_samples <- as.data.frame(rbind(dataSmTP.miR2, dataSmNT.miR2))
NT_TP_samples$patient <- substr(NT_TP_samples$sample_id, 1,12)
sample_data <- as.data.frame(merge(NT_TP_samples, df, by = "patient"))
write.table(sample_data, sep = "\t", file = "miRNAs_sample_data.tab", row.names = F, quote = F, col.names = T)

sample_data <- sample_data[c("sample_id", "sample_type", "BRCA_Subtype_PAM50")]
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
  ggtitle("All samples for miRNAs\n diferential expression analysis") +
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values=c("#56B4E9", "red"))
ggsave(paste(Result_directory, "all_samples_miRNAs_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in")
res <- results(dds, parallel = TRUE)
pdf("all_samples_miRNAs_MA_plotx.pdf", height = 6.5, width = 10)
plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in all samples miRNAs analysis\n", sep = " "), ylim=c(-10,10))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()
bitmap <- pdf_render_page("all_samples_miRNAs_MA_plotx.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, paste(Result_directory, "all_samples_miRNAs_MA_plot.png", sep = ""))
unlink("all_samples_miRNAs_MA_plotx.pdf")
write.csv( res, file=paste(Result_directory, "all_samples_miRNAs_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))


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
      ggtitle(paste(File, "samples for miRNAs\n diferential expression analysis", sep = " ")) +
      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
      geom_point(size=1) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      scale_color_manual(values=c("#56B4E9", "red"))
    ggsave(paste(Result_directory, File, "_samples_all_controls_miRNAs_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in") #####################
    res <- results(dds, parallel = TRUE)
    pdf("all_samples_miRNAs_MA_plotx.pdf", height = 6.5, width = 10)
    plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "samples miRNAs analysis\n", sep = " "), ylim=c(-10,10))
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    bitmap <- pdf_render_page("all_samples_miRNAs_MA_plotx.pdf", page = 1, dpi = 300)
    png::writePNG(bitmap, paste(Result_directory, File, "_samples_all_controls_miRNAs_MA_plot.png", sep = "")) #####################
    unlink("all_samples_miRNAs_MA_plotx.pdf")
    write.csv( res, file=paste(Result_directory, File, "_samples_all_controls_miRNAs_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))#######################

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
      ggtitle(paste(File, "samples for miRNAs\n diferential expression analysis", sep = " ")) +
      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
      geom_point(size=1) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      scale_color_manual(values=c("#56B4E9", "red"))
    ggsave(paste(Result_directory, File, "_samples_subtype_controls_miRNAs_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in") #####################
    res <- results(dds, parallel = TRUE)
    pdf("all_samples_miRNAs_MA_plotx.pdf", height = 6.5, width = 10)
    plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "samples miRNAs analysis\n", sep = " "), ylim=c(-10,10))
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    bitmap <- pdf_render_page("all_samples_miRNAs_MA_plotx.pdf", page = 1, dpi = 300)
    png::writePNG(bitmap, paste(Result_directory, File, "_samples_subtype_controls_miRNAs_MA_plot.png", sep = "")) #####################
    unlink("all_samples_miRNAs_MA_plotx.pdf")
    write.csv( res, file=paste(Result_directory, File, "_samples_subtype_controls_miRNAs_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))#######################    
}
system(paste("tar -czvf resultados_miRNAs.tar.gz", Result_directory, sep = " "))