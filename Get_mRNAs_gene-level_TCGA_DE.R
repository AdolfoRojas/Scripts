#!/usr/bin/env Rscript
library(pdftools)
library(DESeq2)
library(ggplot2)
library(mdp)
library(BiocParallel)
library(telegram.bot)
library(dplyr)
library(SummarizedExperiment)
bot = Bot(token = bot_token("AARH_95_bot"))
chat_id <- bot_token("chat_id")
register(MulticoreParam(50))
Taruca_adolfo_tesis <- "adolfo@200.89.65.156:/media/run-projects/Adolfo/Resultados_Tesis/"
system("mkdir -p Resultados_expresion_diferencial_mRNAs-gene_level")
system("mkdir -p Resultados_expresion_diferencial_miRNAs-gene_level")
system("mkdir -p Corregido_resultados_expresion_diferencial_miRNAs-gene_level")
system("mkdir -p Corregido_resultados_expresion_diferencial_mRNAs-gene_level")
#my.mode <- 2
for (my.mode in 1:2){
print(my.mode)
#######################################################################################################################################
# Contiene el script tesis/3_co-expresion_analysis/1_Get_mRNAs_gene-level_FPKM_TCGA_matrix.R
#######################################################################################################################################
if (my.mode == 6) {
        require(TCGAbiolinks)        
        CancerProject <- "TCGA-BRCA"
        DataDirectory <- paste0("GDC/",gsub("-","_",CancerProject))
        FileNameData <- paste0(DataDirectory, "_","mRNA_gene_quantification",".rda")
        query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts", legacy = FALSE)
        samplesDown.miR <- getResults(query.miR,cols=c("cases"))
        dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
        dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
        queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - Counts", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
        GDCdownload(query = queryDown.miR, directory = DataDirectory)
        dataAssy.miR <- GDCprepare(query = queryDown.miR, save = TRUE, save.filename = FileNameData, summarizedExperiment = TRUE, directory =DataDirectory)
        matrix <- assay(dataAssy.miR,"HTSeq - Counts")
        write.table(matrix, sep = "\t", file = "1_mRNAs_gene-level_Counts_TCGA_matrix.tab", row.names = T, quote = F, col.names = T)
        #######################################################################################################################################
        #              miRNAs  contiene el script tesis/1_expression_data/TCGA_data/5_Get_miRNAs_count_matrix.R
        #######################################################################################################################################
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
        write.table(matrix, sep = "\t", file = "miRNAs_counts.tab", row.names = F, quote = F, col.names = T)} 
#######################################################################################################################################
expr0 <- read.delim("1_mRNAs_gene-level_Counts_TCGA_matrix.tab")
#######################################################################################################################################
sample_annot <- ""
if (my.mode == 1){
        sample_annot <- read.delim("Samples_Subtype_BRCA.tab", sep = "\t")
        sample_annot <- sample_annot[sample_annot$BRCA_Subtype_PAM50 != "Normal",]
        sample_annot <- sample_annot[sample_annot$BRCA_Subtype_PAM50 != "NA",]
        sample_annot$AliquotBarcode <- substring(sample_annot$AliquotBarcode, 1, 15)
        sample_annot <- sample_annot[duplicated(sample_annot$patient) | duplicated(sample_annot$patient,fromLast=T),]
        colnames(sample_annot)[colnames(sample_annot)=="AliquotBarcode"] <- "SampleName"}
if (my.mode == 2){
        sample_annot <- read.delim("2_sample_annot.tab", sep = "\t") #añadir muestras outliers        
        #sample_annot <- sample_annot[duplicated(sample_annot$patient) | duplicated(sample_annot$patient,fromLast=T),]
        cluster_info <- read.delim("cluster_output_miRNAs.tab", sep = "\t", header=FALSE)
        #cluster_info2 <- read.delim("cluster_output_mRNAs.tab", sep = "\t", header=FALSE)
        sample_annot <- sample_annot[!(sample_annot$SampleName %in% cluster_info$V1),]
        #sample_annot <- sample_annot[!(sample_annot$SampleName %in% cluster_info2$V1),]
        sample_annot <- sample_annot[duplicated(sample_annot$patient) | duplicated(sample_annot$patient,fromLast=T),]
        }
colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
colnames(expr0) <- substring(colnames(expr0), 1, 15)
#######################################################################################################################################
expr1 <- read.delim("miRNAs_counts.tab", sep = "\t") 
colnames(expr1) <- gsub("\\.", "-", colnames(expr1))
colnames(expr1) <- substring(colnames(expr1), 1, 15)
expr0 <- expr0[,colnames(expr0)[colnames(expr0) %in% sample_annot$SampleName]]
expr1 <- expr1[,colnames(expr1)[colnames(expr1) %in% sample_annot$SampleName]]
all(colnames(expr0) %in% sample_annot$SampleName)
all(colnames(expr1) %in% sample_annot$SampleName)
sample_annot <- sample_annot[sample_annot$SampleName %in% colnames(expr0),]
sample_annot <- sample_annot[sample_annot$SampleName %in% colnames(expr1),]
sample_annot <- sample_annot[!duplicated(sample_annot$SampleName),]
all(sample_annot$SampleName %in% colnames(expr0))
all(sample_annot$SampleName %in% colnames(expr1))
expr0 <- expr0[sample_annot$SampleName]
expr1 <- expr1[sample_annot$SampleName]
#######################################################################################################################################
library(edgeR)
expr0 <- DGEList(counts = expr0, samples = sample_annot, genes = rownames(expr0), remove.zeros = F, group = sample_annot$sample_type) ### remove zero originalmente en T
expr0 <- calcNormFactors(expr0, method = "TMM")
expr0 <- expr0$counts
expr0 <- cpm(expr0, log= F)
identical(colnames(expr0), sample_annot$SampleName)

expr1 <- DGEList(counts = expr1, samples = sample_annot, genes = rownames(expr1), remove.zeros = F, group = sample_annot$sample_type) ### remove zero originalmente en T
expr1 <- calcNormFactors(expr1, method = "TMM")
expr1 <- expr1$counts
expr1 <- cpm(expr1, log= F)
identical(colnames(expr1), sample_annot$SampleName)

if (identical(colnames(expr0), sample_annot$SampleName) == T & identical(colnames(expr1), sample_annot$SampleName) == T){
        print("guardando tablas y procediendo con anotacion en gene symbol de interacciones proteina-proteina")
        if (my.mode == 1){
                write.table(expr1, file = "2_miRNAs_gene-level_TMM_TCGA_matrix.tab", row.names = T, quote = F, col.names = T, sep = "\t")
                write.table(expr0, file = "2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_id.tab", row.names = T, quote = F, col.names = T, sep = "\t")
                write.table(sample_annot, file = "2_sample_annot.tab", row.names = F, quote = F, col.names = T, sep = "\t")}
        if (my.mode == 2){
                write.table(sample_annot, file = "2_sample_annot_corregido.tab", row.names = F, quote = F, col.names = T, sep = "\t")}}
if (my.mode == 1){
        sample_data <- read.delim("2_sample_annot.tab", sep = "\t")}
if (my.mode == 2){
        sample_data <- read.delim("2_sample_annot_corregido.tab", sep = "\t")}

#######################################################################################################################################
sample_data <- as.data.frame(sample_data)
sample_data <- sample_data[c("SampleName", "sample_type", "BRCA_Subtype_PAM50")]
colnames(sample_data) <- c("column", "Condition", "Subtype")
sample_data <- sample_data[is.na(sample_data$Subtype) != T,]
sample_data$Condition[sample_data$Condition == "NT"] <- "Normal"
sample_data$Condition[sample_data$Condition == "TP"] <- "Tumoral"
sample_data <- as.data.frame(sample_data)
#######################################################################################################################################
for (RNAs in c("mRNAs", "miRNAs")) {
        print(RNAs)
        Result_directory <- paste("Resultados_expresion_diferencial_", RNAs, "-gene_level/", sep = "")
        if (my.mode == 2){Result_directory <- paste("Corregido_resultados_expresion_diferencial_", RNAs, "-gene_level/", sep = "")}
        if (RNAs == "mRNAs") {
           expr0 <- read.delim("1_mRNAs_gene-level_Counts_TCGA_matrix.tab", sep = "\t")
           expr1 <- read.delim("2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_id.tab", sep = "\t")
        } else {
           expr0 <- read.delim("miRNAs_counts.tab", sep = "\t", row.names=1)
           expr1 <- read.delim("2_miRNAs_gene-level_TMM_TCGA_matrix.tab", sep = "\t")
        }
        colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
        colnames(expr0) <- substring(colnames(expr0), 1, 15)        
        matrix <- as.data.frame(expr0)       
        coldata <- sample_data
        cts <- matrix[coldata$column]        

        colnames(expr1) <- gsub("\\.", "-", colnames(expr1))
        colnames(expr1) <- substring(colnames(expr1), 1, 15)        
        matrix2 <- as.data.frame(expr1) 
        mdp_samples <- coldata
        colnames(mdp_samples) <- c("Sample", "Class", "Subtype")
        mdp_samples <- as.data.frame(mdp_samples)
        head(coldata$column)
        cts2 <- matrix2[coldata$column]
        print("antesdelerror?3")
        head(colnames(cts2))
        colnames(cts2) <- gsub("-","\\.",colnames(cts2))
        mdp_samples$Sample <- gsub("-","\\.",mdp_samples$Sample)
        print("antes de mdp")
        mdp(cts2, mdp_samples, control_lab = "Normal", file_name = "all_samples_", measure = "mean", directory = Result_directory) 

        dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
        keep <- rowSums(counts(dds)) > 1 ## Pre-Filtering
        dds <- dds[keep,]
        dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
        dds$Condition <- droplevels(dds$Condition)
        dds <- DESeq(dds, fitType="local", parallel = TRUE)        
        if (my.mode == 2){
                DESeq_norm <- counts(dds, normalized=T)
                write.table(DESeq_norm, file = paste0("Gene-level_",RNAs,"_DESeq_norm_matrix_WT.csv"), row.names = T, quote = F, col.names = T, sep = ",")}  
        vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
        pcaData <- plotPCA(vsd, intgroup=c("Condition", "Subtype"), returnData=TRUE)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Subtype)) +
          ggtitle(paste("All samples for ", RNAs, "\n diferential expression analysis-gene level", sep = "")) +
          theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
          geom_point(size=1) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          coord_fixed() +
          scale_color_manual(values=c("#56B4E9", "red"))
        ggsave(paste(Result_directory, "all_samples_", RNAs, "_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in")
        res <- results(dds, parallel = TRUE)
        pdf(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep= ""), height = 6.5, width = 10)
        DESeq2::plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in all samples", RNAs, "analysis\n", sep = " "), ylim=c(-10,10))
        abline(h=c(-1,1), col="dodgerblue", lwd=2)
        dev.off()
        bitmap <- pdf_render_page(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""), page = 1, dpi = 300)
        png::writePNG(bitmap, paste(Result_directory, "all_samples_", RNAs, "_MA_plot.png", sep = ""))
        unlink(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""))
        write.csv(res, file=paste(Result_directory, "all_samples_", RNAs, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))
        if (my.mode == 2) {          
                TCGA_files <- c("LumA", "LumB", "Her2", "Basal")
                for (File in TCGA_files){
                    coldata <- sample_data[sample_data$Subtype == File | sample_data$Condition == "Normal",]
                    cts <- matrix[coldata$column]
                    mdp_samples <- coldata
                    colnames(mdp_samples) <- c("Sample", "Class", "Subtype")
                    mdp_samples <- as.data.frame(mdp_samples)
                    cts2 <- matrix2[coldata$column]
                    colnames(cts2) <- gsub("-","\\.",colnames(cts2))
                    mdp_samples$Sample <- gsub("-","\\.",mdp_samples$Sample)

                    mdp(cts2, mdp_samples, control_lab = "Normal", file_name = paste(File, "_all_controls_", sep = ""), measure = "mean", directory = Result_directory)
                    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
                    keep <- rowSums(counts(dds)) > 1 ## Pre-Filtering
                    dds <- dds[keep,]
                    dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
                    dds$Condition <- droplevels(dds$Condition)
                    dds <- DESeq(dds, fitType="local", parallel = TRUE)
                    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
                    pcaData <- plotPCA(vsd, intgroup=c("Condition", "Subtype"), returnData=TRUE)
                    percentVar <- round(100 * attr(pcaData, "percentVar"))
                    p <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Subtype)) +
                      ggtitle(paste(File, "samples for", RNAs, "\n diferential expression analysis", sep = " ")) +
                      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
                      geom_point(size=1) +
                      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                     coord_fixed() +
                     scale_color_manual(values=c("#56B4E9", "red"))
                    ggsave(paste(Result_directory, File, "_samples_all_controls_", RNAs, "_PCA_plot.png", sep = ""), plot = p, width = 8.5, height = 8.5, dpi = 300, units = "in") #####################
                   res <- results(dds, parallel = TRUE)
                    pdf(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""), height = 6.5, width = 10)
                    DESeq2::plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "samples", RNAs, "analysis\n", sep = " "), ylim=c(-10,10))
                    abline(h=c(-1,1), col="dodgerblue", lwd=2)
                    dev.off()
                    bitmap <- pdf_render_page(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""), page = 1, dpi = 300)
                   png::writePNG(bitmap, paste(Result_directory, File, "_samples_all_controls_", RNAs, "_MA_plot.png", sep = "")) #####################
                    unlink(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""))
                    write.csv( res, file=paste(Result_directory, File, "_samples_all_controls_", RNAs,"_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))}}
        if (my.mode == 6){
                system(paste("tar -czvf", paste("resultados_", RNAs, "-gene_level.tar.gz", sep =""), Result_directory, sep = " "))}
        if (my.mode == 2){
                system(paste("tar -czvf", paste("corregido_resultados_", RNAs, "-gene_level.tar.gz", sep =""), Result_directory, sep = " "))}}
if (my.mode == 6){
        system(paste("scp -P 1313 ", "resultados*.tar.gz ",Taruca_adolfo_tesis, "Pareados/", sep = ""))}
if (my.mode == 2){
        system(paste("scp -P 1313 ", "corregido*.tar.gz ",Taruca_adolfo_tesis, "Pareados/", sep = ""))}
}
message_to_bot <- 'Script Finalizado:\n"Analisis de muestras y expresion diferencial"'
bot$sendMessage(chat_id, text = message_to_bot)
