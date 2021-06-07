#!/usr/bin/env Rscript
library(pdftools)
library(DESeq2)
library(ggplot2)
library(mdp)
library(BiocParallel)
register(MulticoreParam(60))
Taruca_adolfo_tesis <- "adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/"
system("mkdir -p Resultados_expresion_diferencial_mRNAs-gene_level")
system("mkdir -p Resultados_expresion_diferencial_miRNAs-gene_level")
system("mkdir -p Corregido_resultados_expresion_diferencial_miRNAs-gene_level")
system("mkdir -p Corregido_resultados_expresion_diferencial_mRNAs-gene_level")
my.mode <- 2
#while (my.mode!=1 | my.mode != 2){
#my.mode <- readline(prompt="Enter mode: \n(1 = Descargar y analizar)\n(2 = Ajustar y re-analizar)")
print(my.mode)
#######################################################################################################################################
# Contiene el script tesis/3_co-expresion_analysis/1_Get_mRNAs_gene-level_FPKM_TCGA_matrix.R
#######################################################################################################################################
if (my.mode == 1) {
        require(TCGAbiolinks)
        library(SummarizedExperiment)
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
 #              Tratado de muestras  contiene el script tesis/3_co-expresion_analysis/2_prepare_CEMiTool_data.R
#######################################################################################################################################
expr0 <- read.delim("1_mRNAs_gene-level_Counts_TCGA_matrix.tab")
library(EnsDb.Hsapiens.v86)
ensembl.genes <- rownames(expr0)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
for (duplicado in levels(as.factor(geneIDs1[duplicated(geneIDs1$SYMBOL),]$SYMBOL))){
        count = 0
        for (elementos in rownames(geneIDs1[geneIDs1$SYMBOL == duplicado,])){
                count = count + 1
                geneIDs1$SYMBOL[as.numeric(elementos)] <- paste(duplicado, count, sep = "_")}}
for (ensembl_gene_id in rownames(expr0)){
        rownames(expr0)[rownames(expr0)== ensembl_gene_id] <- geneIDs1[geneIDs1$GENEID == ensembl_gene_id,]$SYMBOL}
if (my.mode == 1) {
        write.table(geneIDs1, file = "2_gene_ID_to_gene_symbol.tab", row.names = F, quote = F, col.names = T, sep = "\t")}
#######################################################################################################################################
sample_annot <- ""
if (my.mode == 1){
        sample_annot <- read.delim("Samples_Subtype_BRCA.tab", sep = "\t")
        sample_annot <- sample_annot[sample_annot$BRCA_Subtype_PAM50 != "Normal",]
        sample_annot <- sample_annot[sample_annot$BRCA_Subtype_PAM50 != "NA",]
        sample_annot$AliquotBarcode <- substring(sample_annot$AliquotBarcode, 1, 15)
        colnames(sample_annot)[colnames(sample_annot)=="AliquotBarcode"] <- "SampleName"}
if (my.mode == 2){
        sample_annot <- read.delim("2_sample_annot.tab", sep = "\t") #añadir muestras outliers        
        #sample_annot <- sample_annot[duplicated(sample_annot$patient) | duplicated(sample_annot$patient,fromLast=T),]
        #mdp <- read.delim("Resultados_expresion_diferencial_mRNAs-gene_level/all_samples_sample_scores.tsv", sep = " ")
        #mdp$allgenes.Sample <- gsub("\\.", "-", mdp$perturbedgenes.Sample)
        #normal_mdp <- mdp[mdp$allgenes.Class == "Normal",]
        #normal_mdp_disturbed <- tail(normal_mdp[order(normal_mdp$allgenes.Score),]$allgenes.Sample, 10)
        #sample_annot <- sample_annot[!(sample_annot$SampleName %in% normal_mdp_disturbed),]
        cluster_info <- read.delim("cluster_output.tab", sep = "\t", header=FALSE)
        sample_annot <- sample_annot[!(sample_annot$SampleName %in% cluster_info$V1),]}
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

if (identical(colnames(expr0), sample_annot$SampleName) == T){
        print("guardando tablas y procediendo con anotacion en gene symbol de interacciones proteina-proteina")
        if (my.mode == 1){
                write.table(expr0, file = "2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_symbol.tab", row.names = T, quote = F, col.names = T, sep = "\t")
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
        } else {
           expr0 <- read.delim("miRNAs_counts.tab", sep = "\t", row.names=1)
        }
        colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
        colnames(expr0) <- substring(colnames(expr0), 1, 15)        
        matrix <- as.data.frame(expr0)
        coldata <- sample_data
        cts <- matrix[coldata$column] + 1
        mdp_samples <- coldata
        colnames(mdp_samples) <- c("Sample", "Class", "Subtype")
        mdp_samples <- as.data.frame(mdp_samples)
        head(coldata$column)
        cts2 <- matrix[coldata$column]
        print("antesdelerror?3")
        head(colnames(cts2))
        colnames(cts2) <- gsub("-","\\.",colnames(cts2))
        mdp_samples$Sample <- gsub("-","\\.",mdp_samples$Sample)
        print("antes de mdp")
        mdp(cts2, mdp_samples, control_lab = "Normal", file_name = "all_samples_", measure = "mean", directory = Result_directory) 

        dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
        keep <- rowSums(counts(dds)) > length(cts) ## Pre-Filtering
        dds <- dds[keep,]
        dds$Condition <- factor(dds$Condition, levels = c("Normal","Tumoral"))
        dds$Condition <- droplevels(dds$Condition)
        dds <- DESeq(dds, fitType="local", parallel = TRUE)
        vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
        pcaData <- plotPCA(vsd, intgroup=c("Condition", "Subtype"), returnData=TRUE)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        if (my.mode == 1 & RNAs == "mRNAs"){
                normal_pca <- pcaData[pcaData$Condition == "Normal",]
                library(dbscan)
                df <- normal_pca[c("PC1","PC2")]
                pdf("Grafico.pdf", height = 6.5, width = 10)
                kNNdistplot(df, k = 3  )
                dev.off()
                cl<-dbscan(df,eps=4,MinPts = 3)
                df2 <- cbind(df,cl$cluster)
                names(df2)[3] <- "cluster"
                df2 <- rownames(df2[df2$cluster != 1,])
                pdf("Grafico2.pdf", height = 6.5, width = 10)
                hullplot(df,cl$cluster, main = "Convex cluster Hulls, eps= 4")
                dev.off()
                write.table(df2, file = "cluster_output.tab", row.names = F, quote = F, col.names = F, sep = "\t")}
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
        plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in all samples", RNAs, "analysis\n", sep = " "), ylim=c(-10,10))
        abline(h=c(-1,1), col="dodgerblue", lwd=2)
        dev.off()
        bitmap <- pdf_render_page(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""), page = 1, dpi = 300)
        png::writePNG(bitmap, paste(Result_directory, "all_samples_", RNAs, "_MA_plot.png", sep = ""))
        unlink(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""))
        write.csv(res, file=paste(Result_directory, "all_samples_", RNAs, "_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))

        TCGA_files <- c("LumA", "LumB", "Her2", "Basal")
        for (File in TCGA_files){
            coldata <- sample_data[sample_data$Subtype == File | sample_data$Condition == "Normal",]
            cts <- matrix[coldata$column] + 1
            mdp_samples <- coldata
            colnames(mdp_samples) <- c("Sample", "Class", "Subtype")
            mdp_samples <- as.data.frame(mdp_samples)
            cts2 <- matrix[coldata$column]
            colnames(cts2) <- gsub("-","\\.",colnames(cts2))
            mdp_samples$Sample <- gsub("-","\\.",mdp_samples$Sample)

            mdp(cts2, mdp_samples, control_lab = "Normal", file_name = paste(File, "_all_controls_", sep = ""), measure = "mean", directory = Result_directory)
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
            plotMA(res, main = paste(levels(dds$Condition)[1], "vs", levels(dds$Condition)[2], "Diferential Expresion in", File, "samples", RNAs, "analysis\n", sep = " "), ylim=c(-10,10))
            abline(h=c(-1,1), col="dodgerblue", lwd=2)
            dev.off()
            bitmap <- pdf_render_page(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""), page = 1, dpi = 300)
           png::writePNG(bitmap, paste(Result_directory, File, "_samples_all_controls_", RNAs, "_MA_plot.png", sep = "")) #####################
            unlink(paste("all_samples_", RNAs, "_MA_plotx.pdf", sep = ""))
            write.csv( res, file=paste(Result_directory, File, "_samples_all_controls_", RNAs,"_", levels(dds$Condition)[1], "_vs_", levels(dds$Condition)[2],"_DE.tab", sep= ""))#######################
        }
        if (my.mode == 1){
                system(paste("tar -czvf", paste("resultados_", RNAs, "-gene_level.tar.gz", sep =""), Result_directory, sep = " "))}
        if (my.mode == 2){
                system(paste("tar -czvf", paste("corregido_resultados_", RNAs, "-gene_level.tar.gz", sep =""), Result_directory, sep = " "))}}
if (my.mode == 1){
        system(paste("scp ", "resultados*.tar.gz ",Taruca_adolfo_tesis, sep = ""))}
if (my.mode == 2){
        system(paste("scp ", "corregido*.tar.gz ",Taruca_adolfo_tesis, sep = ""))}
