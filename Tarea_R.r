#!/usr/bin/env Rscript
#library(pdftools)
#library(DESeq2)
#library(ggplot2)
#library(mdp)
#library(BiocParallel)
#library(telegram.bot)
#library(dplyr)
#library(SummarizedExperiment)

require(TCGAbiolinks)  

CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("GDC/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","miRNA_gene_quantification",".rda")
query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE)

samplesDown.miR <- getResults(query.miR,cols=c("cases"))

dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
# problema diferente numero de muestras de Tumor primario  y Tejido solido normal, se necesitan muestras pareadas (cada paciente con su respectivo control normal y muestra tumoral)
# PISTA = los barcodes de las muestras, etiquetan a los pacientes hasta el 12vo caracter
# Solucion sugerida 
df <- as.data.frame(dataSmTP.miR)
df2 <- as.data.frame(dataSmNT.miR)

df$paciente <- substr(df$dataSmTP.miR, 1, 12)
df[duplicated(df$paciente),]
df2$paciente <- substr(df2$dataSmNT.miR, 1, 12)
df2[duplicated(df2$paciente),]
paired_samples <- merge(df,df2,by= "paciente")

# Cuantos pacientes tienen muestras pareadas?   length(unique(paired_samples$paciente))
# Hay pacientes con mas de una muestra control??    df2[duplicated(df2$paciente),]
# Cuantos pacientes tienen mas de una muestra tumoral?  paired_samples$paciente[duplicated(paired_samples$paciente)]

# Mantener solo una muestra tumoral por paciente y reasignar los vectores dataSmTP.miR y dataSmNT.miR para continuar con la creacion de la matriz de cuentas
paired_samples <- paired_samples[!duplicated(paired_samples$paciente),]
dataSmTP.miR <- paired_samples$dataSmTP.miR
dataSmNT.miR <- paired_samples$dataSmNT.miR

queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
GDCdownload(query = queryDown.miR, directory = DataDirectory) ## aprox 10 megas de descarga, reducir numero de pacientes si esta lento



dataAssy.miR <- GDCprepare(query = queryDown.miR, save = TRUE, save.filename = FileNameData, summarizedExperiment = TRUE, directory =DataDirectory)
rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID
read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))] # using read_count's data
matrix <- dataAssy.miR[, read_countData]
colnames(matrix) <- gsub("read_count_","", colnames(matrix))
names <- rownames(dataAssy.miR)
matrix<- cbind(names,matrix)      
write.table(matrix, sep = "\t", file = "miRNAs_counts.tab", row.names = FALSE, quote = FALSE, col.names = T)