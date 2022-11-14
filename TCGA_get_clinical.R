library(TCGAbiolinks)
#-----------------------------------------------------------------------------------------------------------------------------------------#
Cancer <- "LUAD"            ##Enter Cancer study abbreviation

#if (file.exists(paste("TCGA_ID_MAP_",  Cancer, ".csv", sep = "")) == F){
#  ID_file <- read.table("TCGA_ID_MAP.csv", header=T, sep = ",") #archivo con las muestras en el serve#r
#  Cancer_IDs_file <- ID_file[which(ID_file$Disease == Cancer), ]# extrae las muestras asociadas al cancer desead#o
#  write.table(Cancer_IDs_file, sep = ",", file = paste("TCGA_ID_MAP_",  Cancer, ".csv", sep = ""), row.names = F, quote = F, col#.names = T)
#  remove(ID_file)
#  remove(Cancer_IDs_file)
#}

#Original_file <- read.table(paste("TCGA_ID_MAP_",  Cancer, ".csv", sep = ""), header=T, sep = ",")
#Original_file$patient <- as.factor(substr(Original_file$AliquotBarcode, start = 1 , stop = 12)) # recorta el barcode de la muestra #para obtener la etiqueta del paciente
#length(unique(Original_file$patient))
#Original_file$TSS <- as.factor(substr(Original_file$AliquotBarcode, start = 6 , stop = 7)) # extrae los datos finales del barcode #asociado al TSS

if (file.exists(paste("TCGA-",  Cancer, "_clinical.csv", sep = "")) == F){
  GDCquery_clinic(paste("TCGA-", Cancer, sep = ""), type = "clinical", save.csv = T) # Descarga las caracteristicas clinicas de las muestras
}

Clinical_data <- read.table(paste("TCGA-",  Cancer, "_clinical.csv", sep = ""), header=T, sep = ",")
Reduced_Clinical_data <- Clinical_data[c("submitter_id", "primary_diagnosis", "ethnicity", "race", "gender","age_at_diagnosis")] # Caracteristicas Clinicaas de Interes
length(unique(Reduced_Clinical_data$submitter_id))


# Obtener tipo de muestra, ejemplo TP=tumor primario, NT = tejido normal (generalmente cercano a la zona del tumor, segun lei)
bar <- Original_file$AliquotBarcode
Primary_tumors <- TCGAquery_SampleTypes(bar,"TP")
Normal_solid_tissue <- TCGAquery_SampleTypes(bar,"NT")
Blood_Derived_Normal <- TCGAquery_SampleTypes(bar,"NBM")
Metastatic <- TCGAquery_SampleTypes(bar,"TM")

Sample_type_Primary_tumors <- data.frame(AliquotBarcode = Primary_tumors, sample_type = rep("TP"))
Sample_type_Normal_solid_tissue <- data.frame(AliquotBarcode = Normal_solid_tissue, sample_type = rep("NT"))
Sample_type_Metastatic <- data.frame(AliquotBarcode = Metastatic, sample_type = rep("TM"))
sample_Type <- rbind(Sample_type_Primary_tumors, Sample_type_Metastatic, Sample_type_Normal_solid_tissue)


# Rastrea lugar de procedencia de las muestras, centro de secuenciacion
projects <- getGDCprojects()
TSS_Cancer <- projects$name[grep(Cancer, projects$tumor)]
Tissue_Source_Site <- read.delim("tissueSourceSite.tsv", header=T, sep = "\t")
names(Tissue_Source_Site)[names(Tissue_Source_Site) == "TSS.Code"] <- "TSS"
names(Tissue_Source_Site)[names(Tissue_Source_Site) == "Source.Site"] <- "Tissue_Source_Site"
Tissue_Source_Site <- Tissue_Source_Site[grep(TSS_Cancer, Tissue_Source_Site$Study.Name, ignore.case = T),]
Reduced_Tissue_Source_Site <- Tissue_Source_Site[c(1, 2)]


# Cruzar tablas
Original_file <- as.data.frame(merge(Original_file, Reduced_Tissue_Source_Site, by = "TSS"))
length(unique(Original_file$patient)) # Agrego esto para revisar si hay perdidas de muestras y pacientes al cruzar
Original_file <- as.data.frame(merge(Original_file, sample_Type, by = "AliquotBarcode"))
length(unique(Original_file$patient))
Original_file <- as.data.frame(merge(Original_file, Reduced_Clinical_data, by.x = "patient", by.y = "submitter_id"))
length(unique(Original_file$patient))
Reduced_Original_file <- Original_file[c(1,2,4,7:13)]
Reduced_Original_file$CGHubAnalysisID <- paste(Reduced_Original_file$sample_type, Reduced_Original_file$CGHubAnalysisID, sep = "_")

# Subtipo tumoral, entrega datos del tumor y algunos marcadores moleculares

Subtype <- TCGAquery_subtype(tumor = Cancer) #
Reduced_Subtype <- Subtype[c(1,9,12)] #cambiar por columnas de interes

length(unique(Reduced_Subtype$patient))
length(unique(Original_file$patient))

df <- as.data.frame(merge(Reduced_Original_file, Reduced_Subtype, by = "patient"))
df <- df[c(1:3,9,5,6,10,11,12,4,7,8)]
length(unique(df$patient))
df <- df[ which(df$gender=='female'), ]
df <- df[ which(df$sample_type !='TM'), ]
length(unique(df$patient))

write.table(df, sep = "\t", file = paste("Samples_Subtype_",  Cancer, ".tsv", sep = ""), row.names = F, quote = F, col.names = T)
remove(df)
remove(Original_file)
remove(Reduced_Subtype)
remove(Reduced_Original_file)
remove(Subtype)


        require(TCGAbiolinks)        
        CancerProject <- "TCGA-LUAD"
        DataDirectory <- paste0("GDC/",gsub("-","_",CancerProject))
        FileNameData <- paste0(DataDirectory, "_","mRNA_gene_quantification",".rda")
        query.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", legacy = FALSE) #HTSeq - Counts
        samplesDown.miR <- getResults(query.miR,cols=c("cases"))
        dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
        dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
        queryDown.miR <- GDCquery(project = CancerProject, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", legacy = FALSE, barcode = c(dataSmTP.miR, dataSmNT.miR))
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
        write.table(matrix, sep = "\t", file = "miRNAs_counts.tab", row.names = F, quote = F, col.names = T) 