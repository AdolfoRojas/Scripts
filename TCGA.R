library(TCGAbiolinks)
library(ggplot2)
library(dplyr)
setwd("C:/Users/adolf/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/adolf/Tesis")

#-----------------------------------------------------------------------------------------------------------------------------------------#

file <- "TCGA_ID_MAP_BRCA.csv"
Original_file<- read.table(file, header=T, sep = ",")
Original_file$patient <- as.factor(substr(Original_file$AliquotBarcode, start = 1 , stop = 12 ))
length(unique(Original_file$patient))
Original_file$TSS <- as.factor(substr(Original_file$AliquotBarcode, start = 6 , stop = 7 ))

#Clinical_data_BrCa <- GDCquery_clinic("TCGA-BRCA", type = "clinical", save.csv = T)
Clinical_data_BrCa <- read.table("TCGA-BRCA_clinical.csv", header=T, sep = ",")
Reduced_Clinical_data_BrCa <- Clinical_data_BrCa[c(1, 12, 46, 45, 44)]
names(Reduced_Clinical_data_BrCa)[names(Reduced_Clinical_data_BrCa) == "submitter_id"] <- "patient"

Tissue_Source_Site <- read.table("tissueSourceSite_BrCa.tab", header=T, sep = "\t")
names(Tissue_Source_Site)[names(Tissue_Source_Site) == "TSS.Code"] <- "TSS"
names(Tissue_Source_Site)[names(Tissue_Source_Site) == "Source.Site"] <- "Tissue_Source_Site"
Reduced_Tissue_Source_Site <- Tissue_Source_Site[c(1, 2)]

bar <- Original_file$AliquotBarcode
Primary_tumors <- TCGAquery_SampleTypes(bar,"TP")
Normal_solid_tissue <- TCGAquery_SampleTypes(bar,"NT")
Blood_Derived_Normal <- TCGAquery_SampleTypes(bar,"NBM")
Metastatic <- TCGAquery_SampleTypes(bar,"TM")

Sample_type_Primary_tumors <- data.frame(AliquotBarcode = Primary_tumors, sample_type = rep("TP"))
Sample_type_Normal_solid_tissue <- data.frame(AliquotBarcode = Normal_solid_tissue, sample_type = rep("NT"))
Sample_type_Metastatic <- data.frame(AliquotBarcode = Metastatic, sample_type = rep("TM")) 

sample_Type <- rbind(Sample_type_Primary_tumors, Sample_type_Metastatic, Sample_type_Normal_solid_tissue)

Original_file <- as.data.frame(merge(Original_file, Reduced_Tissue_Source_Site, by = "TSS")) #             Problema
length(unique(Original_file$patient))
Original_file <- as.data.frame(merge(Original_file, sample_Type, by = "AliquotBarcode"))
length(unique(Original_file$patient))
Original_file <- as.data.frame(merge(Original_file, Reduced_Clinical_data_BrCa, by = "patient"))
length(unique(Original_file$patient))

Reduced_Original_file <- Original_file[c(1,2,4,7:12)]
Reduced_Original_file$CGHubAnalysisID <- paste(Reduced_Original_file$sample_type, Reduced_Original_file$CGHubAnalysisID, sep = "_")

Subtype_BCa <- TCGAquery_subtype(tumor = "brca") #Retrieve molecular subtypes for a given tumor
Reduced_Subtype_BCa <- Subtype_BCa[c(1,9,12)]

df <- as.data.frame(merge(Reduced_Original_file, Reduced_Subtype_BCa, by = "patient"))
length(unique(df$patient))
df <- df[c(1:3,9,5,6,10,11,4,7,8)]
df <- df[ which(df$gender=='female'), ]

write.table(df, sep = "\t",
            file = "Samples_Subtype_BCa.tsv", 
            row.names = F, quote = F, col.names = T)
remove(df)
remove(Original_file)
remove(Reduced_Subtype_BCa)
remove(Reduced_Original_file)
remove(Subtype_BCa)

#--------------------------------------------------------------------------------------------------------------------------------------------------#

file <- "Samples_Subtype_BCa.tsv"
df<- read.table(file, header=T, sep = "\t")
df$BRCA_Subtype_PAM50 <- as.factor(df$BRCA_Subtype_PAM50)

# Cuantificar subtipos
Subtype_BCa_df <- data.frame(Subtype = c("Her2", "Basal", "LumB", "LumA", "Normal"), Samples =c(nrow(subset(df, grepl("Her2",df$BRCA_Subtype_PAM50))),
                      nrow(subset(df, grepl("Basal",df$BRCA_Subtype_PAM50))),
                      nrow(subset(df, grepl("LumB",df$BRCA_Subtype_PAM50))),
                      nrow(subset(df, grepl("LumA",df$BRCA_Subtype_PAM50))),
                      nrow(subset(df, grepl("Normal",df$BRCA_Subtype_PAM50)))))

Patient_samples <- as.data.frame(table(df$patient, dnn = list("Patient")), responseName = "Samples") # Ver individuos con mas de una muestra


# Cuentas 

Count_file <- "head_TCGA_BRCA_tpm.tsv"
Counts <- read.table(Count_file, header=T, sep = "\t")
Counts[] <- lapply(Counts[], factor) 

LuminalA <- df[ which(df$BRCA_Subtype_PAM50=='LumA'), ]
LuminalA_Cols <- c()
for (AnalysisID in LuminalA$CGHubAnalysisID) {
  names(Counts)[grep(substr(AnalysisID, start = 4 , stop = 11 ), names(Counts))] <- AnalysisID
  loop <- (grep(AnalysisID, names(Counts)))
  LuminalA_Cols <- append(LuminalA_Cols, loop)
}
LuminalA_counts <- Counts[LuminalA_Cols]
LuminalA_counts <- cbind(Counts$X, LuminalA_counts)
names(LuminalA_counts)[names(LuminalA_counts) == "Counts$X"] <- "Transcripts"
Number_LumA <- data.frame(Total = nrow(subset(LuminalA, grepl("LumA",LuminalA$BRCA_Subtype_PAM50))), 	PRIMARY_SOLID_TUMOR = nrow(subset(LuminalA, grepl("TP",LuminalA$sample_type))), 	Metastatic = nrow(subset(LuminalA, grepl("TM",LuminalA$sample_type))), Solid_Tissue_Normal = nrow(subset(LuminalA, grepl("NT",LuminalA$sample_type))), Unique_Patients =  length(unique(LuminalA$patient)), row.names = "Luminal A")
write.table(LuminalA_counts, sep = "\t",
            file = "LuminalA_counts.tab", 
            row.names = F, quote = F, col.names = T)
remove(LuminalA_counts)
remove(LuminalA)

LuminalB <- df[ which(df$BRCA_Subtype_PAM50=='LumB'), ]
LuminalB_Cols <- c()
for (AnalysisID in LuminalB$CGHubAnalysisID) {
  names(Counts)[grep(substr(AnalysisID, start = 4 , stop = 11 ), names(Counts))] <- AnalysisID
  loop <- (grep(AnalysisID, names(Counts)))
  LuminalB_Cols <- append(LuminalB_Cols, loop)
}
LuminalB_counts <- Counts[LuminalB_Cols]
LuminalB_counts <- cbind(Counts$X, LuminalB_counts)
names(LuminalB_counts)[names(LuminalB_counts) == "Counts$X"] <- "Transcripts"
Number_LumB <- data.frame(Total = nrow(subset(LuminalB, grepl("LumB",LuminalB$BRCA_Subtype_PAM50))), 	PRIMARY_SOLID_TUMOR = nrow(subset(LuminalB, grepl("TP",LuminalB$sample_type))), 	Metastatic = nrow(subset(LuminalB, grepl("TM",LuminalB$sample_type))), Solid_Tissue_Normal = nrow(subset(LuminalB, grepl("NT",LuminalB$sample_type))), Unique_Patients =  length(unique(LuminalB$patient)), row.names = "Luminal B")
write.table(LuminalB_counts, sep = "\t",
            file = "LuminalB_counts.tab", 
            row.names = F, quote = F, col.names = T)
remove(LuminalB_counts)
remove(LuminalB)

Her2 <- df[ which(df$BRCA_Subtype_PAM50=='Her2'), ]
Her2_Cols <- c()
for (AnalysisID in Her2$CGHubAnalysisID) {
  names(Counts)[grep(substr(AnalysisID, start = 4 , stop = 11 ), names(Counts))] <- AnalysisID
  loop <- (grep(AnalysisID, names(Counts)))
  Her2_Cols <- append(Her2_Cols, loop)
}
Her2_counts <- Counts[Her2_Cols]
Her2_counts <- cbind(Counts$X, Her2_counts)
names(Her2_counts)[names(Her2_counts) == "Counts$X"] <- "Transcripts"
Number_Her2 <- data.frame(Total = nrow(subset(Her2, grepl("Her2",Her2$BRCA_Subtype_PAM50))), 	PRIMARY_SOLID_TUMOR = nrow(subset(Her2, grepl("TP",Her2$sample_type))), 	Metastatic = nrow(subset(Her2, grepl("TM",Her2$sample_type))), Solid_Tissue_Normal = nrow(subset(Her2, grepl("NT",Her2$sample_type))), Unique_Patients =  length(unique(Her2$patient)), row.names = "Her2")
write.table(Her2_counts, sep = "\t",
            file = "Her2_counts.tab", 
            row.names = F, quote = F, col.names = T)
remove(Her2_counts)
remove(Her2)

Basal <- df[ which(df$BRCA_Subtype_PAM50=='Basal'), ]
Basal_Cols <- c()
for (AnalysisID in Basal$CGHubAnalysisID) {
  names(Counts)[grep(substr(AnalysisID, start = 4 , stop = 11 ), names(Counts))] <- AnalysisID
  loop <- (grep(AnalysisID, names(Counts)))
  Basal_Cols <- append(Basal_Cols, loop)
}
Basal_counts <- Counts[Basal_Cols]
Basal_counts <- cbind(Counts$X, Basal_counts)
names(Basal_counts)[names(Basal_counts) == "Counts$X"] <- "Transcripts"
Number_Basal <- data.frame(Total = nrow(subset(Basal, grepl("Basal",Basal$BRCA_Subtype_PAM50))), 	PRIMARY_SOLID_TUMOR = nrow(subset(Basal, grepl("TP",Basal$sample_type))), 	Metastatic = nrow(subset(Basal, grepl("TM",Basal$sample_type))), Solid_Tissue_Normal = nrow(subset(Basal, grepl("NT",Basal$sample_type))), Unique_Patients =  length(unique(Basal$patient)), row.names = "Basal")
write.table(Basal_counts, sep = "\t",
            file = "Basal_counts.tab", 
            row.names = F, quote = F, col.names = T)

Normal <- df[ which(df$BRCA_Subtype_PAM50=='Normal'), ]
Normal_Cols <- c()
for (AnalysisID in Normal$CGHubAnalysisID) {
  names(Counts)[grep(substr(AnalysisID, start = 4 , stop = 11 ), names(Counts))] <- AnalysisID
  loop <- (grep(AnalysisID, names(Counts)))
  Normal_Cols <- append(Normal_Cols, loop)
}
Basal_counts <- Counts[Basal_Cols]
Basal_counts <- cbind(Counts$X, Basal_counts)
names(Basal_counts)[names(Basal_counts) == "Counts$X"] <- "Transcripts"
Number_Basal <- data.frame(Total = nrow(subset(Basal, grepl("Basal",Basal$BRCA_Subtype_PAM50))), 	PRIMARY_SOLID_TUMOR = nrow(subset(Basal, grepl("TP",Basal$sample_type))), 	Metastatic = nrow(subset(Basal, grepl("TM",Basal$sample_type))), Solid_Tissue_Normal = nrow(subset(Basal, grepl("NT",Basal$sample_type))), Unique_Patients =  length(unique(Basal$patient)), row.names = "Basal")
write.table(Basal_counts, sep = "\t",
            file = "Basal_counts.tab", 
            row.names = F, quote = F, col.names = T)

# Cuantificar subtipos

Subtype_BCa_df <- rbind(Number_LumA, Number_LumB, Number_Her2, Number_Basal) 
write.table(Subtype_BCa_df, sep = "\t",
            file = "Subtype_BCa_df.info", 
            row.names = F, quote = F, col.names = T)
remove(Subtype_BCa_df)
remove(Basal_counts)
remove(Basal)

remove(df)
remove(Counts)



#---------------------------------------------------------------------------------------------------------------------------------------------#
file3 <- "Tab_counts_procesado.tab"
Counts2 <- as.data.frame(read.table(file3, header=T, sep = "\t"))


summary(Original_file$patient)



# Retrieve multiple tissue types  NOT FROM THE SAME PATIENTS
SS <- TCGAquery_SampleTypes(bar,c("TP","NB"))

# Retrieve multiple tissue types  FROM THE SAME PATIENTS
SSS <- TCGAquery_MatchedCoupledSampleTypes(bar,c("NT","TP"))

metadata <- colDataPrepare(dataSubt$patient)

metadata <- colDataPrepare(Original_file$AliquotBarcode)
summary <- getSampleFilesSummary("TCGA-BRCA") #Retrieve summary of files per sample in a project

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
results <- getResults(query)
Manifest<-getManifest(query)

query2 <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification",
                  experimental.strategy = "miRNA-Seq",
                  barcode = "TCGA-AO-A0JF-01A-11R-A057-13")

results2 <- getResults(query2)


# You can define a list of samples to query and download providing relative TCGA barcodes.
listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07")

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts",
                  barcode = "TCGA-E9-A1NG-11A-52R-A14M-07")

 # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query2)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")


query <- GDCquery(project = "TCGA-ACC",
                  data.category =  "Copy number variation",
                  legacy = TRUE,
                  file.type = "hg19.seg",
                  barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01", "TCGA-OR-A5LJ-10A-01D-A29K-01"))
# data will be saved in  GDCdata/TCGA-ACC/legacy/Copy_number_variation/Copy_number_segmentation
GDCdownload(query, method = "client")





data <- data.frame(
  Group=row.names(Subtype_BCa_df),
  value=Subtype_BCa_df$Total
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(Group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=Group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  #theme(legend.position="none") +
  
  geom_text(aes(y = ypos, label = ypos+Group), color = "white", size=6) +
  scale_fill_brewer(palette="Set1")

