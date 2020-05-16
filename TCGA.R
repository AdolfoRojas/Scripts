library(TCGAbiolinks)
library(ggplot2)
require(moonBook)
require(webr)
library(pdftools)
setwd("../Tesis")
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

Original_file <- as.data.frame(merge(Original_file, sample_Type, by = "AliquotBarcode"))
length(unique(Original_file$patient))
Original_file <- as.data.frame(merge(Original_file, Reduced_Clinical_data_BrCa, by = "patient")) # 1256 M, 1095 P.

Reduced_Original_file <- Original_file[c(1,2,4,7:12)]
Reduced_Original_file$CGHubAnalysisID <- paste(Reduced_Original_file$sample_type, Reduced_Original_file$CGHubAnalysisID, sep = "_")

Subtype_BCa <- TCGAquery_subtype(tumor = "brca") #Retrieve molecular subtypes for a given tumor
Reduced_Subtype_BCa <- Subtype_BCa[c(1,9,12)]

df <- as.data.frame(merge(Reduced_Original_file, Reduced_Subtype_BCa, by = "patient")) # 1243M, 1083 P
df <- df[c(1:3,9,5,6,10,11,4,7,8)]
df <- df[ which(df$gender=='female'), ] # 1241 M, 1082 P
df <- df[ which(df$sample_type !='TM'), ] # 1234 M, 1082 P
length(unique(df$patient))

write.table(df, sep = "\t",
            file = "Samples_Subtype_BCa.tsv", 
            row.names = F, quote = F, col.names = T)
remove(df)
remove(Original_file)
remove(Reduced_Subtype_BCa)
remove(Reduced_Original_file)
remove(Subtype_BCa)

#--------------------------------------------------------------------------------------------------------------------------------------------------#

file2 <- "Samples_Subtype_BCa.tsv"
df<- read.table(file2, header=T, sep = "\t")
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
library(ggplot2)
require(moonBook)
require(webr)
library(pdftools)
setwd("../Tesis")
file2 <- "Samples_Subtype_BCa.tsv"

df_plot1 <- read.table(file2, header=T, sep = "\t")
df_plot1 <- df_plot1[ which(df_plot1$BRCA_Subtype_PAM50!='Normal'), ] # 1192 M, 1042 P
write.table(df_plot1, sep = "\t",
            file = "df_plot1", 
            row.names = F, quote = F, col.names = T)
df_plot1 <- read.table("df_plot1", header=T, sep = "\t")
levels(df_plot1$BRCA_Subtype_PAM50)[levels(df_plot1$BRCA_Subtype_PAM50) == "LumA"] <- "Luminal A"
levels(df_plot1$BRCA_Subtype_PAM50)[levels(df_plot1$BRCA_Subtype_PAM50) == "LumB"] <- "Luminal B"
levels(df_plot1$BRCA_Subtype_PAM50)[levels(df_plot1$BRCA_Subtype_PAM50) == "Her2"] <- "HER2"
length(unique(df_plot1$patient))
pdf("df_plot1.pdf", height = 8.5, width = 8.5)
PieDonut(df_plot1,aes(pies=BRCA_Subtype_PAM50,donuts=sample_type), title="Subtipos tumorales", ratioByGroup=FALSE, showPieName=FALSE, r0=0, labelpositionThreshold = 0.5, titlesize = 15, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.001),  pieLabelSize = 5, donutLabelSize = 5)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, "df_plot1.png")
unlink("df_plot1.pdf")
unlink("df_plot1")
remove(df_plot1)

#---------------------------------------------------------------------------------------------------------------------------------------------#

df_plot2 <- read.table(file2, header=T, sep = "\t")
pdf("df_plot2.pdf", height = 8.5, width = 8.5)
PieDonut(df_plot2,aes(pies=gender,donuts=race),showPieName=FALSE, title ="Distribución según raza", r0=0.2, r1 = 0.2, r2 = 0.4, start=5.9*pi/2, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.0001), labelpositionThreshold=0.5, titlesize = 15, donutLabelSize = 6, addPieLabel = FALSE, pieLabelSize = 0, showRatioPie = F)
dev.off()
df_plot2 <- df_plot2[ which(df_plot2$sample_type=='NT'), ] # 112 pacientes 112 Muestras iguales, son 113 elementos ya que se tiene una muestra con 2 analisis
write.table(df_plot2, sep = "\t",
            file = "df_plot2", 
            row.names = F, quote = F, col.names = T)
df_plot2 <- read.table("df_plot2", header=T, sep = "\t")
bitmap <- pdf_render_page("df_plot2.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, "df_plot2.png")
unlink("df_plot2.pdf")
unlink("df_plot2")
remove(df_plot2)
