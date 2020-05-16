#!/usr/bin/env Rscript
setwd("C:/Users/adolf/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/adolf/Tesis")
Subtype_file <- "Samples_Subtype_BCa.tsv"
df<- read.table(Subtype_file, header=T, sep = "\t")
df$BRCA_Subtype_PAM50 <- as.factor(df$BRCA_Subtype_PAM50)


Patient_samples <- as.data.frame(table(df$patient, dnn = list("Patient")), responseName = "Samples") # ver individuos con mas de una muestra
write.table(Patient_samples, sep = "\t",
            file = "Patient_samples.info", 
            row.names = F, quote = F, col.names = T)
remove(Patient_samples)


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
#ftgty