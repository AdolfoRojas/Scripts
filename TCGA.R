library(TCGAbiolinks)
library(ggplot2)
require(moonBook)
require(webr)
library(pdftools)
setwd("../Tesis")
#-----------------------------------------------------------------------------------------------------------------------------------------#

file <- "TCGA_ID_MAP_BRCA.csv"
file.exists("TCGA_ID_MAP_BRCA.csv")
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
PieDonut(df_plot1,aes(pies=BRCA_Subtype_PAM50,donuts=sample_type), title="Subtipos tumorales", ratioByGroup=T, showPieName=FALSE, r0=0, labelpositionThreshold = 0.5, titlesize = 15, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.001),  pieLabelSize = 5, donutLabelSize = 5)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, "df_plot1.png")
unlink("df_plot1.pdf")
unlink("df_plot1")
remove(df_plot1)

#---------------------------------------------------------------------------------------------------------------------------------------------#

df_plot2 <- read.table(file2, header=T, sep = "\t")
pdf("df_plot2.pdf", height = 8.5, width = 8.5)
PieDonut(df_plot2,aes(pies=gender,donuts=race),showPieName=FALSE, title ="Distribuci�n seg�n raza", r0=0.2, r1 = 0.2, r2 = 0.4, start=5.9*pi/2, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.0001), labelpositionThreshold=0.5, titlesize = 15, donutLabelSize = 6, addPieLabel = FALSE, pieLabelSize = 0, showRatioPie = F)
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

#---------------------------------------------------------------------------------------------------------------------------------------------#

file2 <- "Samples_Subtype_BCa.tsv"
df_plot3 <- read.table(file2, header=T, sep = "\t")
data <- data.frame(individual=character(),
                   group=character(), 
                   value=character(), 
                   stringsAsFactors=FALSE) 

for (Race in levels(df_plot3$race)) {
  df_loop1 <- df_plot3[ which(df_plot3$race==Race), ]
  for (Etnia in levels(df_loop1$ethnicity)) {
    df_loop2 <- df_loop1[ which(df_loop1$ethnicity==Etnia), ]
    freq <- nrow(df_loop2)
    data = rbind(data, data.frame(Etnia,Race,freq, fix.empty.names = F))
  }
  for (Type_sample in levels(df_loop1$sample_type)) {
    df_loop2 <- df_loop1[ which(df_loop1$sample_type==Type_sample), ]
    freq <- nrow(df_loop2)
    data = rbind(data, data.frame(Type_sample,Race,freq, fix.empty.names = F))
  }
}
Total = data.frame(individual=character())
for (n in levels(df_plot3$ethnicity)) {
  df_loop3 <- df_plot3[ which(df_plot3$ethnicity==n), ]
  Total = rbind(Total, data.frame(n, nrow(df_loop3), fix.empty.names = F))
}
for (n in levels(df_plot3$sample_type)) {
  df_loop3 <- df_plot3[ which(df_plot3$sample_type==n), ]
  Total = rbind(Total, data.frame(n, nrow(df_loop3), fix.empty.names = F))
}
colnames(Total) <- c("Etnia o tipo de muestra", "n")
colnames(data) <- c("Etnia o tipo de muestra", "Raza", "Proporci�n")
data <- as.data.frame(merge(data, Total, by = "Etnia o tipo de muestra"))

pdf("df_plot3.pdf", height = 8.5, width = 11)
ggplot(data, aes(fill=Raza, y=Proporci�n, x=`Etnia o tipo de muestra`)) + 
  geom_bar(position="fill", stat="identity")+ 
  ggtitle("Distribuci�n de Raza, seg�n etnia y tipo de muestra") +
  geom_text(data=data, aes(x=`Etnia o tipo de muestra`, y= rep(1.02, length(`Etnia o tipo de muestra`)), label=paste("n =", as.factor(n), sep = " "))) 
dev.off()
bitmap <- pdf_render_page("df_plot3.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, "df_plot3.png")
unlink("df_plot3.pdf")
unlink("df_plot3")

#---------------------------------------------------------------------------------------------------------------------------------------------#

file2 <- "Samples_Subtype_BCa.tsv"
df_plot3 <- read.table(file2, header=T, sep = "\t")
data <- data.frame(individual=character(),
                   group=character(), 
                   value=character(), 
                   stringsAsFactors=FALSE) 

for (Race in levels(df_plot3$primary_diagnosis)) {
  df_loop1 <- df_plot3[ which(df_plot3$primary_diagnosis==Race), ]
  for (Etnia in levels(df_loop1$pathologic_stage)) {
    df_loop2 <- df_loop1[ which(df_loop1$pathologic_stage==Etnia), ]
    freq <- nrow(df_loop2)
    data = rbind(data, data.frame(Etnia,Race,freq, fix.empty.names = F))
  }
  for (Type_sample in levels(df_loop1$sample_type)) {
    df_loop2 <- df_loop1[ which(df_loop1$sample_type==Type_sample), ]
    freq <- nrow(df_loop2)
    data = rbind(data, data.frame(Type_sample,Race,freq, fix.empty.names = F))
  }
}
Total = data.frame(individual=character())
for (n in levels(df_plot3$pathologic_stage)) {
  df_loop3 <- df_plot3[ which(df_plot3$pathologic_stage==n), ]
  Total = rbind(Total, data.frame(n, nrow(df_loop3), fix.empty.names = F))
}
for (n in levels(df_plot3$sample_type)) {
  df_loop3 <- df_plot3[ which(df_plot3$sample_type==n), ]
  Total = rbind(Total, data.frame(n, nrow(df_loop3), fix.empty.names = F))
}
colnames(Total) <- c("Etapa o tipo de muestra", "n")
colnames(data) <- c("Etapa o tipo de muestra", "Diagn�stico", "Proporci�n")
data <- as.data.frame(merge(data, Total, by = "Etapa o tipo de muestra"))

pdf("df_plot4.pdf", height = 8.5, width = 13)
ggplot(data, aes(fill=Diagn�stico, y=Proporci�n, x=`Etapa o tipo de muestra`)) + 
  geom_bar(position="fill", stat="identity")+ 
  ggtitle("Distribuci�n de diagnostico, seg�n Etapa y tipo de muestra") +
  geom_text(data=data, aes(x=`Etapa o tipo de muestra`, y= rep(1.02, length(`Etapa o tipo de muestra`)), label=paste("n =", as.factor(n), sep = " "))) 
dev.off()
bitmap <- pdf_render_page("df_plot4.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, "df_plot4.png")
unlink("df_plot4.pdf")
unlink("df_plot4")



