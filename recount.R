library(recount)
load(file.path('SRP012682', 'rse_tx_breast.Rdata'))  TCGA_breast.RData
rse <- scale_counts(rse_tx)
GTEx <- rse_tx[,grepl("female", rse_tx$bigwig_path) == TRUE]

##TCGA
assays(se)$counts ##########################
TCGA <- rse_tx[,grepl("female", rse$gdc_cases.demographic.gender) == TRUE]
library(DESeq2)
TCGA$gdc_cases.samples.sample_type <- gsub("_type != "Metastatic"]
 ", "_", TCGA$gdc_cases.samples.sample_type)
TCGA <- TCGA[,TCGA$gdc_cases.samples.sample


colData(rse)$gdc_cases.demographic.gender or $xml_gender
colData(rse)$cgc_sample_sample_type



colnames(colData(rse))[grep("molecular", colnames(colData(rse)))]