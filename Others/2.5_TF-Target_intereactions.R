#!/usr/bin/env Rscript
Taruca_adolfo_tesis <- "adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/"
## Cargar librerya de interacciones
library(dorothea) 
## Cargar archivo con elementos que se expresan en tejido mamario normalmente
matrix<- read.delim("../1_expression_data/GTEx_data/2_mRNAs_gene-level_TMM_GTEx_matrix_with_gene_symbol.tab")
reg_info <-  dorothea_hs
tf_in_breast <- reg_info[reg_info$tf %in% rownames(matrix),]
tf_in_breast <- tf_in_breast[tf_in_breast$target %in% rownames(matrix),]
## Obtener informacion de cantidades de interacciones totales y segun grado de confianza
dim(tf_in_breast)
dim(tf_in_breast[tf_in_breast$confidence == "A",])
dim(tf_in_breast[tf_in_breast$confidence == "B",])
dim(tf_in_breast[tf_in_breast$confidence == "C",])
dim(tf_in_breast[tf_in_breast$confidence == "D",])
dim(tf_in_breast[tf_in_breast$confidence == "E",])
dim(reg_info)
## Descartar interraciones de mas baja confianza, que solo son respaldadas por prediccion en base a motivos de union de TF
tf_in_breast <- tf_in_breast[tf_in_breast$confidence != "E",]
dim(tf_in_breast)
write.table(tf_in_breast, sep = "\t", file ="5_Interacciones_TF-Target_tejido_mamario.tsv", row.names = F, quote = F, col.names = T)
system(paste("scp ", "5_Interacciones_TF-Target_tejido_mamario.tsv ",Taruca_adolfo_tesis, "5_Interacciones_TF-Target_tejido_mamario.tsv", sep = ""))
