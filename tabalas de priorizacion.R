#!/usr/bin/env Rscript
genesymbol_to_geneid <- read.delim("gencode_genes.v38.annotation.tab")
genesymbol_to_geneid <- genesymbol_to_geneid[c("gene_id", "gene_type", "gene_name")]

DE_data_mRNA <- read.delim("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_mRNAs-gene_level/all_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
dim(DE_data_mRNA)
DE_data_mRNA <- merge(DE_data_mRNA, genesymbol_to_geneid, by.x = "X", by.y = "gene_id")
DE_data_mRNA <- DE_data_mRNA[!duplicated(DE_data_mRNA$X),]
dim(DE_data_mRNA)
DE_data_mRNA$S_DE_all <- 0
DE_data_mRNA[DE_data_mRNA$gene_type == "lncRNA" & (DE_data_mRNA$log2FoldChange >= 1 | DE_data_mRNA$log2FoldChange <= -1),]$S_DE_all <- 1
DE_data_mRNA[DE_data_mRNA$gene_type != "lncRNA" & (DE_data_mRNA$log2FoldChange >= 2 | DE_data_mRNA$log2FoldChange <= -2),]$S_DE_all <- 1
DE_data_mRNA_final <- DE_data_mRNA[c("X","gene_type","S_DE_all")]

subtipos <- c("Basal", "Her2", "LumB", "LumA")

for (i in subtipos) {
    print(i)
    DE_data_mRNA <- read.delim(paste("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", i, "_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
    print(dim(DE_data_mRNA))
    DE_data_mRNA <- merge(DE_data_mRNA, genesymbol_to_geneid, by.x = "X", by.y = "gene_id")
    DE_data_mRNA <- DE_data_mRNA[!duplicated(DE_data_mRNA$X),]
    print(dim(DE_data_mRNA))
    DE_data_mRNA$S_DE <- 0
    DE_data_mRNA[DE_data_mRNA$gene_type == "lncRNA" & (DE_data_mRNA$log2FoldChange >= 1 | DE_data_mRNA$log2FoldChange <= -1),]$S_DE <- 1
    DE_data_mRNA[DE_data_mRNA$gene_type != "lncRNA" & (DE_data_mRNA$log2FoldChange >= 2 | DE_data_mRNA$log2FoldChange <= -2),]$S_DE <- 1
    DE_data_mRNA <- DE_data_mRNA[c("X","gene_type","S_DE")]
    names(DE_data_mRNA)[names(DE_data_mRNA) == "S_DE"] <- paste("S_DE_", i, sep = "")
    DE_data_mRNA_final <- merge(DE_data_mRNA_final,DE_data_mRNA, by = intersect(names(DE_data_mRNA_final), names(DE_data_mRNA)), all = T)}
DE_data_mRNA_final[is.na(DE_data_mRNA_final)] <- 0
DE_data_mRNA_final[DE_data_mRNA_final$S_DE_Basal == 1,]$S_DE_Basal <- 2
DE_data_mRNA_final$DE_Score <- rowSums(DE_data_mRNA_final[,3:7])


DE_data_miRNA <- read.delim("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_miRNAs-gene_level/all_samples_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
dim(DE_data_miRNA)
DE_data_miRNA <- DE_data_miRNA[!duplicated(DE_data_miRNA$X),]
dim(DE_data_miRNA)
DE_data_miRNA$S_DE_all <- 0
DE_data_miRNA[DE_data_miRNA$log2FoldChange >= 1 | DE_data_miRNA$log2FoldChange <= -1,]$S_DE_all <- 1
DE_data_miRNA_final <- DE_data_miRNA[c("X","S_DE_all")]
for (i in subtipos) {
    print(i)
    DE_data_miRNA <- read.delim(paste("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_miRNAs-gene_level/", i, "_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
    print(dim(DE_data_miRNA))    
    DE_data_miRNA <- DE_data_miRNA[!duplicated(DE_data_miRNA$X),]
    print(dim(DE_data_miRNA))
    DE_data_miRNA$S_DE <- 0
    DE_data_miRNA[DE_data_miRNA$log2FoldChange >= 1 | DE_data_miRNA$log2FoldChange <= -1,]$S_DE <- 1    
    DE_data_miRNA <- DE_data_miRNA[c("X","S_DE")]
    names(DE_data_miRNA)[names(DE_data_miRNA) == "S_DE"] <- paste("S_DE_", i, sep = "")
    DE_data_miRNA_final <- merge(DE_data_miRNA_final,DE_data_miRNA, by = intersect(names(DE_data_miRNA_final), names(DE_data_miRNA)), all = T)}
DE_data_miRNA_final[is.na(DE_data_miRNA_final)] <- 0
DE_data_miRNA_final[DE_data_miRNA_final$S_DE_Basal == 1,]$S_DE_Basal <- 2
DE_data_miRNA_final$DE_Score <- rowSums(DE_data_miRNA_final[,2:6])
DE_data_miRNA_final$gene_type <- "miRNA"

final_table <- rbind(DE_data_mRNA_final, DE_data_miRNA_final)

###########################################################################

modulos_info <- read.delim("co-expression/normal_vs_tumoral/Tables/module.tsv", sep = "\t")
names(modulos_info) <- c("X","Module")
change_id <- read.delim("co-expression/2_gene_ID_to_gene_symbol_TCGA.tab", sep = "\t")
modulos_info <- merge(modulos_info, change_id, by.x= "X", by.y= "gene_name")
modulos_info <- modulos_info[c("gene_id","Module")]
colnames(modulos_info)[colnames(modulos_info) == "gene_id"] <- "X"
final_table2 <- merge(final_table,modulos_info, by = intersect(names(final_table), names(modulos_info)), all.x = T)
final_table2[is.na(final_table2)] <- "-"
matrix <- read.delim("GTEx_TFs_PPI/2_mRNAs_gene-level_TMM_GTEx_matrix_with_geneID.tab")
final_table3 <- final_table[final_table$X %in% rownames(matrix),]
expressed_not_in_final_table <- rownames(matrix)[!(rownames(matrix) %in% final_table$X)]
expressed_not_in_final_table <-genesymbol_to_geneid[genesymbol_to_geneid$gene_id %in% expressed_not_in_final_table,]