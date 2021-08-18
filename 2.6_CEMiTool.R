#!/usr/bin/env Rscript
library("CEMiTool")
system("mkdir -p normal_vs_tumoral")
gmt_info <- read_gmt("c2.cp.v7.4.symbols.gmt")
PP_int_df <- read.delim("../GTEx_TFs_PPI/2_interacciones_sobre_700_gene_symbol_in_GTEx.tab", sep = "\t")
expr0 <- read.delim("../analisis_muestras_y_DE/2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_id.tab", sep = "\t")
gencode_v38 <- read.delim("../gencode_genes.v38.annotation.tab", sep = "\t")
gencode_v38 <- gencode_v38[c("gene_id","gene_name")]
gencode_v38 <- gencode_v38[!duplicated(gencode_v38),]
gencode_v38 <- gencode_v38[gencode_v38$gene_id %in% rownames(expr0), ]
gencode_v38_sorted <- gencode_v38[order(gencode_v38$gene_id), ]

for (duplicado in levels(as.factor(gencode_v38_sorted[duplicated(gencode_v38_sorted$gene_name),]$gene_name))){
        count = 0
        for (elementos in rownames(gencode_v38_sorted[gencode_v38_sorted$gene_name == duplicado,])){ 
                if(count != 0) {
                        gencode_v38_sorted[elementos,]$gene_name <- paste(duplicado, count, sep = "_")
                        print(gencode_v38_sorted[elementos, ]$gene_name)}
                count = count + 1}}
write.table(gencode_v38_sorted, file = "2_gene_ID_to_gene_symbol_TCGA.tab", row.names = F, quote = F, col.names = T, sep = "\t")

for (index in rownames(expr0)){
        rownames(expr0)[rownames(expr0) == index] <- gencode_v38_sorted[gencode_v38_sorted$gene_id == index, ]$gene_name
}

PP_int_df <- merge(PP_int_df, gencode_v38_sorted, by.x = "protein1", by.y = "gene_id")
PP_int_df <- PP_int_df[c("gene_name","protein2")]
colnames(PP_int_df)[colnames(PP_int_df) == "gene_name"] <- "protein1"
PP_int_df <- merge(PP_int_df, gencode_v38_sorted, by.x = "protein2", by.y = "gene_id")
PP_int_df <- PP_int_df[c("gene_name","protein1")]
colnames(PP_int_df)[colnames(PP_int_df) == "gene_name"] <- "protein2"

colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
sample_annot <- read.delim("../analisis_muestras_y_DE/2_sample_annot_corregido.tab", sep = "\t")
expr0 <- expr0[sample_annot$SampleName]
sample_annot$Class <- sample_annot$sample_type
sample_annot$Class[sample_annot$Class == "NT"] <- "Normal control"
sample_annot$Class[sample_annot$Class == "TP"] <- "Tumoral"
sample_annot <- sample_annot[c("SampleName","Class")]
if (identical(colnames(expr0), sample_annot$SampleName) == T){
        cem <- cemitool(expr0, sample_annot, gmt_info, interactions=PP_int_df, filter= T, plot = TRUE, verbose=T, apply_vst = T)
        
        write_files(cem, directory="normal_vs_tumoral/Tables", force = T) # write analysis results into files
        save_plots(cem, "all", directory="normal_vs_tumoral/Plots", force = T)# save all plots
} else {
        print("problema de ejecucion")}

system("mkdir -p subtipos_vs_normal")

gmt_info <- read_gmt("c5.go.bp.v7.4.symbols.gmt")
sample_annot <- read.delim("../analisis_muestras_y_DE/2_sample_annot_corregido.tab", sep = "\t")
sample_annot$Class <- sample_annot$BRCA_Subtype_PAM50
sample_annot[sample_annot$sample_type == "NT",]$Class <- "Normal control"
sample_annot <- sample_annot[c("SampleName","Class")]

if (identical(colnames(expr0), sample_annot$SampleName) == T){
        cem <- cemitool(expr0, sample_annot, gmt_info, interactions=PP_int_df, filter= T, plot=TRUE, verbose=T, apply_vst = T)
        write_files(cem, directory="subtipos_vs_normal/Tables", force = T)# write analysis results into files
        save_plots(cem, "all", directory="subtipos_vs_normal/Plots", force = T)# save all plots
        print(cem)
} else {
        print("problema de ejecucion")}
system("scp -r subtipos_vs_normal adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/")
system("scp -r normal_vs_tumoral adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/")