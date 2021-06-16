#!/usr/bin/env Rscript
library("CEMiTool")
system("mkdir -p normal_vs_tumoral")
gmt_info <- read_gmt("c5.all.v7.2.symbols.gmt")
PP_int_df <- read.delim("../GTEx_TFs_PPI/2_interacciones_sobre_700_gene_symbol_in_GTEx.tab", sep = "\t")
expr0 <- read.delim("../analisis_muestras_y_DE/2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_symbol.tab", sep = "\t")
colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
sample_annot <- read.delim("../analisis_muestras_y_DE/2_sample_annot_corregido.tab", sep = "\t")
expr0 <- expr0[sample_annot$SampleName]
sample_annot$Class <- sample_annot$sample_type
sample_annot$Class[sample_annot$Class == "NT"] <- "Normal control"
sample_annot$Class[sample_annot$Class == "TP"] <- "Tumoral"
sample_annot <- sample_annot[c("SampleName","Class")]


if (identical(colnames(expr0), sample_annot$SampleName) == T){

        cem <- cemitool(expr0, sample_annot, gmt_info, interactions=PP_int_df, filter= T, plot=TRUE, verbose=T, apply_vst = T)


        # write analysis results into files
        write_files(cem, directory="normal_vs_tumoral/Tables", force = T)

        # save all plots
        save_plots(cem, "all", directory="normal_vs_tumoral/Plots", force = T)

} else {
        print("problema de ejecucion")}

system("mkdir -p subtipos_vs_normal")

#gmt_info <- read_gmt("c5.all.v7.2.symbols.gmt")
#expr0 <- read.delim("2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_symbol.tab", sep = "\t")
#colnames(expr0) <- gsub("\\.", "-", colnames(expr0))
sample_annot <- read.delim("../analisis_muestras_y_DE/2_sample_annot_corregido.tab", sep = "\t")
sample_annot$Class <- sample_annot$BRCA_Subtype_PAM50
sample_annot[sample_annot$sample_type == "NT",]$Class <- "Normal control"
sample_annot <- sample_annot[c("SampleName","Class")]

if (identical(colnames(expr0), sample_annot$SampleName) == T){

        cem <- cemitool(expr0, sample_annot, gmt_info, interactions=PP_int_df, filter= T, plot=TRUE, verbose=T, apply_vst = T)


        # write analysis results into files
        write_files(cem, directory="subtipos_vs_normal/Tables", force = T)

        # save all plots
        save_plots(cem, "all", directory="subtipos_vs_normal/Plots", force = T)
        print(cem)

} else {
        print("problema de ejecucion")}
system("scp -r subtipos_vs_normal adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/")
system("scp -r normal_vs_tumoral adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/")