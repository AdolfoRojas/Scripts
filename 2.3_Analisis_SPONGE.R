#!/usr/bin/env Rscript
library(edgeR)
library(tidyverse)
library(pdftools)
library(SPONGE)
N_CORES <- 60

interaction_matrix <- read.delim("matrix_anotada_full.txt", h =  T)
RNA_exp <- read.delim("../analisis_muestras_y_DE/1_mRNAs_gene-level_Counts_TCGA_matrix.tab", h = T, row.names = 1, stringsAsFactors = F)
colnames(RNA_exp) <- gsub("\\.", "-", colnames(RNA_exp))
colnames(RNA_exp) <- substring(colnames(RNA_exp), 1, 15)
table(interaction_matrix$Gene.stable.ID%in%rownames(RNA_exp))
interaction_matrix <- interaction_matrix[interaction_matrix$Gene.stable.ID%in%rownames(RNA_exp),]
RNA_exp <- RNA_exp[rownames(RNA_exp)%in%interaction_matrix$Gene.stable.ID,]
table(Biobase::isUnique(interaction_matrix$Gene.stable.ID))
interaction_matrix2 <- aggregate(interaction_matrix[,4:ncol(interaction_matrix)],list(interaction_matrix$Gene.stable.ID), sum)
rownames(interaction_matrix2) <- interaction_matrix2$Group.1
interaction_matrix2$Group.1 <- NULL
RNA_exp <- RNA_exp[as.character(rownames(interaction_matrix2)),]
table(rownames(RNA_exp)==rownames(interaction_matrix2))

miRNA <- read.delim("matrix_miRNAs_isoforms_with_names.tab", sep = "\t", h = T,stringsAsFactors = F)
rownames(miRNA) <- miRNA$TargetName
miRNA$TargetName <- NULL
colnames(miRNA) <- gsub("\\.", "-", colnames(miRNA))
colnames(miRNA) <- substring(colnames(miRNA), 1, 15)
RNA_exp <- RNA_exp[,colnames(RNA_exp)%in%colnames(miRNA)]
miRNA <- miRNA[,colnames(miRNA)%in%colnames(RNA_exp)]
miRNA <- miRNA[,as.character(colnames(RNA_exp))]
table(colnames(miRNA)==colnames(RNA_exp))
colnames(interaction_matrix2) <- gsub("\\.", "-", colnames(interaction_matrix2))
interaction_matrix2 <- interaction_matrix2[, colnames(interaction_matrix2)%in%rownames(miRNA)]
miRNA <- miRNA[rownames(miRNA)%in%colnames(interaction_matrix2),]
miRNA <- miRNA[as.character(colnames(interaction_matrix2)),]

metadata <- read.delim("../analisis_muestras_y_DE/2_sample_annot_corregido.tab", h = T, stringsAsFactors = F)
metadata$SampleName <- substring(metadata$SampleName, 1, 15)
metadata <- metadata[!duplicated(metadata$SampleName),]
rownames(metadata) <- metadata$SampleName
metadata <- metadata[metadata$SampleName%in%colnames(miRNA),]
table(metadata$SampleName%in%colnames(miRNA))
miRNA <- miRNA[, colnames(miRNA)%in%metadata$SampleName]
metadata <- metadata[as.character(colnames(miRNA)),]
table(rownames(metadata)==colnames(miRNA))
RNA_exp <- RNA_exp[colnames(miRNA)]
table(rownames(metadata)==colnames(RNA_exp))

metadata$group <- metadata$BRCA_Subtype_PAM50
metadata$group[metadata$sample_type=="NT"] <- "NT"

raw_miR <- DGEList(counts = miRNA, samples = metadata, genes = rownames(miRNA), remove.zeros = T, group = metadata$group)
cpm_miR <- cpm(raw_miR$counts)
lcpm_miR <- cpm(raw_miR, log=TRUE)
mean(raw_miR$samples$lib.size) * 1e-6
median(raw_miR$samples$lib.size) * 1e-6
table(rowSums(cpm(raw_miR)>=4) >=50)
keep <- rowSums(cpm(raw_miR)>=4) >= 50  ### equivale a 10
filtered_mir <- raw_miR[keep, keep.lib.sizes=FALSE]
dim(filtered_mir)
norm_mir <- calcNormFactors(filtered_mir, method = "TMM")

raw_RNA <- DGEList(counts = RNA_exp, samples = metadata, genes = rownames(RNA_exp), remove.zeros = T, group = metadata$group)
cpm_RNA <- cpm(raw_RNA$counts)
lcpm_RNA <- cpm(raw_RNA, log=TRUE)
mean(raw_RNA$samples$lib.size) * 1e-6
median(raw_RNA$samples$lib.size) * 1e-6
table(rowSums(cpm(raw_RNA)>=0.263) >=50)  ### para mRNA 15 
keep <- rowSums(cpm(raw_RNA)>=0.263) >= 50
filtered_RNA <- raw_RNA[keep, keep.lib.sizes=FALSE]
dim(filtered_RNA)
norm_RNA <- calcNormFactors(filtered_RNA, method = "TMM")

interaction_matrix3 <- interaction_matrix2[, colnames(interaction_matrix2)%in%rownames(norm_mir)]
norm_mir <- norm_mir[rownames(norm_mir)%in%colnames(interaction_matrix3),]
norm_mir <- norm_mir[as.character(colnames(interaction_matrix3)),]
interaction_matrix3 <- interaction_matrix3[rownames(interaction_matrix3)%in%rownames(norm_RNA),]
norm_RNA <- norm_RNA[rownames(norm_RNA)%in%rownames(interaction_matrix3),]
norm_RNA <- norm_RNA[as.character(rownames(interaction_matrix3)),]
table(rownames(norm_RNA)==rownames(interaction_matrix3))
table(rownames(norm_mir)==colnames(interaction_matrix3))

head(colnames(norm_mir))
print("convertir a log")
norm_mir$counts <- cpm(norm_mir$counts, log=TRUE)
norm_RNA$counts <- cpm(norm_RNA$counts, log=TRUE)
print("trasponer")
miR_exp <- t(norm_mir$counts)
rna_exp <- t(norm_RNA$counts)
interaction_matrix3 <- data.matrix(interaction_matrix3)
library(doParallel)
library(foreach)
num.of.cores <- N_CORES
logging.file <- "./BRCA_interactions_v1.log"
cl <- makeCluster(num.of.cores, outfile=logging.file) 
registerDoParallel(cl)

genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(gene_expr = rna_exp, mir_expr = miR_exp, mir_predicted_targets = interaction_matrix3)
ceRNA_interactions <- sponge(gene_expr = rna_exp, mir_expr = miR_exp, mir_interactions = genes_miRNA_candidates)
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100000, number_of_samples = nrow(rna_exp))

pdf("df_plot1.pdf", height = 16, width = 9)
sponge_plot_simulation_results(mscor_null_model)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 500)
png::writePNG(bitmap, "df_plot1.png")

ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, null_model = mscor_null_model)
stopCluster(cl) 
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < 0.01),]

save.image(file = "1_entorno_completo_sponge2.RData")

###################################################
#              Analizar
###################################################

library(limma)
#load("entorno_completo_sponge2.RData")

pdf("df_plot1.pdf", height = 9, width = 16)
boxplot(norm_RNA$counts[,1:100], las = 2)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 500)
png::writePNG(bitmap, "norm_rna_boxplot.png")
pdf("df_plot1.pdf", height = 9, width = 16)
boxplot(norm_mir$counts[,1:100], las = 2)
dev.off()
bitmap <- pdf_render_page("df_plot1.pdf", page = 1, dpi = 500)
png::writePNG(bitmap, "norm_mir_boxplot.png")

gencode <- read.delim("gencodeV35_fixed.txt")
ids <- as.data.frame(unique(append(ceRNA_interactions_fdr$geneA,ceRNA_interactions_fdr$geneB)))
colnames(ids) <- "gene_id"
gencode <- merge(gencode, ids, by = "gene_id")
levels(as.factor(gencode$gene_type))

gene_type_protein_coding <- gencode[gencode$gene_type=="protein_coding",]
levels(as.factor(gene_type_protein_coding$transcript_type))

gencode <- unique(gencode[c("gene_id", "gene_type", "gene_name")])
gene_type_lncRNA <- gencode[gencode$gene_type=="lncRNA",]
gene_type_lncRNA <- gene_type_lncRNA["gene_id"]


########## Inicio Loop de subtipos ##########
subtipos <- c("all", "LumA", "LumB", "Her2", "Basal")
for (subtype in subtipos){
        if (subtype == "all"){
                DE_info <- read.delim("../analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_mRNAs-gene_level/all_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",")
        } else{
        DE_info <- read.delim(paste("../analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_mRNAs-gene_level/", subtype,"_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ""), sep = ",")
        }
        DE_info <- DE_info[DE_info$padj < 0.05,]
        DE_info <- DE_info[DE_info$log2FoldChange < -1 | DE_info$log2FoldChange > 1,]

        colnames(DE_info)[1] <- "gene_id"

        DE_lncRNA <- merge(gene_type_lncRNA,DE_info, by = "gene_id")
        DE_lncRNA_ids <- as.data.frame(DE_lncRNA$gene_id)
        DE_lncRNA_up <- DE_lncRNA[DE_lncRNA$log2FoldChange < 1,]
        DE_lncRNA_down <- DE_lncRNA[DE_lncRNA$log2FoldChange > -1,]
        colnames(DE_lncRNA_ids) <- "gene_id"
        
        new_table <- merge(gencode, ceRNA_interactions_fdr, by.x = "gene_id", by.y = "geneB")
        
        colnames(new_table)[1:3] <- c("GeneB", "gene_typeB", "gene_nameB")
        new_table <- merge(gencode, new_table, by.x = "gene_id", by.y = "geneA")
        colnames(new_table)[1:3] <- c("GeneA", "gene_typeA", "gene_nameA")
        
        new_table2 <- merge(gencode, ceRNA_interactions_fdr, by.x = "gene_id", by.y = "geneA")

        colnames(new_table2)[1:3] <- c("GeneA", "gene_typeA", "gene_nameA")
        new_table2 <- merge(gencode, new_table2, by.x = "gene_id", by.y = "geneB")
        colnames(new_table2)[1:3] <- c("GeneB", "gene_typeB", "gene_nameB")     
        new_table <- unique(rbind(new_table, new_table2))
        head(new_table)

        lncRNA_ceRNA_interactions <- new_table[grepl("lncRNA", new_table$gene_typeA) | grepl("lncRNA", new_table$gene_typeB),]
        A <- merge(DE_lncRNA_ids, lncRNA_ceRNA_interactions, by.x = "gene_id", by.y = "GeneA")
        colnames(A)[1] <- "GeneA"
        B <- merge(DE_lncRNA_ids, lncRNA_ceRNA_interactions, by.x = "gene_id", by.y = "GeneB")

        names(B)[names(B) =="GeneA"] <- "GeneC"
        names(B)[names(B) =="gene_typeA"] <- "gene_typeC"
        names(B)[names(B) =="gene_nameA"] <- "gene_nameC"
        colnames(B)[1] <- "GeneA"
        names(B)[names(B) =="gene_typeB"] <- "gene_typeA"
        names(B)[names(B) =="gene_nameB"] <- "gene_nameA"
        names(B)[names(B) =="GeneC"] <- "GeneB"
        names(B)[names(B) =="gene_typeC"] <- "gene_typeB"
        names(B)[names(B) =="gene_nameC"] <- "gene_nameB"
        B <- B[B$gene_typeB != "lncRNA",]

        DE_lncRNA_ceRNA_interactions <- rbind(A,B) #####

        save(DE_lncRNA_ceRNA_interactions, file = paste("4_", subtype, "_DE_lncRNA_ceRNA_interactions.Rdata", sep = ""))

        DE_lncRNA_ceRNA_interactions <- DE_lncRNA_ceRNA_interactions[DE_lncRNA_ceRNA_interactions$gene_typeB == "protein_coding"| DE_lncRNA_ceRNA_interactions$gene_typeB == "lncRNA",]

        write.table(DE_lncRNA_ceRNA_interactions, file = paste("4_DE_", subtype, "_lncRNA_ceRNA_interactions.tab", sep = ""), row.names = F, quote = F, col.names = T, sep = "\t")
        save(DE_lncRNA, file = paste("4_DE_", subtype, "_lncRNA.Rdata", sep = ""))
        write.table(DE_lncRNA, file = paste("4_DE_lncRNA_", subtype, ".tab", sep = ""), row.names = F, quote = F, col.names = T, sep = "\t")
        save(genes_miRNA_candidates, file = paste("4_", subtype, "_genes_miRNA_candidates.Rdata", sep = ""))
        genes_miRNA_candidates2 <- NULL
        for (interactor in names(genes_miRNA_candidates)){
                if (length(genes_miRNA_candidates[[interactor]][,1]) != 0){
                        for(i in 1:length(genes_miRNA_candidates[[interactor]][,1])){
                                genes_miRNA_candidates2 <- rbind(genes_miRNA_candidates2,as.data.frame(cbind(interactor, genes_miRNA_candidates[[interactor]][i,])))
                        }               
                }
        }

        genes_miRNA_candidates3 <- genes_miRNA_candidates2[genes_miRNA_candidates2$interactor %in% DE_lncRNA_ceRNA_interactions$GeneA,]
        genes_miRNA_candidates4 <- genes_miRNA_candidates2[genes_miRNA_candidates2$interactor %in% DE_lncRNA_ceRNA_interactions$GeneB,]
        genes_miRNA_candidates4 <- genes_miRNA_candidates4[genes_miRNA_candidates4$mirna %in% genes_miRNA_candidates3$mirna,]
        genes_miRNA_candidates3 <- genes_miRNA_candidates3[genes_miRNA_candidates3$mirna %in% genes_miRNA_candidates4$mirna,]
        genes_miRNA_candidates2 <- rbind(genes_miRNA_candidates3,genes_miRNA_candidates4)
        genes_miRNA_candidates2 <- unique(genes_miRNA_candidates2)

        genes_miRNA_candidates2 <- merge(gencode, genes_miRNA_candidates2, by.x = "gene_id", by.y = "interactor")
        print(dim(genes_miRNA_candidates2))

        names(genes_miRNA_candidates2)[names(genes_miRNA_candidates2) =="gene_id"] <- "GeneB"
        names(genes_miRNA_candidates2)[names(genes_miRNA_candidates2) =="gene_type"] <- "gene_typeB"
        names(genes_miRNA_candidates2)[names(genes_miRNA_candidates2) =="gene_name"] <- "gene_nameB"
        names(genes_miRNA_candidates2)[names(genes_miRNA_candidates2) =="mirna"] <- "gene_nameA"
        genes_miRNA_candidates2 <- genes_miRNA_candidates2[genes_miRNA_candidates2$gene_typeB == "lncRNA" | genes_miRNA_candidates2$gene_typeB == "protein_coding",]
        genes_miRNA_candidates2 <- genes_miRNA_candidates2[genes_miRNA_candidates2$gene_nameB %in% DE_lncRNA_ceRNA_interactions$gene_nameB | genes_miRNA_candidates2$gene_nameB %in% DE_lncRNA_ceRNA_interactions$gene_nameA,]
        write.table(genes_miRNA_candidates2, file = paste("4_genes_miRNA_candidates_", subtype,".tab", sep = ""), row.names = F, quote = F, col.names = T, sep = "\t")

        DE_lncRNA_ceRNA_interactions <- DE_lncRNA_ceRNA_interactions[c("GeneA","gene_nameA","gene_typeA","GeneB","gene_typeB","gene_nameB")]
        genes_miRNA_candidates2$gene_typeA <- "miRNA"
        genes_miRNA_candidates2$GeneA <- NA
        genes_miRNA_candidates2$coefficient <- NULL

        cytoscape <- rbind(DE_lncRNA_ceRNA_interactions,genes_miRNA_candidates2)

        write.table(cytoscape, file = paste("4_to_cytoscape_", subtype,".tab", sep = ""), row.names = F, quote = F, col.names = T, sep = "\t")
        mirRNAs_unicos <- sapply(genes_miRNA_candidates, function(x){x[1]})
        mirRNAs_unicos <- unlist(mirRNAs_unicos, use.names = F)
        mirRNAs_unicos <- unique(mirRNAs_unicos)
        length(mirRNAs_unicos)}