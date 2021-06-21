#!/usr/bin/env Rscript
## PPI
protein_int <- read.delim("GTEx_TFs_PPI/2_interacciones_sobre_700_gene_symbol_in_GTEx.tab")
names(protein_int) <- c("element1","element2")
protein_int <- protein_int[!duplicated(protein_int),]
protein_int$int_type <- "PPI"
## TFs
tfs_int <- read.delim("GTEx_TFs_PPI/5_Interacciones_TF-Target_tejido_mamario.tsv")
tfs_int <- tfs_int[c("tf","target")]
names(tfs_int) <- c("element1","element2")
tfs_int <- tfs_int[!duplicated(tfs_int),]
tfs_int2 <- tfs_int
names(tfs_int2) <- c("element2","element1")
tfs_int <- rbind(tfs_int,tfs_int2)
tfs_int <- tfs_int[!duplicated(tfs_int),]
tfs_int$int_type <- "TF-Target"
## Co-Expression
coex_modules <- read.delim("co-expression/subtipos_vs_normal/Tables/module.tsv")
coex_int <- data.frame()
for (module in levels(as.factor(coex_modules$modules))) {
   print(module)
   module_genes <- coex_modules[coex_modules$modules == module,]$genes
   for (gene in module_genes) {
      df_loop <- data.frame(element1 = rep(gene,length(module_genes)-1), element2 = module_genes[module_genes != gene], int_type = rep(module, length(module_genes)-1))
      coex_int <- rbind(coex_int,df_loop)}}
coex_int <- coex_int[!duplicated(coex_int),]
## ceRNA
load("Sponge/1_entorno_completo_sponge2.RData")
ceRNA_int <- ceRNA_interactions_fdr[c("geneA","geneB")]
names(ceRNA_int) <- c("element2","element1")
ceRNA_int <- ceRNA_int[!duplicated(ceRNA_int),]
ceRNA_int2 <- ceRNA_int
names(ceRNA_int2) <- c("element1","element2")
ceRNA_int <- rbind(ceRNA_int2,ceRNA_int)

ceRNA_int$int_type <- "ceRNA"
ceRNA_int <- ceRNA_int[!duplicated(ceRNA_int),]


head(genes_miRNA_candidates)
for (interactor in names(genes_miRNA_candidates)){
   if (length(genes_miRNA_candidates[[interactor]][,1]) != 0){
      for(i in 1:length(genes_miRNA_candidates[[interactor]][,1])){
         genes_miRNA_candidates2 <- rbind(genes_miRNA_candidates2,as.data.frame(cbind(interactor, genes_miRNA_candidates[[interactor]][i,])))}}}
genes_miRNA_candidates2$coefficient <- NULL
names(genes_miRNA_candidates2) <- c("element2","element1")
genes_miRNA_candidates3 <- genes_miRNA_candidates2
names(genes_miRNA_candidates3) <- c("element1","element2")
genes_miRNA_candidates2 <- rbind(genes_miRNA_candidates2,genes_miRNA_candidates3)
genes_miRNA_candidates2 <- genes_miRNA_candidates2[!duplicated(genes_miRNA_candidates2),]
genes_miRNA_candidates2$int_type <- "miRNA-mRNA"

genes_miRNA_candidates2 <- genes_miRNA_candidates2[!duplicated(genes_miRNA_candidates2),]

final_int <- rbind(ceRNA_int,coex_int, tfs_int, protein_int, genes_miRNA_candidates2)
write.table(final_int, sep = "\t", file = "final_interaction.tab", row.names = F, quote = F, col.names = T)

protein_int[protein_int$element1 == "KIF13B" & protein_int$element2 == "ARF5",]