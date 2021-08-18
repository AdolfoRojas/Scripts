#!/usr/bin/env Rscript

genesymbol_to_geneid <- read.delim("gencode_genes.v38.annotation.tab")
genesymbol_to_geneid <- genesymbol_to_geneid[c("gene_name","gene_id")]

## PPI
protein_int <- read.delim("GTEx_TFs_PPI/2_interacciones_sobre_700_gene_symbol_in_GTEx.tab")
names(protein_int) <- c("element1","element2")
protein_int <- protein_int[!duplicated(protein_int),]
protein_int$int_type <- "PPI"

## TFs
tfs_int <- read.delim("GTEx_TFs_PPI/5_Interacciones_TF-Target_tejido_mamario.tsv")
tfs_int <- tfs_int[c("tf","target")]
names(tfs_int) <- c("element1","element2")
tfs_int <- tfs_int[!duplicated(tfs_int),] # 21666     2
alias <- tfs_int[!(tfs_int$element2 %in% genesymbol_to_geneid$gene_name),]
tfs_int <- merge(tfs_int, genesymbol_to_geneid, by.x = "element1", by.y = "gene_name")
tfs_int <- tfs_int[!duplicated(tfs_int),] #
head(tfs_int)
tfs_int$element1 <- NULL 
head(tfs_int)
names(tfs_int)[names(tfs_int) == "gene_id"] <- "element1"
tfs_int <- merge(tfs_int, genesymbol_to_geneid, by.x = "element2", by.y = "gene_name")
tfs_int <- tfs_int[!duplicated(tfs_int),] #
head(tfs_int)
tfs_int$element2 <- NULL
head(tfs_int) 
names(tfs_int)[names(tfs_int) == "gene_id"] <- "element2"
tfs_int <- tfs_int[!duplicated(tfs_int),] #
head(tfs_int) 

library(httr)
library(jsonlite)
library(xml2)

alias$gene_id <- "desconocido"
alias$n_id_encontrados <- 0
server <- "https://rest.ensembl.org"

for (identifier in unique(alias$element2)) {
   
   ext <- paste("/xrefs/symbol/homo_sapiens/", identifier, "?", sep = "")
   r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
   stop_for_status(r)    
   alias[alias$element2== identifier,]$n_id_encontrados <- length(fromJSON(toJSON(content(r)))$id)
   print(alias[alias$element2== identifier,]$n_id_encontrados)
   # use this if you get a simple nested list back, otherwise inspect its structure
   # head(data.frame(t(sapply(content(r),c))))
   ens_gene_id <- as.character(fromJSON(toJSON(content(r)))$id[[1]])
   if (length(ens_gene_id) > 0) {             
      alias[alias$element2== identifier,]$gene_id <- ens_gene_id}}

length(unique(alias[alias$n_id_encontrados== 0,]$element2))
alias <- alias[alias$n_id_encontrados!= 0,]
length(unique(alias[alias$n_id_encontrados == 2,]$element2))
length(unique(alias[alias$n_id_encontrados > 2,]$element2))
alias$element2 <- NULL
names(alias)[names(alias) == "gene_id"] <- "element2"
alias$n_id_encontrados <- NULL
tfs_int <- rbind(tfs_int,alias)
tfs_int <- tfs_int[!duplicated(tfs_int),]
tfs_int$int_type <- "TF-Target"

## Co-Expression
#coex_modules <- read.delim("co-expression/subtipos_vs_normal/Tables/module.tsv")
#coex_int <- data.frame()
#for (module in levels(as.factor(coex_modules$modules))) {
#   print(module)
#   module_genes <- coex_modules[coex_modules$modules == module,]$genes
#   for (gene in module_genes) {
#      df_loop <- data.frame(element1 = rep(gene,length(module_genes)-1), element2 = module_genes[module_genes != gene], int_type = rep(module, length(module_genes)-1))
#      coex_int <- rbind(coex_int,df_loop)}}
#coex_int <- coex_int[!duplicated(coex_int),]

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

genes_miRNA_candidates2 <- NULL
head(genes_miRNA_candidates)
for (interactor in names(genes_miRNA_candidates)){
   if (length(genes_miRNA_candidates[[interactor]][,1]) != 0){
      for(i in 1:length(genes_miRNA_candidates[[interactor]][,1])){
         genes_miRNA_candidates2 <- rbind(genes_miRNA_candidates2,as.data.frame(cbind(interactor, genes_miRNA_candidates[[interactor]][i,])))}}}
genes_miRNA_candidates2$coefficient <- NULL
names(genes_miRNA_candidates2) <- c("element2","element1")
genes_miRNA_candidates2 <- genes_miRNA_candidates2[!duplicated(genes_miRNA_candidates2),]
genes_miRNA_candidates2$int_type <- "miRNA-mRNA"
genes_miRNA_candidates2 <- genes_miRNA_candidates2[!duplicated(genes_miRNA_candidates2),]

final_int <- rbind(ceRNA_int, tfs_int, protein_int, genes_miRNA_candidates2)
write.table(final_int, sep = "\t", file = "final_interaction.tab", row.names = F, quote = F, col.names = T)