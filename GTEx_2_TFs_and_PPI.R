#!/usr/bin/env Rscript
library(edgeR)
library(ggpubr)
library(EnsDb.Hsapiens.v86)

Taruca_adolfo_tesis <- "adolfo@200.89.65.156:/media/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/"
expr0 <- read.delim("1_Female_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads_without_PAR_Y.gct.gz_gene_level.tab", row.names = 1)
rownames(expr0) <- sapply(strsplit(as.character(rownames(expr0)),'\\.'), "[", 1)
expr0$Description <- NULL

raw_RNA <- DGEList(counts = expr0, samples = NULL, genes = rownames(expr0), remove.zeros = F, group = NULL) ### remove zero originalmente en T
cpm_RNA <- cpm(raw_RNA$counts)
lcpm_RNA <- cpm(raw_RNA, log=TRUE)
mean_lcpm <- rowMeans(lcpm_RNA)

mean_log_cpm <- mean_lcpm #aveLogCPM(raw_RNA$counts)
df_mean_log_cpm <- as.data.frame(mean_log_cpm)
filter_threshold <- log2(0.035) # 0.03159702
p1 <- ggplot() + aes(x=mean_log_cpm) +
    geom_histogram(binwidth=0.2) +
    geom_vline(xintercept=filter_threshold, colour="red") +
    ggtitle("Distribution of mean gene expression values") + theme_minimal()
p2 <- ggplot(df_mean_log_cpm, aes(sample = mean_log_cpm)) + stat_qq() + geom_hline(yintercept=filter_threshold, colour="red") + theme_minimal() + ggtitle("Progression of mean gene expression values")
ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("Expression_threshold.png", width = 16, height = 9, dpi = 600, units = "in")
#system("scp -P 1313 Expression_threshold.png adolfo@200.89.65.156:/media/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/")

print(mean(raw_RNA$samples$lib.size) * 1e-6)
print(median(raw_RNA$samples$lib.size) * 1e-6)
table(rowSums(cpm(raw_RNA)>=0.035) >= 0.1*ncol(expr0))  ### para mRNA 15
keep <- rowSums(cpm(raw_RNA)>=0.035) >= 0.1*ncol(expr0) #0.263
filtered_RNA <- raw_RNA[keep, keep.lib.sizes=FALSE]
dim(filtered_RNA)
norm_RNA <- calcNormFactors(filtered_RNA, method = "TMM")
expr0 <- norm_RNA$counts
expr0 <- cpm(expr0, log= F)
#expr0 <- expr0[rownames(expr0) %in% geneIDs1$GENEID,] 

#for (ensembl_gene_id in rownames(expr0)){
#        rownames(expr0)[rownames(expr0)== ensembl_gene_id] <- geneIDs1[geneIDs1$GENEID == ensembl_gene_id,]$SYMBOL}

#print("guardando tablas y procediendo con anotacion en gene symbol de interacciones proteina-proteina")
write.table(expr0, file = "2_mRNAs_gene-level_TMM_GTEx_matrix_with_geneID.tab", row.names = T, quote = F, col.names = T, sep = "\t") ###############################################################################

p <- read.delim("interactions_above_400_combined_score.txt", sep=" ")
p <- p[p$combined_score >= 400,]
p$protein1 <- sapply(strsplit(as.character(p$protein1),'\\.'), "[", 2)
p$protein2 <- sapply(strsplit(as.character(p$protein2),'\\.'), "[", 2)
ensembl.genes <- rownames(expr0)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID","PROTEINID"))
as <- as.data.frame(levels(as.factor(p$protein1)))
az <- merge(geneIDs1, as, by.y ="levels(as.factor(p$protein1))", by.x = "PROTEINID")
qw <- p[p$protein1%in%az$PROTEINID,]
qw <- qw[qw$protein2%in%az$PROTEINID,]
count1 = 0
count2 = 0
for (Protein1 in levels(as.factor(qw$protein1))){
        count1 = count1 + 1
        qw[qw$protein1 == Protein1,]$protein1 <- az[az$PROTEINID == Protein1,]$GENEID}
for (Protein2 in levels(as.factor(qw$protein2))){
        count2 = count2 + 1
        qw[qw$protein2 == Protein2,]$protein2 <- az[az$PROTEINID == Protein2,]$GENEID}

int_df <- as.data.frame(qw)
int_df <- int_df[c("protein1","protein2")]
print(dim(int_df))
int_df <- int_df[int_df$protein1 %in% rownames(expr0),]
int_df <- int_df[int_df$protein2 %in% rownames(expr0),]
int_df <- int_df[!duplicated(int_df),]
print(dim(int_df))
write.table(int_df, file = "2_interacciones_sobre_700_gene_symbol_in_GTEx.tab", row.names = F, quote = F, col.names = T, sep = "\t") 
#system(paste("scp -P 1313 ", "2_interacciones_sobre_700_gene_symbol_in_GTEx.tab ",Taruca_adolfo_tesis, "2_interacciones_sobre_700_gene_symbol_in_GTEx.tab", sep = ""))
## Cargar libreria de interacciones
library(dorothea) 
## Cargar archivo con elementos que se expresan en tejido mamario normalmente
matrix<- read.delim("2_mRNAs_gene-level_TMM_GTEx_matrix_with_geneID.tab")

ensembl.genes <- rownames(matrix)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
for (duplicado in levels(as.factor(geneIDs1[duplicated(geneIDs1$SYMBOL),]$SYMBOL))){
        count = 0
        for (elementos in rownames(geneIDs1[geneIDs1$SYMBOL == duplicado,])){
                count = count + 1
                geneIDs1$SYMBOL[as.numeric(elementos)] <- paste(duplicado, count, sep = "_")}}

matrix <- matrix[rownames(matrix) %in% geneIDs1$GENEID,]

for (ensembl_gene_id in rownames(matrix)){
        rownames(matrix)[rownames(matrix)== ensembl_gene_id] <- geneIDs1[geneIDs1$GENEID == ensembl_gene_id,]$SYMBOL}

reg_info <-  dorothea_hs
tf_in_breast <- reg_info
tf_in_breast <- tf_in_breast[tf_in_breast$target %in% rownames(matrix),]
tf_in_breast <- tf_in_breast[tf_in_breast$tf %in% rownames(matrix),]
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
length(unique(tf_in_breast$tf))
dim(tf_in_breast)
gencode <- read.delim("../gencode_genes.v38.annotation.tab", sep="\t", header=T)
gencode <- gencode[c("gene_name", "gene_id")]

tf_in_breast <- merge(tf_in_breast,gencode, by.x ="tf",by.y="gene_name")
length(unique(tf_in_breast$tf))
tf_in_breast$tf <- NULL
colnames(tf_in_breast)[colnames(tf_in_breast)== "gene_id"] <- "tf"
tf_in_breast <- merge(tf_in_breast,geneIDs1, by.x ="target",by.y="SYMBOL")
tf_in_breast$target <- NULL
colnames(tf_in_breast)[colnames(tf_in_breast)== "GENEID"] <- "target"
tf_in_breast <- tf_in_breast[!duplicated(tf_in_breast),]
dim(tf_in_breast)
dim(tf_in_breast[tf_in_breast$confidence == "A",])
dim(tf_in_breast[tf_in_breast$confidence == "B",])
dim(tf_in_breast[tf_in_breast$confidence == "C",])
dim(tf_in_breast[tf_in_breast$confidence == "D",])
dim(tf_in_breast[tf_in_breast$confidence == "E",])
write.table(tf_in_breast, sep = "\t", file ="5_Interacciones_TF-Target_tejido_mamario.tsv", row.names = F, quote = F, col.names = T)
#system(paste("scp -P 1313 ", "5_Interacciones_TF-Target_tejido_mamario.tsv ",Taruca_adolfo_tesis, "5_Interacciones_TF-Target_tejido_mamario.tsv", sep = ""))
print("finalizado")
datos.loc[((datos[1].isin(GTEx.index))   | (datos[1].str.contains("hsa"))) & (datos[0].isin(GTEx.index))]