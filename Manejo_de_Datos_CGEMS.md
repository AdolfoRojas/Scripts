## **Se uso ligate_minimac.py para unir los archivos chunk en que estaban los datos imputados** ##


## **Convertir de archivos .info y .dose de minimac a plink por cromosoma**
~~~
for i in {1..23} 
do    
	./gcta_1.93.2beta/gcta64  --thread-num 3 --dosage-mach-gz imputed_chr$i.dose.gz imputed_chr$i.info.gz --maf 0.01 --imput-rsq 0.3 --make-bed --out imputed_chr$i
done
~~~

## **Unir archivos formato plink por cromosoma a genoma completo**

~~~
for i in {2..23} 
do
	echo imputed_chr$i.bed imputed_chr$i.bim imputed_chr$i.fam >> unir_cromosomas.txt
done

plink1.9 --noweb --make-bed --bfile imputed_chr1 --merge-list unir_cromosomas.txt --out CGEMS
~~~

## **Añadir info caso-control**
~~~
grep -v "#" ../phg000032.v1.CGEMS_BreastCancer.sample-info.MULTI.txt| cut -d, -f 1,4,8 > phenotypes_CGEMS.txt

sed -i 's/,/\t/g' phenotypes_CGEMS.txt 

grep -E "2$" phenotypes_CGEMS.txt | cut -d$'\t' -f 1,2> cases_pheno_CGEMS.txt

plink1.9 --bfile CGEMS --make-pheno cases_pheno_CGEMS.txt '*' --make-bed --out CGEMS_pheno_added
~~~

## **Añadir info de sexo de individuos**
~~~
grep -v "#" ../phg000032.v1.CGEMS_BreastCancer.sample-info.MULTI.txt| cut -d, -f 1,4,7 | sort | uniq > sex_CGEMS.txt

sed -i 's/,/\t/g' sex_CGEMS.txt 

plink1.9 --bfile CGEMS_pheno_added --update-sex sex_CGEMS.txt --make-bed --out CGEMS_pheno_sex_added
~~~

## **Arreglar informacion referente a cromosomas**

~~~
R
chr <- read.delim("CGEMS_pheno_sex_added.bim", header= F)
chr2 <- chr[c("V2","V1")]
chr2$V1 <- sapply(strsplit(as.character(chr2$V2),'\\:'), "[", 1)
write.table(chr2, file = "reemplazo_cromosomas.txt", row.names= F, quote= F, col.names= F, sep = "\t")
q()
plink1.9 --bfile CGEMS_pheno_sex_added --update-chr reemplazo_cromosomas.txt 2 1 0 --make-bed --out CGEMS_pheno_sex_chr_added
~~~

## **Arreglar Variantes geneticas del tipo InDel y busqueda de genoma de referencia** 

~~~
indels <- read.delim("CGEMS_pheno_sex_chr_added.bim", header= F)
indels2 <- indels[c("V2", "V5", "V6")]
indels2 <- indels2[grepl("\\:[^0-9]", indels2$V2),]
indels2$alelos <- sapply(strsplit(as.character(indels2$V2),'\\:'), "[", 3)
indels2$Ref <- sapply(strsplit(as.character(indels2$alelos),'\\_'), "[", 1)
indels2$Alt <- sapply(strsplit(as.character(indels2$alelos),'\\_'), "[", 2)
indels2[is.na(indels2$Alt),]$Alt <- "-"
indels2$alelos <- NULL
indels2 <- indels2[c("V2", "V5", "V6", "Alt", "Ref")]
write.table(indels2, file = "reemplazo_alelos_InDels.txt", row.names= F, quote= F, col.names= F, sep = "\t")
# se determino mediante Genome Browser buscando variantes como X:99974528:ACCT_ en GRCH 36, 37 y 38, encontrando coincidencia con grch37
# preguntar el Lunes por la implicancia de que el alelo de referencia sea menor en los casos de InDels donde ocurre que en V5 esta R y en V6 la letra que designa el tipo de variante (I/D)
system("plink1.9 --bfile CGEMS_pheno_sex_chr_added --update-alleles  reemplazo_alelos_InDels.txt --make-bed --out CGEMS_pheno_sex_chr_added_fixed_alleles") 
~~~


## **Añadir info de posicion de las variantes** 

~~~
coordenadas <- read.delim("CGEMS_pheno_sex_chr_added_fixed_alleles.bim", header= F)
coordenadas2 <- coordenadas["V2"]
coordenadas2$pos <- sapply(strsplit(as.character(coordenadas2$V2),'\\:'), "[", 2)
write.table(coordenadas2, file = "añadir_posicion.txt", row.names= F, quote= F, col.names= F, sep = "\t")
system("plink1.9 --bfile CGEMS_pheno_sex_chr_added_fixed_alleles --update-map añadir_posicion.txt --make-bed --out CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles")
q()
~~~

## **Trabajo en resumen estadistico GWAS Pan cancer UK BB** 
chr pos ref alt af_cases_meta af_controls_meta beta_meta se_meta pval_meta pval_heterogeneity af_cases_AFR af_cases_CSA af_cases_EAS    af_cases_EUR af_controls_AFR af_controls_CSA af_controls_EAS af_cont
~~~
# En /home/venus/mar/colab/ahidalgo/Cancer-PRS-decrypted/Cancer-PRS/CGEMS/79292/PhenoGenotypeFiles/RootStudyConsentSet_phs000147.CGEMS_BreastCancer.v3.p1.c1.GRU/GenotypeFiles/joined_chrN_minimac/PRS_inicial
head -2000 phecode-174.1-females.tsv | cut -d$'\t' -f1-4,11,14,15,18,19,22,23,26,27,30,31,34 # AFR y EUR
cut -d$'\t' -f 1-4,14,18,22,26,30,34 phecode-174.1-females.tsv > Pan_UK_GWAS_EUR_only.tsv
head full_variant_qc_metrics.txt | cut -d$'\t' -f1-11,33-38
cut -d$'\t' -f1-11,33-38 full_variant_qc_metrics.txt > Variant_info_EUR.tsv
# N casos EUR 13257; N controles EUR 205913
R
eur <- read.delim("Pan_UK_GWAS_EUR_only_NA_rm.tsv", header= T)
eur$low_confidence_EUR <- NULL
cases <- 13257
controles <- 205913
eur$frqGWAS <- ((eur$af_controls_EUR*controles*2)+(eur$af_cases_EUR*cases*2))/(cases*2 + controles*2)
eur2 <- eur[eur$frqGWAS > 0.01,]
write.table(eur2, file = "Pan_UK_GWAS_EUR_only_NA_rm_freq_added.tsv", row.names= F, quote= F, col.names= T, sep = "\t")

variant_info <- read.delim("Variant_info_EUR.tsv", header= T)
filtered_variants <- variant_info[variant_info$info > 0.8,]
write.table(filtered_variants, file = "Variant_info_EUR_filtered.tsv", row.names= F, quote= F, col.names= T, sep = "\t")

end_file <- merge(eur2, filtered_variants, by.x= c("chr","pos","ref","alt"), by.y= c("chrom","pos","ref","alt"))
Nsamples <- cases + controles
end_file$N <- Nsamples
end_file2 <- end_file[c("chr","pos","rsid","alt","ref","N","se_EUR","pval_EUR","beta_EUR","info","frqGWAS")]
names(end_file2) <- c("CHR","BP","SNP","A1","A2","N","SE","P","BETA","INFO","MAF")
write.table(end_file2, file = "GWAS_PanUKBB_final.tsv", row.names= F, quote= F, col.names= T, sep = "\t")
~~~

## **Control de calidad Base Data** 

#### **Standard GWAS QC** 
~~~
cat GWAS_PanUKBB_final.tsv | awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}'> GWAS_PanUKBB_final ### realizado en pasos anteriores
~~~
#### **Ambiguous SNPs**
cat GWAS_PanUKBB_final.tsv | awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' > GWAS_PanUKBB_final_not_ambiguous.tsv
#### Duplicate SNPs
cat GWAS_PanUKBB_final_not_ambiguous.tsv | awk '{ print $3}' | sort | uniq -d > duplicated.snp ## wc -l duplicated.snp -> 7665 ## snp2 es del archivo no ambiguo
cat GWAS_PanUKBB_final_not_ambiguous.tsv | grep -vf duplicated.snp2 > GWAS_PanUKBB_final_not_ambiguous.nodup.tsv


## **Control de calidad Target Data** 
~~~
plink1.9 --bfile CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out EUR_CGEMS.QC
~~~
#plink1.9 --bfile CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles --keep EUR_CGEMS.QC.fam --extract EUR_CGEMS.QC.snplist --indep-pairwise 200 50 0.25 --out EUR_CGEMS.QC
~~~
plink1.9 --bfile CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles --keep EUR_CGEMS.QC.fam --extract EUR_CGEMS.QC.snplist --het --out EUR_CGEMS.QC

R
dat <- read.table("EUR_CGEMS.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "EUR_CGEMS.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R
~~~
##### **--EUR_CGEMS.valid.sample**

#### **Merge variants**
~~~
R
# Read in bim file
bim <- read.table("../CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles.bim")
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
# Read in QCed SNPs
qc <- read.table("../EUR_CGEMS.QC.snplist", header = F, stringsAsFactors = F)
# Read in the GWAS data
height <- read.table("GWAS_PanUKBB_final.tsv", header = T, stringsAsFactors = F, sep="\t")
# Change all alleles to upper case for easy comparison
height$A1 <- toupper(height$A1)
height$A2 <- toupper(height$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)
#height[height$CHR == "X",]$CHR <- "23" ### mas adelante se analizan por cromosomas y el 23 (X) no es considerado
# Merge summary statistic with target
info <- merge(bim, height, by.x = c("CHR", "BP", "B.A1", "B.A2"), by.y = c("CHR", "BP", "A1", "A2"))
info <- info[info$SNP.x %in% qc$V1,]# Filter QCed SNPs
match_variants <- bim$SNP[bim$SNP %in% info$SNP.x]
write.table(match_variants, "EUR_CGEMS_GWAS.shared_variants", quote = F, row.names = F, col.names = F)
q() # exit R
~~~
##### **--EUR_CGEMS_GWAS.shared_variants**

~~~
##plink1.9 --bfile CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles --keep EUR_CGEMS.valid.sample --extract PRS_inicial/EUR_CGEMS_GWAS.shared_variants --check-sex --out EUR_CGEMS.QC

# Read in file
#valid <- read.table("EUR_CGEMS.valid.sample", header=T)
#dat <- read.table("EUR_CGEMS.QC.sexcheck", header=T)
#valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
#write.table(valid[,c("FID", "IID")], "EUR_CGEMS.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
#q() # exit R

plink1.9 --bfile CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles --extract PRS_inicial/EUR_CGEMS_GWAS.shared_variants --keep EUR_CGEMS.valid.sample --rel-cutoff 0.125 --out EUR_CGEMS.QC ## EUR_CGEMS.QC.valid reemplazado por --keep EUR_CGEMS.valid.sample  por haberse eliminado variantes cromosoma X
~~~


~~~
plink1.9 --bfile CGEMS_pheno_sex_chr_pos_added_fixed_indels_alleles --make-bed --keep EUR_CGEMS.QC.rel.id --out EUR_CGEMS.QC --extract PRS_inicial/EUR_CGEMS_GWAS.shared_variants
~~~

# **LDPred**

~~~
R

library(bigsnpr)
#Sys.setenv(LANG="en_US.UTF-8")
system("rm *.QC.rds")
system("rm *.QC.bk")
snp_readBed("EUR_CGEMS.QC.bed") # Read from bed/bim/fam, it generates .bk and .rds files.
obj.bigSNP <- snp_attach("EUR_CGEMS.QC.rds") # Attach the "bigSNP" object in R session
# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- 60 # nb_cores()
sumstats <- bigreadr::fread2("GWAS_PanUKBB_final.tsv") # Read external summary statistics

set.seed(1)
ind.val <- sample(nrow(G), 1500)
ind.test <- setdiff(rows_along(G), ind.val)
#~~~
## **Matching variants between genotype data and summary statistics**

#~~~
sumstats$n_case <- 13257
sumstats$n_control <- 205913
sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
sumstats$n_case <- sumstats$n_control <- NULL
names(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "N", "beta_se", "p", "beta", "INFO", "MAF", "n_eff")
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
sumstats[sumstats$chr == "X",]$chr <- "23"
sumstats$chr <- as.integer(sumstats$chr)
info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
save.image(file='myEnvironment_before_loop.RData')
~~~

~~~
load('myEnvironment_before_loop.RData')
library(bigsnpr)
Sys.setenv(LANG="en_US.UTF-8")
info_snp <- snp_match(sumstats[sumstats$p < 0.5,], map, strand_flip = FALSE)

for (chr in 1:22) {

    print(chr)

    ## indices in 'info_snp'
    ind.chr <- which(info_snp$chr == chr)
    ## indices in 'G'
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    
    corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)

    if (chr == 1) {
      df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp)
    } else {
      df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
}
}
save.image(file='corr_done.RData')
~~~

~~~
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]
~~~


# Infinitesimal model
~~~
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
AUCBoot(pred_inf, y[ind.test])
##########################################################################################
Modificar Beta y numero de variantes a trabajar 
##########################################################################################

betas_index_0.001 <- which(beta_inf > 0.001| beta_inf < -0.001)

~~~

# Grid of models
~~~
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = info_snp$`_NUM_ID_`)####################################
params$score <- big_univLogReg(as_FBM(pred_grid[ind.val, ]), y[ind.val])$score
library(ggplot2)
ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
  library(dplyr)
  params %>%
  mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
  arrange(desc(score)) %>%
  mutate_at(c("score", "sparsity"), round, digits = 3) %>%
  slice(1:10)
  best_grid_nosp <- params %>%
  mutate(id = row_number()) %>%
  filter(!sparse) %>% 
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(id) %>% 
  beta_grid[, .]

pred_nosp <- big_prodVec(G, best_grid_nosp, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)############################
AUCBoot(pred_nosp, y[ind.test])
best_grid_sp <- params %>%
  mutate(id = row_number()) %>%
  filter(sparse) %>% 
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(id) %>% 
  beta_grid[, .]

pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)################################
AUCBoot(pred_sp, y[ind.test])
~~~

# Automatic model
~~~
# takes a few minutes
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
                               ncores = NCORES)
str(multi_auto)
auto <- multi_auto[[1]]
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)#######################################
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])
final_pred_auto <- rowMeans(pred_auto[, keep])
AUCBoot(final_pred_auto, y[ind.test])
~~~