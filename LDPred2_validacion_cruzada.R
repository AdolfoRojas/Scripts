#!/usr/bin/env Rscript
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(telegram.bot)
bot = Bot(token = bot_token("AARH_95_bot"))
chat_id <- bot_token("chat_id")
library(ggpubr)
#################################################################################################
system("rm *.rds")
system("rm *.bk")
snp_readBed("EUR_CGEMS.QC.bed") # Read from bed/bim/fam, it generates .bk and .rds files.
obj.bigSNP <- snp_attach("EUR_CGEMS.QC.rds") # Attach the "bigSNP" object in R session
# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- 40
sumstats <- bigreadr::fread2("GWAS_PanUKBB_final.tsv") # Read external summary statistics

set.seed(1)
ind.final.val <- sample(nrow(G), 258)
poblacion_restante <- setdiff(rows_along(G), ind.final.val)
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
#################################################################################
message_to_bot <- 'Script Iniciado:\n"Validacion cruzada LDPred2"'
bot$sendMessage(chat_id, text = message_to_bot)
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
message_to_bot <- 'Informacion cargada iniciando Loop'
bot$sendMessage(chat_id, text = message_to_bot)

p_value <- 1
message_to_bot <- paste("P-Value < ", p_value, ":", sep = "")
bot$sendMessage(chat_id, text = message_to_bot)

full_list <- data.frame(Iteraccion=character(),                
                Infinitesimal=character(), 
                Grid_NOSP=character(), 
                Grid_SP=character(), 
                Automatico=character(),                 
                stringsAsFactors=FALSE)
for(semilla in 1:10){
message_to_bot <- paste('K = ', semilla, sep = " ")
bot$sendMessage(chat_id, text = message_to_bot)
set.seed(semilla)
#ind.val <- sample(poblacion_restante, 1500)
#ind.test <- setdiff(poblacion_restante, ind.val)
ind.test <- poblacion_restante[(1+200*(semilla-1)):(200+200*(semilla-1))]
ind.val <- setdiff(poblacion_restante, ind.test)

info_snp <- snp_match(sumstats[sumstats$p < p_value,], map, strand_flip = FALSE)
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)
tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
for (chr in 1:22) {    
    ind.chr <- which(info_snp$chr == chr) ## indices in 'info_snp'
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr] ## indices in 'G'
    corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)
    if (chr == 1) {
      df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
      corr <- as_SFBM(corr0, tmp)
      ld <- Matrix::colSums(corr0^2)
    } else {
      df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
      corr$add_columns(corr0, nrow(corr))
      ld <- c(ld, Matrix::colSums(corr0^2))
}}
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
AUC_inf <- AUCBoot(pred_inf, y[ind.test], seed = 1)
#message_to_bot <- paste("AUCinf =", AUC_inf[1], "SD =", AUC_inf[4], sep = " ")
#bot$sendMessage(chat_id, text = message_to_bot)
################################################################################################
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
params$score <- big_univLogReg(as_FBM(pred_grid[ind.val, ]), y[ind.val])$score

Grid_plot <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
 
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

pred_nosp <- big_prodVec(G, best_grid_nosp, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
AUC_pred_nosp <- AUCBoot(pred_nosp, y[ind.test], seed = 1)
#message_to_bot <- paste("AUCgrid-nosp =", AUC_pred_nosp[1], "SD =", AUC_pred_nosp[4], sep = " ")
#bot$sendMessage(chat_id, text = message_to_bot)
best_grid_sp <- params %>%
  mutate(id = row_number()) %>%
  filter(sparse) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(id) %>%
  beta_grid[, .]

pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
AUC_pred_sp <- AUCBoot(pred_sp, y[ind.test], seed = 1)
#message_to_bot <- paste("AUCgrid-sp =", AUC_pred_sp[1], "SD =", AUC_pred_sp[4], sep = " ")
#bot$sendMessage(chat_id, text = message_to_bot)
################################################################################################
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES), ncores = NCORES)
auto <- multi_auto[[1]]
Auto_plot <- plot_grid(
  qplot(y = auto$path_p_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv")
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])
final_pred_auto <- rowMeans(pred_auto[, keep])
AUC_Auto <- AUCBoot(final_pred_auto, y[ind.test], seed = 1)
#message_to_bot <- paste("AUCauto =", AUC_Auto[1], "SD =", AUC_Auto[4], sep = " ")
#bot$sendMessage(chat_id, text = message_to_bot)
nam_best_grid_sp <- paste("best_grid_sp_K", semilla, sep = "")
assign(nam_best_grid_sp, best_grid_sp)

nam_best_grid_nosp <- paste("best_grid_nosp_K", semilla, sep = "")
assign(nam_best_grid_nosp, best_grid_nosp)

nam_beta_inf <- paste("beta_inf_K", semilla, sep = "")
assign(nam_beta_inf, beta_inf)

nam_final_beta_auto <- paste("final_beta_auto_K", semilla, sep = "")
assign(nam_final_beta_auto, final_beta_auto)
if (semilla == 1){
  K_betas <- as.data.frame(cbind(get(nam_beta_inf), get(nam_best_grid_nosp), get(nam_best_grid_sp), get(nam_final_beta_auto)))
  names(K_betas) <- c(nam_beta_inf,
                      nam_best_grid_nosp,
                      nam_best_grid_sp,
                      nam_final_beta_auto)
} else {
   K_betas <- as.data.frame(cbind(K_betas, get(nam_beta_inf), get(nam_best_grid_nosp), get(nam_best_grid_sp), get(nam_final_beta_auto)))
   colnames(K_betas)[ncol(K_betas)-3] <- nam_beta_inf
   colnames(K_betas)[ncol(K_betas)-2] <- nam_best_grid_nosp
   colnames(K_betas)[ncol(K_betas)-1] <- nam_best_grid_sp
   colnames(K_betas)[ncol(K_betas)] <- nam_final_beta_auto
}
full_list <- rbind(full_list, data.frame(semilla,
                    paste(round(AUC_inf[1],6)," (", round(AUC_inf[4],4),")", sep = ""),
                    paste(round(AUC_pred_nosp[1],6)," (", round(AUC_pred_nosp[4],4),")", sep = ""),
                    paste(round(AUC_pred_sp[1],6)," (", round(AUC_pred_sp[4],4),")", sep = ""), 
                    paste(round(AUC_Auto[1],6)," (", round(AUC_Auto[4],4),")", sep = ""),
                    fix.empty.names = F))
system("rm tmp-data/*.sbk")
}
write.table(K_betas, sep = "\t", file ="Efectos_beta_validacion_cruzada.tsv", row.names = F, quote = F, col.names = T)
system("scp Efectos_beta_validacion_cruzada.tsv adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Efectos_beta_validacion_cruzada.tsv")

info_snp$best_grid_nosp_mean_K <- rowMeans(K_betas[grepl("best_grid_nosp_K",names(K_betas))])
info_snp$best_grid_sp_mean_K <- rowMeans(K_betas[grepl("best_grid_sp_K",names(K_betas))])
info_snp$beta_inf_mean_K <- rowMeans(K_betas[grepl("beta_inf_K",names(K_betas))])
info_snp$final_beta_auto_mean_K <- rowMeans(K_betas[grepl("final_beta_auto_K",names(K_betas))])

pred_inf <- big_prodVec(G, info_snp$beta_inf_mean_K, ind.row = ind.final.val, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
AUC_inf <- AUCBoot(pred_inf, y[ind.final.val], seed = 1)

pred_nosp <- big_prodVec(G, best_grid_nosp, ind.row = ind.final.val, ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)
AUC_pred_nosp <- AUCBoot(pred_nosp, y[ind.final.val], seed = 1)

pred_sp <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = ind.final.val, ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)
AUC_pred_sp <- AUCBoot(pred_sp, y[ind.final.val], seed = 1)

final_pred_auto <- big_prodVec(G, info_snp$final_beta_auto_mean_K , ind.row = ind.final.val, ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)
AUC_Auto <- AUCBoot(final_pred_auto, y[ind.final.val], seed = 1)

full_list <- rbind(full_list, data.frame("AUC estimado promedio betas",
                    paste(round(AUC_inf[1],6)," (", round(AUC_inf[4],4),")", sep = ""),
                    paste(round(AUC_pred_nosp[1],6)," (", round(AUC_pred_nosp[4],4),")", sep = ""),
                    paste(round(AUC_pred_sp[1],6)," (", round(AUC_pred_sp[4],4),")", sep = ""), 
                    paste(round(AUC_Auto[1],6)," (", round(AUC_Auto[4],4),")", sep = ""),
                    fix.empty.names = F))

colnames(full_list) <- c("K","Infinitesimal", "Grid-NOSP", "Grid-SP", "Automatico")
write.table(full_list, sep = "\t", file ="Rendimientos_validacion_cruzada.tsv", row.names = F, quote = F, col.names = T)
bot$sendDocument(chat_id, document = "Rendimientos_validacion_cruzada.tsv")
save.image(file = "2_entorno_validacion_cruzada.RData")
message_to_bot <- "Analisis Finalizado"
bot$sendMessage(chat_id, text = message_to_bot)

##################################################################################################################################