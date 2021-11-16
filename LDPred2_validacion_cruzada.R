#!/usr/bin/env Rscript
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(telegram.bot)
bot = Bot(token = bot_token("AARH_95_bot"))
chat_id <- bot_token("chat_id")
library(ggpubr)
Taruca_adolfo_tesis <- "adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_1/"
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
## **Matching variants between genotype data and summary statistics**
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

best_grid_sp <- params %>%
  mutate(id = row_number()) %>%
  filter(sparse) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(id) %>%
  beta_grid[, .]

pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
AUC_pred_sp <- AUCBoot(pred_sp, y[ind.test], seed = 1)
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
system("scp Efectos_beta_validacion_cruzada.tsv adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_1/Efectos_beta_validacion_cruzada.tsv")

info_snp$best_grid_nosp_mean_K <- rowMeans(K_betas[grepl("best_grid_nosp_K",names(K_betas))])
info_snp$best_grid_sp_mean_K <- rowMeans(K_betas[grepl("best_grid_sp_K",names(K_betas))])
info_snp$beta_inf_mean_K <- rowMeans(K_betas[grepl("beta_inf_K",names(K_betas))])
info_snp$final_beta_auto_mean_K <- rowMeans(K_betas[grepl("final_beta_auto_K",names(K_betas))])

pred_inf <- big_prodVec(G, info_snp$beta_inf_mean_K, ind.row = ind.final.val, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
AUC_inf <- AUCBoot(pred_inf, y[ind.final.val], seed = 1)

pred_nosp <- big_prodVec(G, info_snp$best_grid_nosp_mean_K, ind.row = ind.final.val, ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)
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

full_list2 <- full_list[1:(nrow(full_list)-1),]
full_list2$Infinitesimal <- sapply(strsplit(full_list2$Infinitesimal," \\("), `[`, 1)
full_list2$'Grid-NOSP' <- sapply(strsplit(full_list2$'Grid-NOSP'," \\("), `[`, 1)
full_list2$'Grid-SP' <- sapply(strsplit(full_list2$'Grid-SP'," \\("), `[`, 1)
full_list2$Automatico <- sapply(strsplit(full_list2$Automatico," \\("), `[`, 1)
full_list2$Infinitesimal <- as.numeric(full_list2$Infinitesimal)
full_list2$'Grid-NOSP' <- as.numeric(full_list2$'Grid-NOSP')
full_list2$'Grid-SP' <- as.numeric(full_list2$'Grid-SP')
full_list2$Automatico <- as.numeric(full_list2$Automatico)
crossval_rendimiento <- colMeans(full_list2[,2:5])
crossval_rendimiento <- data.frame("AUC rendimiento promedio",crossval_rendimiento[[1]],crossval_rendimiento[[2]],crossval_rendimiento[[3]],crossval_rendimiento[[4]],fix.empty.names = F)
colnames(crossval_rendimiento) <- c("K","Infinitesimal", "Grid-NOSP", "Grid-SP", "Automatico")

full_list2 <- rbind(full_list2, crossval_rendimiento)
rendimiento_258 <- data.frame("AUC estimado promedio betas",
                    paste(round(AUC_inf[1],6)," (", round(AUC_inf[4],4),")", sep = ""),
                    paste(round(AUC_pred_nosp[1],6)," (", round(AUC_pred_nosp[4],4),")", sep = ""),
                    paste(round(AUC_pred_sp[1],6)," (", round(AUC_pred_sp[4],4),")", sep = ""), 
                    paste(round(AUC_Auto[1],6)," (", round(AUC_Auto[4],4),")", sep = ""),
                    fix.empty.names = F)
colnames(rendimiento_258) <- c("K","Infinitesimal", "Grid-NOSP", "Grid-SP", "Automatico")
full_list2 <- rbind(full_list2, rendimiento_258)
variantes_geneticas <- info_snp$rsid
write.table(variantes_geneticas, sep = "\t", file ="ids_rs.txt", row.names = F, quote = F, col.names = F)

write.table(full_list2, sep = "\t", file ="Rendimientos_validacion_cruzada.tsv", row.names = F, quote = F, col.names = T)
bot$sendDocument(chat_id, document = "Rendimientos_validacion_cruzada.tsv")
save.image(file = "2_entorno_validacion_cruzada.RData")
###################################################################################################################################
# Analisis resultados promedio
## eleccion mejor modelo en base a AUC mayor
if (AUC_pred_sp[[1]] > AUC_inf[[1]] & AUC_pred_sp[[1]] > AUC_pred_nosp[[1]] & AUC_pred_sp[[1]] > AUC_Auto[[1]]) {
   best_auc_name <- "Grid SP model"
   best_auc_name2 <- "Grid_SP_model"
   best_auc_value <- AUC_pred_sp[[1]]
   best_auc_pred <- pred_sp
} else if (AUC_pred_nosp[[1]] > AUC_inf[[1]] & AUC_pred_nosp[[1]] > AUC_pred_sp[[1]] & AUC_pred_nosp[[1]] > AUC_Auto[[1]]) {
   best_auc_name <- "Grid NOSP model"
   best_auc_name2 <- "Grid_SP_model"
   best_auc_value <- AUC_pred_nosp[[1]]
   best_auc_pred <- pred_nosp
} else if (AUC_Auto[[1]] > AUC_inf[[1]] & AUC_Auto[[1]] > AUC_pred_nosp[[1]] & AUC_Auto[[1]] > AUC_pred_sp[[1]]) {
   best_auc_name <- "Automatic model"
   best_auc_name2 <- "Automatic_model"
   best_auc_value <- AUC_Auto[[1]]
   best_auc_pred <- final_pred_auto
} else if (AUC_inf[[1]] > AUC_pred_sp[[1]] & AUC_inf[[1]] > AUC_pred_nosp[[1]] & AUC_inf[[1]] > AUC_Auto[[1]]) {
   best_auc_name <- "Infinitesimal model"
   best_auc_name2 <- "Infinitesimal_model"
   best_auc_value <- AUC_inf[[1]]
   best_auc_pred <- pred_inf
}  
## VEP Analisis
non_zero_grid_sp <- as.data.frame(info_snp[info_snp$best_grid_sp_mean_K != 0,]$rsid)
names(non_zero_grid_sp) <- "rsid"
df2 <- info_snp
DF3 <- df2[c("chr","pos","a1","a0","rsid")]
DF3$alelos <- paste(DF3$a1, DF3$a0, sep = "/")
DF3$pos2 <- DF3$pos
DF3$strand <- NA
DF3 <- DF3[c("chr","pos","pos2","alelos","strand","rsid")]
DF3 <- DF3[order(DF3$chr,DF3$pos),]
DF2 <- merge(DF3, non_zero_grid_sp, by = "rsid")
DF2 <- DF2[c("chr","pos","pos2","alelos","strand","rsid")]
DF2 <- DF2[order(DF2$chr,DF2$pos),]
write.table(DF3,file = paste("p-Value_threshold_", p_value, "_hapmap3_all_variant_effect.txt", sep = ""), quote = F, col.names=F, row.names= F)
write.table(DF2,file = paste("p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero.txt", sep = ""), quote = F, col.names=F, row.names= F)
system(paste("perl /media/storage2/software/ensembl-vep/vep --input_file ", paste("p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero.txt", sep = ""), " --output_file ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = ""), " --cache --regulatory --no_headers --biotype --dir_cache /media/storage2/software/ensembl-vep/cache/ --force_overwrite --verbose --per_gene --symbol --domains --fork 15 --offline --assembly GRCh37", sep = ""))
system(paste("perl /media/storage2/software/ensembl-vep/vep --input_file ", paste("p-Value_threshold_", p_value, "_hapmap3_all_variant_effect.txt", sep = ""), " --output_file ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect", sep = ""), " --cache --regulatory --no_headers --biotype --dir_cache /media/storage2/software/ensembl-vep/cache/ --force_overwrite --verbose --per_gene --symbol --domains --fork 15 --offline --assembly GRCh37", sep = ""))
system("rm *_warnings.txt")
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = "")," ", Taruca_adolfo_tesis, paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = ""), sep = ""))
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_summary.html", sep = "")," ", Taruca_adolfo_tesis, paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_summary.html", sep = ""), sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_summary.html", sep = ""))
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = "")," ", Taruca_adolfo_tesis, paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = ""), sep = ""))
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect", sep = "")," ", Taruca_adolfo_tesis, paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect", sep = ""), sep = ""))
#bot$sendDocument(chat_id, document = paste("Grid_model_plot_p-value_", p_value,".png", sep = ""))
#bot$sendDocument(chat_id, document = paste("Auto_model_plot_p-value_", p_value,".png", sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = ""))
# Odds Ratio Analisis
Prediccion <- as.data.frame(cbind(obj.bigSNP$fam$sample.ID[ind.final.val],y[ind.final.val],best_auc_pred))
Prediccion$V2 <- as.character(Prediccion$V2)
Prediccion[Prediccion$V2 == 0,]$V2 <- "Control"
Prediccion[Prediccion$V2 == 1,]$V2 <- "Case"
Prediccion$V2 <- as.factor(Prediccion$V2)
Prediccion$best_auc_pred <- as.character(Prediccion$best_auc_pred)
Prediccion$best_auc_pred <- as.numeric(Prediccion$best_auc_pred)
less_risk <- head(Prediccion[order(Prediccion$best_auc_pred),],length(ind.test)*0.95)
low_risk <- head(Prediccion[order(Prediccion$best_auc_pred),],length(ind.test)*0.05)
odds_less_risk <- length(less_risk[less_risk$V2 == "Case",]$V2)/length(less_risk[less_risk$V2 == "Control",]$V2)
greater_risk <- tail(Prediccion[order(Prediccion$best_auc_pred),],length(ind.test)*0.05)
odds_greater_risk <- length(greater_risk[greater_risk$V2 == "Case",]$V2)/length(greater_risk[greater_risk$V2 == "Control",]$V2)
odds_ratio <- odds_greater_risk/odds_less_risk
ggdensity(Prediccion, x = "best_auc_pred", add = "mean", rug = TRUE, color = "V2", fill = "V2", palette = c("#00AFBB", "#E7B800")) +
    ggtitle("Distribution of test samples PRS", subtitle = paste(best_auc_name, " = ", round(best_auc_value,4), sep = ""))+ 
    theme(legend.position = c(0.85,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), legend.title = element_text(size = 15), legend.text = element_text(size = 12), legend.direction = "horizontal") +
    ylab("Density")+ xlab("PRS") + labs(fill = "Sample type", color ="Sample type")+ geom_vline(xintercept = greater_risk$best_auc_pred[1], na.rm = FALSE, show.legend = NA)

ggsave(paste("Density_plot_",best_auc_name2,"_p-Value_threshold_", p_value,".png", sep = ""), width = 16, height = 9, dpi = 600, units = "in")
system(paste("scp ", paste("Density_plot_",best_auc_name2,"_p-Value_threshold_", p_value,".png", sep = "")," ", Taruca_adolfo_tesis, paste("Density_plot_",best_auc_name2,"_p-Value_threshold_", p_value,".png", sep = ""), sep = ""))

ggdotchart(Prediccion, x = "V1", y = "best_auc_pred", color = "V2", dot.size = 0.5, palette = c("#00AFBB", "#E7B800"), sorting = "ascending", ggtheme = theme_pubr()) + theme(axis.text.x=element_blank()) + 
           ylab("PRS")+ xlab("Test Samples")+ geom_vline(xintercept = greater_risk$V1[1], na.rm = FALSE, show.legend = NA)

ggsave(paste("Scores_vs_samples_p-Value_", p_value,".png", sep = ""), width = 16, height = 9, dpi = 600, units = "in")
system(paste("scp ", paste("Scores_vs_samples_p-Value_", p_value,".png", sep = "")," ", Taruca_adolfo_tesis, paste("Scores_vs_samples_p-Value_", p_value,".png", sep = ""), sep = ""))

#save.image(file = paste("p-value_", p_value, ".RData"),sep = "")

###################################################################################################################################
message_to_bot <- "Analisis Finalizado"
bot$sendMessage(chat_id, text = message_to_bot)
##################################################################################################################################