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
NCORES <- 60
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
#################################################################################
message_to_bot <- 'Script Iniciado:\n"P-Value Threshold y Anotación VEP"'
bot$sendMessage(chat_id, text = message_to_bot)
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
message_to_bot <- 'Informacion cargada iniciando Loop'
bot$sendMessage(chat_id, text = message_to_bot)
for (iter in c(100)) {
p_value <- 0.01*iter
message_to_bot <- paste("P-Value < ", p_value, ":", sep = "")
bot$sendMessage(chat_id, text = message_to_bot)
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
message_to_bot <- paste("AUCinf =", AUC_inf[1], "SD =", AUC_inf[4], sep = " ")
bot$sendMessage(chat_id, text = message_to_bot)
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
ggsave(plot = Grid_plot, filename = paste("Grid_model_plot_p-value_", p_value,".png", sep = ""), height = 9, width = 16, dpi = 600)

pred_nosp <- big_prodVec(G, best_grid_nosp, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
AUC_pred_nosp <- AUCBoot(pred_nosp, y[ind.test], seed = 1)
message_to_bot <- paste("AUCgrid-nosp =", AUC_pred_nosp[1], "SD =", AUC_pred_nosp[4], sep = " ")
bot$sendMessage(chat_id, text = message_to_bot)
best_grid_sp <- params %>%
  mutate(id = row_number()) %>%
  filter(sparse) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(id) %>%
  beta_grid[, .]

pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
AUC_pred_sp <- AUCBoot(pred_sp, y[ind.test], seed = 1)
message_to_bot <- paste("AUCgrid-sp =", AUC_pred_sp[1], "SD =", AUC_pred_sp[4], sep = " ")
bot$sendMessage(chat_id, text = message_to_bot)
################################################################################################
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES), ncores = NCORES)
save(multi_auto,file = "multi_auto.respaldo.Rdata")
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
message_to_bot <- paste("AUCauto =", AUC_Auto[1], "SD =", AUC_Auto[4], sep = " ")
bot$sendMessage(chat_id, text = message_to_bot)
ggsave(plot = Auto_plot, filename = paste("Auto_model_plot_p-value_", p_value,".png", sep = ""), height = 9, width = 16, dpi = 600)
################################################################################################
info_snp$best_grid_sp <- best_grid_sp
info_snp$best_grid_nosp <- best_grid_nosp
info_snp$beta_inf <- beta_inf
info_snp$final_beta_auto <- final_beta_auto
################################################################################################
system("rm tmp-data/*.sbk")
non_zero_grid_sp <- as.data.frame(info_snp[info_snp$best_grid_sp != 0,]$rsid)
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
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = "")," adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = ""), sep = ""))
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_summary.html", sep = "")," adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_summary.html", sep = ""), sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_summary.html", sep = ""))
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = "")," adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = ""), sep = ""))
system(paste("scp ", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect", sep = "")," adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/", paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect", sep = ""), sep = ""))
bot$sendDocument(chat_id, document = paste("Grid_model_plot_p-value_", p_value,".png", sep = ""))
bot$sendDocument(chat_id, document = paste("Auto_model_plot_p-value_", p_value,".png", sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero_summary.html", sep = ""))
bot$sendDocument(chat_id, document = paste("VEP_p-Value_threshold_", p_value, "_hapmap3_all_variant_effect_non_zero", sep = ""))
prediccion_SP <- as.data.frame(cbind(obj.bigSNP$fam$sample.ID[ind.test],y[ind.test],pred_sp))
prediccion_SP$V2 <- as.character(prediccion_SP$V2)
prediccion_SP[prediccion_SP$V2 == 0,]$V2 <- "Control"
prediccion_SP[prediccion_SP$V2 == 1,]$V2 <- "Case"
prediccion_SP$V2 <- as.factor(prediccion_SP$V2)
prediccion_SP$pred_sp <- as.character(prediccion_SP$pred_sp)
prediccion_SP$pred_sp <- as.numeric(prediccion_SP$pred_sp)
less_risk <- head(prediccion_SP[order(prediccion_SP$pred_sp),],length(ind.test)*0.95)
low_risk <- head(prediccion_SP[order(prediccion_SP$pred_sp),],length(ind.test)*0.05)
odds_less_risk <- length(less_risk[less_risk$V2 == "Case",]$V2)/length(less_risk[less_risk$V2 == "Control",]$V2)
greater_risk <- tail(prediccion_SP[order(prediccion_SP$pred_sp),],length(ind.test)*0.05)
odds_greater_risk <- length(greater_risk[greater_risk$V2 == "Case",]$V2)/length(greater_risk[greater_risk$V2 == "Control",]$V2)
odds_ratio <- odds_greater_risk/odds_less_risk
ggdensity(prediccion_SP, x = "pred_sp", add = "mean", rug = TRUE, color = "V2", fill = "V2", palette = c("#00AFBB", "#E7B800")) +
    ggtitle("Distribution of test samples PRS", subtitle = paste("AUC Grid-SP = ", round(AUC_pred_sp[[1]],4), sep = ""))+ 
    theme(legend.position = c(0.85,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), legend.title = element_text(size = 15), legend.text = element_text(size = 12), legend.direction = "horizontal") +
    ylab("Density")+ xlab("PRS") + labs(fill = "Sample type", color ="Sample type")+ geom_vline(xintercept = greater_risk$pred_sp[1], na.rm = FALSE, show.legend = NA)

ggsave(paste("Density_plot_grid_SP_p-Value_threshold_", p_value,".png", sep = ""), width = 16, height = 9, dpi = 600, units = "in")

ggdotchart(prediccion_SP, x = "V1", y = "pred_sp", color = "V2", dot.size = 0.5, palette = c("#00AFBB", "#E7B800"), sorting = "ascending", ggtheme = theme_pubr()) + theme(axis.text.x=element_blank()) + 
           ylab("PRS")+ xlab("Test Samples")+ geom_vline(xintercept = greater_risk$V1[1], na.rm = FALSE, show.legend = NA)

ggsave(paste("Scores_vs_samples_p-Value_", p_value,".png", sep = ""), width = 16, height = 9, dpi = 600, units = "in")
save.image(file = paste("p-value_", p_value, ".RData"),sep = "")
}
message_to_bot <- "Analisis Finalizado"
bot$sendMessage(chat_id, text = message_to_bot)

##################################################################################################################################
less_risk <- head(prediccion_SP[order(prediccion_SP$pred_sp),],length(ind.test)*0.95)
low_risk <- head(prediccion_SP[order(prediccion_SP$pred_sp),],length(ind.test)*0.05)
odds_less_risk <- length(less_risk[less_risk$V2 == "Case",]$V2)/length(less_risk[less_risk$V2 == "Control",]$V2)
greater_risk <- tail(prediccion_SP[order(prediccion_SP$pred_sp),],length(ind.test)*0.05)
odds_greater_risk <- length(greater_risk[greater_risk$V2 == "Case",]$V2)/length(greater_risk[greater_risk$V2 == "Control",]$V2)
odds_ratio <- odds_greater_risk/odds_less_risk

prediccion_SP$decile <- ntile(prediccion_SP$pred_sp, 10)

  ggdotplot(prediccion_SP, x = "decile", y = "pred_sp",
            #color = "V2", 
            #fill = "V2",
            #palette =c("#00AFBB", "#E7B800"),
            #add = "jitter", 
            #shape = "V2",
            size = 0.1,
            #merge = TRUE
            ) +
    ggtitle("prueba")+
    theme(legend.position = "none")+
    #stat_pvalue_manual(stat.test, label = "{p.signif}\np ={p} {method} test", y.position = max(mdp$allgenes.Score)+0.5, vjust = 0.8)+
    ylab("All genes score")+ xlab("Sample type")

ggdotchart(prediccion_SP, x = "decile", y = "pred_sp", color = "V2", dot.size = 0.5, palette = c("#00AFBB", "#E7B800"), sorting = "ascending", ggtheme = theme_pubr()) + theme(axis.text.x=element_blank()) + 
           ylab("PRS")+ xlab("Test Samples")+ geom_vline(xintercept = greater_risk$V1[1], na.rm = FALSE, show.legend = NA)


ggsave("prueba.png", width = 16, height = 9, dpi = 600, units = "in")

